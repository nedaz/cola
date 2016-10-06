#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "base/StringUtil.h"
#include "extern/logger/log.h"
#include "base/RandomStuff.h"
#include "src/FastAlign/SeedingThreads.h"
#include "src/FastAlign/AlignmentThreads.h"
#include "src/FastAlign/FastAlignUnit.h"

//======================================================

void FastAlignUnit::findSeeds(int targetSeqIdx) {
    m_queryUnit.findSeeds(m_targetSeqs[targetSeqIdx], m_params.getSeedSize(), m_seeds[targetSeqIdx]);
}

void FastAlignUnit::findSyntenicBlocks(int targetSeqIdx, svec<SyntenicSeeds>& maxSynts) const {
    const SeedArray& seeds = m_seeds[targetSeqIdx];
    int currQIdx     = -1;
    int prevIdx      =  0;  // Keep track of the previous index where the query index changed
    for(int seedIdx=0; seedIdx<seeds.getNumSeeds(); seedIdx++) {
        int qIdx = seeds[seedIdx].getQueryIdx();
        if((currQIdx != qIdx && currQIdx != -1) || seedIdx==seeds.getNumSeeds()-1) { //New query sequence 
            // Find the best synteny & save if seed coverage of sequence passes acceptance threshold
            const SyntenicSeeds& ss = searchDPSynteny(seeds, prevIdx, seedIdx);
            if(ss.getSeedCoverage(m_params.getSeedSize()*2) > m_params.getMinSeedCover()){
                maxSynts.push_back(ss);
            }
            prevIdx = seedIdx;
        }
        currQIdx = qIdx; 
    }
} 

SyntenicSeeds FastAlignUnit::searchDPSynteny(const SeedArray& seeds, int startQIdx, int endQIdx) const {
    SyntenicSeedFinder sFinder(SeedsSubset(seeds, startQIdx, endQIdx)); 
    return sFinder.searchDP();
}

/*
SyntenicSeeds FastAlignUnit::searchDPSynteny(const SeedArray& seeds, int startQIdx, int endQIdx) const {
    SyntenicSeeds synteny(seeds[startQIdx]); // Initialize new synteny with the latest seed and add to existing candidates
    for(int seedIdx=startQIdx; seedIdx<endQIdx; seedIdx++) { 
        synteny.addSeed(seeds[seedIdx]); 
    }
    return synteny;
}
*/
/* 
SyntenicSeeds FastAlignUnit::searchDPSynteny(const SeedArray& seeds, int startQIdx, int endQIdx) const {
    svec<SyntenicSeeds> syntenies;
    for(int seedIdx=startQIdx; seedIdx<endQIdx; seedIdx++) { 
        syntenies.push_back(SyntenicSeeds(seeds[seedIdx])); // Initialize new synteny with the latest seed and add to existing candidates
        for(int candIdx=0; candIdx<syntenies.isize()-1; candIdx++) { //For all the previously added synteny candidates try concatanating new seed
            syntenies[candIdx].addSeed(seeds[seedIdx]); 
        }
    }
    svec<SyntenicSeeds>::iterator ssIter = max_element(syntenies.begin(), syntenies.end());
    return *ssIter;
}
*/
void FastAlignUnit::findAllSeeds(int numThreads, double identThresh) {
    int totSize   = m_targetSeqs.getNumSeqs();
    if(numThreads>totSize) { numThreads = totSize; }

    FILE_LOG(logDEBUG1) << "Finding Seeds";
    cout << "Finding All Seeds..." << endl;

    m_seeds.resize(totSize);             // Make sure enough memory is declared 
    ThreadQueueVec threadQueue(totSize); // Use for queueing instances threads should handle
    ThreadHandler th;

    for (int i=0; i<numThreads; i++) {
        char tmp[256];
        sprintf(tmp, "%d", i);
        string init = "init_findSeeds";
        init += tmp;
        th.AddThread(new FindSeedsSingleThread(threadQueue, *this, i));
        th.Feed(i, init);
    }
    while (!th.AllDone()) {
        usleep(10000);
    }

    m_seeds.sortSeeds(); //Sorts Seeds

    cout << "\r===================== " << "100.0% " << flush; 
    cout << "Completed finding Seeds." << endl;
}

void FastAlignUnit::alignSequence(int targetSeqIdx, ostream& sOut, ThreadMutex& mtx) const {
    svec<SyntenicSeeds> candidSynts;
    findSyntenicBlocks(targetSeqIdx, candidSynts); 
    for(int i=0; i<candidSynts.isize(); i++) {
        FILE_LOG(logDEBUG4) << candidSynts[i].toString();
//cout<<"Indel size: " << candidSynts[i].getMaxCumIndelSize()<<endl;
//cout<<"Seed Count: " << candidSynts[i].getNumSeeds()<<endl;
//cout<<"Seed Cover: " << candidSynts[i].getSeedCoverage(30)<<endl;
        int queryIdx = candidSynts[i].getQueryIdx();
        if(m_targetSeqs[targetSeqIdx].Name() == getQuerySeq(queryIdx).Name()) { continue; }
        Cola cola1 = Cola();
        DNAVector target = m_targetSeqs[targetSeqIdx];
        DNAVector query   = getQuerySeq(queryIdx);
//cout<<"Alignment Range: "<<candidSynts[i].getInitTargetOffset()<<"   "<<candidSynts[i].getInitQueryOffset()<<"  "<<candidSynts[i].getLastTargetIdx()<< "   " << candidSynts[i].getLastQueryIdx() << endl;
//cout<<"Syntenic seeds:  " <<candidSynts[i].toString()<<endl;
        int colaIndent = candidSynts[i].getMaxCumIndelSize();
        if(colaIndent>500) { colaIndent = 500; } //This is a temporary hack until the bandwidth function in Cola is fixed so that it doesn't require cumulative indel size
//cout<<"Cola Indent: " << colaIndent<<endl;
        cola1.createAlignment(target, query, AlignerParams(colaIndent, SWGA), candidSynts[i].getInitTargetOffset(), 
                              candidSynts[i].getInitQueryOffset(), candidSynts[i].getLastTargetIdx(), candidSynts[i].getLastQueryIdx());
        Alignment& cAlgn = cola1.getAlignment();
        if(cAlgn.getIdentityScore()>=m_params.getMinIdentity()) {
//cout<<"Aligned **********************************************"<<endl;
            mtx.Lock();
            sOut << target.Name() << " vs " << query.Name() << endl;
            cAlgn.print(0,1,sOut,100);
            mtx.Unlock();
        } else {
            FILE_LOG(logDEBUG4) <<"No Alignment at given significance threshold";
        }
    }   
}


void FastAlignUnit::alignAllSeqs(int numThreads, ostream& sOut) {
    int totSize   = m_targetSeqs.getNumSeqs();
    if(numThreads>totSize) { numThreads = totSize; }

    FILE_LOG(logDEBUG1) << "Aligning "; 
    cout << "Finding Syntenic seeds and aligning sequences..." << endl;

    ThreadQueueVec threadQueue(totSize); // Use for queueing instances threads should handle
    ThreadHandler th;
    ThreadMutex mtx;

    for (int i=0; i<numThreads; i++) {
        char tmp[256];
        sprintf(tmp, "%d", i);
        string init = "init_alignment_";
        init += tmp;
        th.AddThread(new AlignmentSingleThread(threadQueue, *this, i, sOut, mtx));
        th.Feed(i, init);
    }
    while (!th.AllDone()) {
        usleep(10000);
    }

    cout << "\r===================== " << "100.0% " << flush; 
    cout << "Completed aligning sequences." << endl;
}
//======================================================

//======================================================
int FastAlignQueryUnit::findSeeds(const DNAVector& targetSeq, int seedSizeThresh, SeedArray& seedArray) const { 
    map<unsigned long, int> stringsUsed_curr;                    // Flagset for strings that have been searched for a given extension
    const svec<SuffixArrayElement>& suffixes = m_suffixes.getSuffixes();
    for(int targetIterPos=0; targetIterPos<=targetSeq.isize()-seedSizeThresh; targetIterPos++) {
        FILE_LOG(logDEBUG4)  << "Iterating position in string: "<< targetIterPos;
        svec<SuffixArrayElement>::const_iterator fIt = lower_bound(suffixes.begin(), suffixes.end(), 
                                                                   SequenceWithOffset<DNAVector>(targetSeq, targetIterPos), 
                                                                   SuffixArray<DNASeqs, DNAVector>::CmpSuffixArrayElementOL(m_suffixes));
        svec<SuffixArrayElement>::const_reverse_iterator rIt(fIt);
        bool keepLooking = true;
        for (;keepLooking && fIt!=suffixes.end(); fIt++) {
            keepLooking = handleIterInstance<svec<SuffixArrayElement>::const_iterator>(fIt, stringsUsed_curr, targetSeq, targetIterPos, seedSizeThresh, seedArray);
        }
        keepLooking = true;
        for (;keepLooking && rIt!=suffixes.rend(); rIt++) {
            keepLooking = handleIterInstance<svec<SuffixArrayElement>::const_reverse_iterator>(rIt, stringsUsed_curr, targetSeq, targetIterPos, seedSizeThresh, seedArray);
        }
    }
    return seedArray.getNumSeeds();
}

template<class IterType>
bool FastAlignQueryUnit::handleIterInstance(IterType iter, map<unsigned long, int>& stringsUsed_curr, const DNAVector& targetSeq,
                                    int targetIterPos, 
                                    int seedSizeThresh, SeedArray& seedArray) const {
    map<unsigned long, int>::iterator usedIter = stringsUsed_curr.find((*iter).getIndex());
    if(usedIter!=stringsUsed_curr.end() && usedIter->second>targetIterPos) {  //Check if string has already been used
        FILE_LOG(logDEBUG4)  << usedIter->first << "    " << (*iter).getOffset() << " has already been found as seed";
        return true; //continue 
    }
    FILE_LOG(logDEBUG4)  << "Checking for match: "; 
    int matchLength = checkInitMatch(targetSeq, targetIterPos, *iter, seedSizeThresh);
    if(matchLength  < seedSizeThresh) { return false; } //Break out of looping! Suitable seed was not found
    int contactPos = targetIterPos; 
    // Synchronized version of overlap adding that locks so no iterference occurs with other threads
    seedArray.addSeed((*iter).getIndex(), (*iter).getOffset(), contactPos, matchLength);

    FILE_LOG(logDEBUG3)  << "Adding seed: " << "\t" << (*iter).getIndex() << "\t" << contactPos
                         << "\t" << (*iter).getOffset() << "\t" << matchLength;
    stringsUsed_curr[(*iter).getIndex()] = contactPos + matchLength; //Record that this sequence has been covered with seeds upto this index
    return true;
}


int FastAlignQueryUnit::checkInitMatch(const DNAVector& targetSeq, int targetOffset, const SuffixArrayElement& extSeqSA, int seedSizeThresh) const {
    int origSize = targetSeq.isize();
    int extSize  = m_suffixes.getStringSize(extSeqSA.getIndex());
    if(origSize < seedSizeThresh || extSize < seedSizeThresh) { return -2; }  // Pre-check 

    int idx2    = extSeqSA.getIndex();
    int offset2 = extSeqSA.getOffset();

    const DNAVector& d2 = m_suffixes.getString(idx2);

    int limit = min(targetSeq.isize()-targetOffset, m_suffixes.getStringSize(idx2)-offset2);
    int seedSize=0;
    for(; seedSize<limit; seedSize++) {
        if(targetSeq[targetOffset+seedSize]!=d2[offset2+seedSize]) {
            break;
        }
    }
    FILE_LOG(logDEBUG4) << "Check seed match size: " << seedSize;  
    return seedSize;
}
//======================================================

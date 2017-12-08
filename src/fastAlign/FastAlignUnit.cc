#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "ryggrad/src/base/StringUtil.h"
#include "ryggrad/src/base/Logger.h"
#include "ryggrad/src/base/RandomStuff.h"
#include "SeedingThreads.h"
#include "AlignmentThreads.h"
#include "FastAlignUnit.h"

//======================================================

void FastAlignUnit::findSeeds(int querySeqIdx) {
    m_targetUnit.findSeeds(m_querySeqs[querySeqIdx], m_params.getSeedSize(), m_seeds[querySeqIdx]);
}

void FastAlignUnit::findSyntenicBlocks(int querySeqIdx, svec<SyntenicSeeds>& maxSynts) const {
    const SeedArray& seeds = m_seeds[querySeqIdx];
    int currTIdx     = -1;
    int prevIdx      =  0;  // Keep track of the previous index where the target index changed
    for(int seedIdx=0; seedIdx<seeds.getNumSeeds(); seedIdx++) {
        int tIdx = seeds[seedIdx].getTargetIdx();
        if((currTIdx != tIdx && currTIdx != -1) || seedIdx==seeds.getNumSeeds()-1) { //New target sequence 
            // Find the best synteny & save if seed coverage of sequence passes acceptance threshold
            int startIdx = prevIdx;
            int endIdx   = seedIdx -1; //inclusive index
            if(seedIdx==seeds.getNumSeeds()-1) { endIdx++; } //Last Item should is covered for
            const SyntenicSeeds& ss = searchDPSynteny(seeds, startIdx, endIdx);
            if(ss.getSeedCoverage(m_params.getSeedSize()*2) > m_params.getMinSeedCover()) {
                FILE_LOG(logDEBUG2) << "Syntenic seeds where seed coverage passes thereshold: ";
                maxSynts.push_back(ss);
            } else {
                FILE_LOG(logDEBUG2) << "Syntenic seeds where seed coverage doesn't pass thereshold: " 
                                    << ss.getSeedCoverage(m_params.getSeedSize()*2);
            }
            FILE_LOG(logDEBUG3) << ss.toString();
            prevIdx = seedIdx;
        }
        currTIdx = tIdx; 
    }
} 

SyntenicSeeds FastAlignUnit::searchDPSynteny(const SeedArray& seeds, int startTIdx, int endTIdx) const {
    SyntenicSeedFinder sFinder(SeedsSubset(seeds, startTIdx, endTIdx)); 
    return sFinder.searchDP();
}

void FastAlignUnit::findAllSeeds(int numThreads, double identThresh) {
    int totSize   = m_querySeqs.getNumSeqs();
    if(numThreads>totSize) { numThreads = totSize; }

    FILE_LOG(logINFO) << "Finding Seeds";
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

void FastAlignUnit::alignSequence(int querySeqIdx, svec<AlignmentInfo>& cAlignmentInfos) const { 
    ThreadMutex mtx;
    ofstream sOut;
    alignSequence(querySeqIdx, cAlignmentInfos, 0, 1,  sOut, mtx);
}

void FastAlignUnit::alignSequence(int querySeqIdx, ostream& sOut , ThreadMutex& mtx) const {
    svec<AlignmentInfo>  cAlgnInfos;
    alignSequence(querySeqIdx, cAlgnInfos, 1, 0,  sOut, mtx); 
} 

void FastAlignUnit::alignSequence(int querySeqIdx, svec<AlignmentInfo>& cAlignmentInfos,
                                  int printResults, int storeAlignmentInfo, ostream& sOut, ThreadMutex& mtx) const {
    svec<SyntenicSeeds> candidSynts;
    findSyntenicBlocks(querySeqIdx, candidSynts); 
    cAlignmentInfos.reserve(candidSynts.isize());
    for(int i=0; i<candidSynts.isize(); i++) {
        FILE_LOG(logDEBUG3) << " Aligning based on candidate syntenic seed set: " << candidSynts[i].toString();
        FILE_LOG(logDEBUG3) << "Indel size: " << candidSynts[i].getMaxCumIndelSize() << "  Seed Count: " 
                            << candidSynts[i].getNumSeeds() << " Seed Cover: " << candidSynts[i].getSeedCoverage(m_params.getSeedSize()*2);
        int targetIdx = candidSynts[i].getTargetIdx();
        if(m_querySeqs[querySeqIdx].Name() == getTargetSeq(targetIdx).Name()) { continue; }
        Cola cola1 = Cola();
        DNAVector query, target;
        query.SetToSubOf(m_querySeqs[querySeqIdx], candidSynts[i].getInitQueryOffset());
        target.SetToSubOf(getTargetSeq(targetIdx), candidSynts[i].getInitTargetOffset());
        query.SetName(m_querySeqs[querySeqIdx].Name());
        target.SetName(getTargetSeq(targetIdx).Name());
        FILE_LOG(logDEBUG3) << "Alignment Range: "<<candidSynts[i].getInitQueryOffset()<<"   "<<candidSynts[i].getInitTargetOffset()
                            <<"  "<<candidSynts[i].getLastQueryIdx()<< "   " << candidSynts[i].getLastTargetIdx() << endl;
        int colaIndent = candidSynts[i].getMaxCumIndelSize();
        FILE_LOG(logDEBUG3) << " Aligning " << query.Name() << " vs. " << target.Name();
        FILE_LOG(logDEBUG3) << " with cola Indent: " << colaIndent << " capped at 500 and inital query offset: " << candidSynts[i].getInitQueryOffset() 
                            << " initial target offset: " << candidSynts[i].getInitTargetOffset();
        if(colaIndent>500) { colaIndent = 500; }
        cola1.createAlignment(query, target, AlignerParams(colaIndent, SWGA));
        if(storeAlignmentInfo) {
          cAlignmentInfos.push_back(cola1.getAlignment().getInfo());
          cAlignmentInfos.back().setSeqAuxInfo(candidSynts[i].getInitQueryOffset(), candidSynts[i].getInitTargetOffset(), true, true); //TODO pass in the strand from function calling alignSequence
        }
        if(printResults) {
          Alignment& tempAlgn = cola1.getAlignment();
          tempAlgn.setSeqAuxInfo(candidSynts[i].getInitQueryOffset(), candidSynts[i].getInitTargetOffset(), true, true); 
          writeAlignment(tempAlgn, sOut, mtx);
        }
    }   
}
void FastAlignUnit::writeAlignment(const Alignment& algn, ostream& sOut, ThreadMutex& mtx) const {
    if(algn.getIdentityScore()>=m_params.getMinIdentity()) {
        FILE_LOG(logDEBUG2) << " Aligned " << algn.getQueryName() << " vs. " << algn.getTargetName();
        mtx.Lock();
        if(m_revCmp) {
            sOut << algn.getQueryName() << "_RC" << " vs " << algn.getTargetName() << endl;
        } else {
            sOut << algn.getQueryName() << " vs " << algn.getTargetName() << endl;
        }
        algn.print(0,1,sOut,100);
        mtx.Unlock();
    } else {
        FILE_LOG(logDEBUG2) <<"No Alignment at given significance threshold";
    }
}   


void FastAlignUnit::alignAllSeqs(int numThreads, ostream& sOut) {
    int totSize   = m_querySeqs.getNumSeqs();
    if(numThreads>totSize) { numThreads = totSize; }

    FILE_LOG(logINFO) << "Aligning "; 
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
int FastAlignTargetUnit::findSeeds(const DNAVector& querySeq, int seedSizeThresh, SeedArray& seedArray) const { 
    map<unsigned long, int> stringsUsed_curr;                    // Flagset for strings that have been searched for a given extension
    const svec<SuffixArrayElement>& suffixes = m_suffixes.getSuffixes();
    for(int queryIterPos=0; queryIterPos<=querySeq.isize()-seedSizeThresh; queryIterPos++) {
        FILE_LOG(logDEBUG4)  << "Iterating position in string: "<< queryIterPos;
        svec<SuffixArrayElement>::const_iterator fIt = lower_bound(suffixes.begin(), suffixes.end(), 
                                                                   SequenceWithOffset<DNAVector>(querySeq, queryIterPos), 
                                                                   SuffixArray<DNASeqs, DNAVector>::CmpSuffixArrayElement(m_suffixes));
        FILE_LOG(logDEBUG4)  << "Searching for suffix - found lower-bound: " << (*fIt).toString() << endl;
        svec<SuffixArrayElement>::const_reverse_iterator rIt(fIt);
        bool keepLooking = true;
        for (;keepLooking && fIt!=suffixes.end(); fIt++) {
            keepLooking = handleIterInstance<svec<SuffixArrayElement>::const_iterator>(fIt, stringsUsed_curr, querySeq, queryIterPos, seedSizeThresh, seedArray);
        }
        keepLooking = true;
        for (;keepLooking && rIt!=suffixes.rend(); rIt++) {
            keepLooking = handleIterInstance<svec<SuffixArrayElement>::const_reverse_iterator>(rIt, stringsUsed_curr, querySeq, queryIterPos, seedSizeThresh, seedArray);
        }
    }
    return seedArray.getNumSeeds();
}

template<class IterType>
bool FastAlignTargetUnit::handleIterInstance(IterType iter, map<unsigned long, int>& stringsUsed_curr, const DNAVector& querySeq,
                                    int queryIterPos, 
                                    int seedSizeThresh, SeedArray& seedArray) const {
    map<unsigned long, int>::iterator usedIter = stringsUsed_curr.find((*iter).getIndex());
    if(usedIter!=stringsUsed_curr.end() && usedIter->second>queryIterPos) {  //Check if string has already been used
        FILE_LOG(logDEBUG4)  << usedIter->first << "    " << (*iter).getOffset() << " has already been found as seed";
        return true; //continue 
    }
    int matchLength = checkInitMatch(querySeq, queryIterPos, *iter, seedSizeThresh);
    FILE_LOG(logDEBUG4) << "Check seed match size: " << matchLength;
    if(matchLength  < seedSizeThresh) { return false; } //Break out of looping! Suitable seed was not found
    int contactPos = queryIterPos; 
    seedArray.addSeed((*iter).getIndex(), (*iter).getOffset(), contactPos, matchLength);
    FILE_LOG(logDEBUG3)  << "Adding seed: " << "\t" << (*iter).getIndex() << "\t" << contactPos
                         << "\t" << (*iter).getOffset() << "\t" << matchLength;
    stringsUsed_curr[(*iter).getIndex()] = contactPos + matchLength; //Record that this sequence has been covered with seeds upto this index
    return true;
}


int FastAlignTargetUnit::checkInitMatch(const DNAVector& querySeq, int queryOffset, const SuffixArrayElement& extSeqSA, int seedSizeThresh) const {
    int origSize = querySeq.isize();
    int extSize  = m_suffixes.getStringSize(extSeqSA.getIndex());
    if(origSize < seedSizeThresh || extSize < seedSizeThresh) { return -2; }  // Pre-check 

    int idx2    = extSeqSA.getIndex();
    int offset2 = extSeqSA.getOffset();

    const DNAVector& d2 = m_suffixes.getString(idx2);

    int limit = min(querySeq.isize()-queryOffset, m_suffixes.getStringSize(idx2)-offset2);
    int seedSize=0;
    for(; seedSize<limit; seedSize++) {
        if(querySeq[queryOffset+seedSize]!=d2[offset2+seedSize]) {
            break;
        }
    }
    return seedSize;
}
//======================================================

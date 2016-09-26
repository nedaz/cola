#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#if defined(OPEN_MP)
  #include <parallel/algorithm>
#else
  #include <algorithm>
#endif

#include <string>
#include <sstream>
#include <stdint.h>
#include "base/SVector.h"
#include "src/DNAVector.h"
#include "base/FileParser.h"
#include "extern/logger/log.h"
#include "src/FastAlign/SeedingObjects.h"

//======================================================
void SeedCandid::set(int qI, int qO, int tO, int l) {
    m_queryIdx     = qI;
    m_queryOffset  = qO;
    m_targetOffset = tO;
    m_length       = l;
}


/** Sorted based on the primary ordering of query Index and secondary ordering of target offset **/
bool SeedCandid::operator < (const SeedCandid & sC) const {
    if (getQueryIdx() != sC.getQueryIdx()) {
        return (getQueryIdx() < sC.getQueryIdx());
    }
    if (getTargetOffset() != sC.getTargetOffset()) {
        return (getTargetOffset() < sC.getTargetOffset());
    }
    return (getQueryOffset() < sC.getQueryOffset()); 
}

string SeedCandid::toString() const {
    stringstream ss;
    ss << getQueryIdx() << "\t" 
       << getTargetOffset() << "\t" 
       << getQueryOffset() << "\t" 
       << getSeedLength() << "\t"; 
    return ss.str();
}
//======================================================

//======================================================
string SeedArray::toString() const {
    stringstream ss;
    for(int j=0; j<getNumSeeds(); j++) {
        ss << m_seeds[j].toString() << endl;  
    }
    return ss.str();
}
//======================================================

//======================================================
void AllSeedCandids::write(const string& seedFile, int mode) const {
    if(mode==0) { writeBin(seedFile);   } 
    if(mode==1) { writeAsc(seedFile);   } 
    if(mode==2) { writeStats(seedFile); }
}

void AllSeedCandids::writeAsc(const string& seedFile) const {
    ofstream sout;
    sout.open(seedFile.c_str());
    int totNumOfSeqs = m_seeds.isize();
    sout << totNumOfSeqs << endl;
    for(int i=0; i<totNumOfSeqs; i++) {
        sout << getSeedString(i);
    }
    sout.close();
}

void AllSeedCandids::writeBin(const string& seedFile) const {
    CMWriteFileStream fs;
    fs.Open(seedFile.c_str());
    int totNumOfSeqs = m_seeds.isize();
    fs.Write(totNumOfSeqs); // To use for allocating memory when reading
    for(int i=0; i<totNumOfSeqs; i++) {
        const SeedArray seeds = m_seeds[i];
        for(int j=0; j<seeds.getNumSeeds(); j++) {
            fs.Write(i);
            fs.Write(seeds[j].getQueryIdx());
            fs.Write(seeds[j].getQueryOffset());
            fs.Write(seeds[j].getTargetOffset());
            fs.Write(seeds[j].getSeedLength());
        }
    }
    fs.Close();
}

string AllSeedCandids::getSeedString(int index) const {
    stringstream ss;
    const SeedArray seeds = m_seeds[index];
    for(int j=0; j<seeds.getNumSeeds(); j++) {
        ss << index  << "\t" << seeds[j].toString() << endl;  
    }
    return ss.str();
}

void AllSeedCandids::writeStats(const string& statFile) const {
    ofstream sout;
    sout.open(statFile.c_str());
    int totNumOfSeqs = m_seeds.isize();
    for(int i=0; i<totNumOfSeqs; i++) {
        sout << getNumSeeds(i) << endl;
    }
    sout.close();
}


void AllSeedCandids::addSeedFromString(const string& strIn){
    CMTokenizer tokenizer;
    tokenizer.AddDelimiter("\t");
    CMPtrStringList tokens;
    tokenizer.Tokenize(tokens, strIn.c_str());
    if(tokens.length()<6) { 
        FILE_LOG(logERROR) << "Wrong overlap file format - six columns required"; 
        return;
    }
    addSeed(atoi((const char*)*tokens[0]), atoi((const char*)*tokens[1]),
            atoi((const char*)*tokens[2]), atoi((const char*)*tokens[3]), 
            atoi((const char*)*tokens[4]));
}

void AllSeedCandids::sortSeeds() {
    for (int i=0; i<m_seeds.isize(); i++) {
        m_seeds[i].sortSeeds();
    }
}
//======================================================

//======================================================
void SyntenicSeeds::initSeed(const SeedCandid& sC) {
    this->m_seeds.push_back(sC);
    m_totalSize     = sC.getSeedLength();
    m_maxIndelSize  = 0;
    m_cumIndelSize  = 0;
}

bool SyntenicSeeds::checkSeedSynt(const SeedCandid& sC) { 
    const SeedCandid& currSeed = getLatestSeed();
    if(sC.getQueryOffset()>currSeed.getQueryOffset() && sC.getTargetOffset()>currSeed.getTargetOffset()) {
        return true;
    } 
    return false;
}

bool SyntenicSeeds::addSeed(const SeedCandid& sC) { 
    if(checkSeedSynt(sC)) {
        addSeedOrdered(sC);
        return true; //Candidate has been added
    } else {
        return false; //Candidate has not been added
    }
}

void SyntenicSeeds::addSeedOrdered(const SeedCandid& sC) { 
    const SeedCandid& prevSeed = getLatestSeed();
    int currSeedSize           = sC.getSeedLength();
    m_totalSize                += currSeedSize; //Add the seed length of the latest seed to the existing total
    //Update the max indent with the latest added seed 
    int latestIndent = (sC.getTargetOffset()-prevSeed.getTargetOffset()) 
                          -(sC.getQueryOffset()-prevSeed.getQueryOffset());
    if(abs(latestIndent) > m_maxIndelSize) { m_maxIndelSize = abs(latestIndent); }
    m_cumIndelSize += latestIndent;
    this->m_seeds.push_back(sC);
} 

float SyntenicSeeds::getSeedCoverage(int singleSeedThresh) { 
    if( getNumSeeds()<1)                                                      { return 0; } 
    if(getNumSeeds()==1 && this->m_seeds[0].getSeedLength()<singleSeedThresh) { return 0; }

    float queryRange  = getLatestSeed().getQueryOffset() - this->m_seeds[0].getQueryOffset();
    float targetRange =  getLatestSeed().getTargetOffset() - this->m_seeds[0].getTargetOffset();
    return getTotalSeedLength()/(max(queryRange, targetRange)+getLatestSeed().getSeedLength());
}

int SyntenicSeeds::getMaxCumIndelSize() const { 
    if(m_maxIndelSize>abs(m_cumIndelSize)) { return m_maxIndelSize;         } 
    else                                   { return abs(m_cumIndelSize);    }
}

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
#include "ryggrad/src/base/SVector.h"
#include "ryggrad/src/general/DNAVector.h"
#include "ryggrad/src/base/FileParser.h"
#include "ryggrad/src/base/Logger.h"
#include "SeedingObjects.h"

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
string SeedsSubset::toString() const {
    stringstream ss;
    for(int i=m_startIdx; i<m_endIdx; i++) {
        ss << m_seeds[i].toString() << endl;  
    }
    return ss.str();
}
//======================================================



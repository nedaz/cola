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
#include "SyntenicSeeds.h"

//======================================================
void SyntenicSeeds::initSeed(const SeedCandid& sC) {
    setTargetIndex(sC.getTargetIdx());
    this->m_seeds.push_back(sC);
    m_totalSize     = sC.getSeedLength();
    m_maxIndelSize  = 0;
    m_cumIndelSize  = 0;
}

bool SyntenicSeeds::checkSeedSynt(const SeedCandid& sC) { 
    const SeedCandid& currSeed = getLatestSeed();
    if(sC.getTargetOffset()>currSeed.getTargetOffset() && sC.getQueryOffset()>currSeed.getQueryOffset()) {
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
    if(!hasSeed()) { 
        initSeed(sC); 
        return;
    }
    const SeedCandid& prevSeed = getLatestSeed();
    int currSeedSize           = sC.getSeedLength();
    m_totalSize                += currSeedSize; //Add the seed length of the latest seed to the existing total
    //Update the max indent with the latest added seed 
    int latestIndent = (sC.getQueryOffset()-prevSeed.getQueryOffset()) 
                          -(sC.getTargetOffset()-prevSeed.getTargetOffset());
    if(abs(latestIndent) > m_maxIndelSize) { m_maxIndelSize = abs(latestIndent); }
    m_cumIndelSize += latestIndent;
    this->m_seeds.push_back(sC);
} 

float SyntenicSeeds::getSeedCoverage(int singleSeedThresh) const { 
    if( getNumSeeds()<1)                                            { return 0;  } 
    if(getNumSeeds()==1 && getTotalSeedLength()<singleSeedThresh)   { return -1; }
    if(getNumSeeds()==2 && getTotalSeedLength()<singleSeedThresh) { return -2; }

    float targetRange  = getLatestSeed().getTargetOffset() - this->m_seeds[0].getTargetOffset();
    float queryRange   =  getLatestSeed().getQueryOffset() - this->m_seeds[0].getQueryOffset();
    return (float) getTotalSeedLength()/(max(targetRange, queryRange)+getLatestSeed().getSeedLength());
}

int SyntenicSeeds::getMaxCumIndelSize() const { 
    if(m_maxIndelSize>abs(m_cumIndelSize)) { return m_maxIndelSize;         } 
    else                                   { return abs(m_cumIndelSize);    }
}
//======================================================

//======================================================
string SSSearchNode::toString() const {
    stringstream ss;
    ss << "Seed Idx: " << m_currSeedIdx << " Best predecessor seed Idx: " << m_bestPredSeedIdx
       << " Score: "   <<  m_score; 
    return ss.str();
}
//======================================================

//======================================================
void SyntenicSeedFinder::setupAdjTable(svec< svec<int> >& adjInvTable) const { 
    adjInvTable.resize(m_seeds.isize());
    for(int i=0; i<m_seeds.isize(); i++) {
        const SeedCandid& sC1 = m_seeds[i];
        for(int j=0; j<i; j++) {
            const SeedCandid& sC2 = m_seeds[j];
            if(sC1.getTargetOffset()>=sC2.getTargetOffset()+sC2.getSeedLength() 
               && sC1.getQueryOffset()>=sC2.getQueryOffset()+sC2.getSeedLength()) {
                adjInvTable[i].push_back(j);
            }
        }
    }
}

SyntenicSeeds SyntenicSeedFinder::searchDP() const {
    if(m_seeds.isize()==0) { return SyntenicSeeds(); } 
    //1a. Declare the search Nodes
    svec<SSSearchNode>   searchNodes;   
    searchNodes.resize(m_seeds.isize());
    //1b. Create synteny graph (adjacency table)
    svec< svec<int> > adjInvTable;
    setupAdjTable(adjInvTable);
    //2. Start the dynamic programming search
    SSSearchNode bestNode = SSSearchNode();
    SSSearchNode currNode;
    for(int seedIdx=0; seedIdx<m_seeds.isize(); seedIdx++) { 
        //2a. Find best predecessor and hence score for current node
        int bestPredIdx   = -1;
        int bestPredScore =  0;
        for(int i=0; i<adjInvTable[seedIdx].isize(); i++) {
            int currIdx   = adjInvTable[seedIdx][i];
            int currScore = searchNodes[currIdx].getScore(); 
            if(currScore > bestPredScore) {
                bestPredIdx   = currIdx; 
                bestPredScore = currScore;
            }
        }
        currNode = SSSearchNode(seedIdx, bestPredIdx, bestPredScore+m_seeds[seedIdx].getSeedLength());
        searchNodes[seedIdx] =  currNode;
        FILE_LOG(logDEBUG4) << "Adding synteny search node at: " << seedIdx << " " << currNode.toString(); 
        // 2b. Set the best node and continue 
        if(currNode.getScore() > bestNode.getScore()) {
            bestNode = currNode;
        }
    }
    // 3. Backtrack from the best node to find the best path of nodes and return as a SyntenicSeeds object
    return( backtrack(bestNode, searchNodes) );
}

SyntenicSeeds SyntenicSeedFinder::backtrack(SSSearchNode& bestNode, svec<SSSearchNode>& searchNodes) const {
    svec<int> revSeedIdxs; // Reverse list of seeds
    int currNodeIdx = bestNode.getSeedIdx();
    while(currNodeIdx != -1) {
        revSeedIdxs.push_back(currNodeIdx);
        currNodeIdx = searchNodes[currNodeIdx].getBestPredSeedIdx();
    }
    // Create the syntenicSeeds from the reverse list of seed ids
    SyntenicSeeds ss;
    for (svec<int>::reverse_iterator revSeedIter = revSeedIdxs.rbegin(); revSeedIter != revSeedIdxs.rend(); ++revSeedIter) { 
        ss.addSeedOrdered(SeedCandid(m_seeds[*revSeedIter]));
    }  
    return ss;
}
//======================================================

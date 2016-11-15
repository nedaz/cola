
#ifndef _SYNTENIC_SEEDS_H_
#define _SYNTENIC_SEEDS_H_

#include "ryggrad/src/base/ThreadHandler.h"
#include "SeedingObjects.h" 

//======================================================
/** Class specifying set of seeds that 
    are in order and syntenic (i.e.) They can form a valid alignment 
*/
class SyntenicSeeds: public SeedArray 
{
public:
    SyntenicSeeds() : m_targetIdx(-1), m_totalSize(0), m_maxIndelSize(0), m_cumIndelSize(0) {}
    SyntenicSeeds(const SeedCandid& sC) : m_targetIdx(sC.getTargetIdx()), m_totalSize(0), m_maxIndelSize(0), m_cumIndelSize(0) { 
        initSeed(sC);
    }

    int  getTargetIdx() const          { return m_targetIdx;     }
    int  getTotalSeedLength() const    { return m_totalSize;     }
    int  getMaxIndelSize() const       { return m_maxIndelSize;  } 
     int getMaxCumIndelSize() const;  // Return the maximum between the cumulative and none cumulative indel size
    void setTargetIndex(int tIdx)      { m_targetIdx = tIdx;      }
    
    /** Returns the start and end of the seeds for target or query, whichever is smaller, 
        returns 0 for single seeds that are shorter than the given threshold */ 
    float getSeedCoverage(int singleSeedThresh) const; 

    /** Functions that return the offset of the first seed in the synteny */
    int getInitTargetOffset() const     { return this->m_seeds[0].getTargetOffset();  }
    int getInitQueryOffset() const    { return this->m_seeds[0].getQueryOffset(); }

    /** Functions that return the offset plus length of the last seed in the synteny */
    int getLastTargetIdx() const     { return getLatestSeed().getTargetOffset() + getLatestSeed().getSeedLength();  }
    int getLastQueryIdx() const    { return getLatestSeed().getQueryOffset() + getLatestSeed().getSeedLength(); }

    /** Functions that return the offset plus length of the last seed in the synteny */
    int getTargetCoverLength() const     { return (getLastTargetIdx()-getInitTargetOffset());   }
    int getQueryCoverLength() const    { return (getLastQueryIdx()-getInitQueryOffset()); } 

    bool operator<(SyntenicSeeds other) const { return getTotalSeedLength() < other.getTotalSeedLength(); } 

    /** Initialize with the first Seed. Note that this function must be called with first seed before adding any otheres */
    void initSeed(const SeedCandid& sc);
    /** Check to see if the given seed is syntenic with the existing seeds */
    bool checkSeedSynt(const SeedCandid& sc);
    /** Add seed that has been checked for synteny and is in increasing offset order */
    void addSeedOrdered(const SeedCandid& sC); 
    /** Add seed & check for synteny before adding, returns true if added & false otherwise */
    bool addSeed(const SeedCandid& sC); 
    /** Add seed & check for synteny before adding */
    void addSeed(int targetIndex, int targetOffset, int queryOffset, int length) {
        this->addSeed(SeedCandid(targetIndex, targetOffset, queryOffset, length));
    } 

private: 
    /** Must call hasSeed & make sure seeds exist before using this function */
    const SeedCandid& getLatestSeed() const { return this->m_seeds[this->m_seeds.isize()-1];   }
    bool hasSeed()                          { return !this->m_seeds.empty();                   } 

    int     m_targetIdx;              /// Index of target to which these seeds belong
    int     m_totalSize;             /// The sum of the length of all seeds 
    int     m_maxIndelSize;          /// The longest indel size that the list of seeds entail 
    int     m_cumIndelSize;          /// The accumulative indel length
};
//======================================================

//======================================================
//======================================================
/** Syntenic seed finder search node */
class SSSearchNode 
{
public:
    SSSearchNode() : m_currSeedIdx(-1), m_bestPredSeedIdx(-1), m_score(-1) {}
    SSSearchNode(int currId, int bestPre, int score) : m_currSeedIdx(currId), m_bestPredSeedIdx(bestPre), m_score(score) {}
    int getSeedIdx() const            { return m_currSeedIdx;      }
    int getBestPredSeedIdx() const    { return m_bestPredSeedIdx;  }
    int getScore() const              { return m_score;            }
    string toString() const; 

private:
   int m_currSeedIdx;               /// The index of the seed to which this node relates
   int m_bestPredSeedIdx;           /// Index of the best scoring predecessor
   int m_score;                     /// Score at this node

    
};

//======================================================
/** 
  Finding the highest scoring (largest total seed length count) subset of syntenic seeds from a given set
  is equivalent to finding the longest path in a directed acyclic graph where the nodes are topologically
  ordered. This requires setting up the adjacency table and its reverse and carrying out a dyanmic program.
  These are the operations that this class handles.
*/
class SyntenicSeedFinder
{
public:
    SyntenicSeedFinder(const SeedsSubset& sS) : m_seeds(sS) {}
    SyntenicSeeds searchDP() const;

private:
    /** function to set up the adjacency table for the seed relations */ 
    void setupAdjTable(svec< svec<int> >& adjInvTable) const;
    /** Backtrack from the best node to find the best path of nodes and return as a SyntenicSeeds object */
    SyntenicSeeds backtrack(SSSearchNode& bestNode, svec<SSSearchNode>& searchNodes) const; 

    SeedsSubset          m_seeds;   /// Seeds that highest scoring syntenic subset will be chosen from 
};
//======================================================

#endif //_SYNTENIC_SEEDS_H_


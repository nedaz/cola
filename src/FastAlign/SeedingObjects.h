
#ifndef _SEEDING_OBJECTS_H_
#define _SEEDING_OBJECTS_H_

#include "base/ThreadHandler.h"

//======================================================
/** Class specifying coordinates for a possible seed to use in anchoring alignments */
class SeedCandid 
{
public:
    SeedCandid() : m_queryIdx(-1), m_queryOffset(-1), m_targetOffset(-1), m_length(-1) {};
    SeedCandid(int qI, int qO, int tO, int l) {
      set(qI, qO, tO, l);
    }

    void set(int qI, int qO, int tO, int l); 

    string toString() const;

    int     getQueryIdx() const           { return m_queryIdx;      }  
    int     getQueryOffset() const        { return m_queryOffset;   }  
    int     getTargetOffset() const       { return m_targetOffset;  }  
    int     getSeedLength() const         { return m_length;        }  

    bool operator < (const SeedCandid & rO) const; 

private: 
    int     m_queryIdx;       /// The index of the sequence to which this seeding match refers to
    int     m_queryOffset;    /// The position in the sequence where this overlap occurs from 
    int     m_targetOffset;   /// The position in the target sequence where this seeding match occurs from 
    int     m_length;         /// The length of the seeding match 
};
//======================================================

//======================================================
class SeedArray
{
public:
    SeedArray():m_seeds() {}

    const SeedCandid& operator[](int i) const { return m_seeds[i]; }
    SeedCandid& operator[](int i)             { return m_seeds[i]; }

    void addSeed(int queryIndex, int queryOffset, int targetOffset, int length) {
        SeedCandid sC = SeedCandid(queryIndex, queryOffset, targetOffset, length);
        m_seeds.push_back(sC);
    }

    int getNumSeeds() const { return m_seeds.isize(); }

    void sortSeeds() { sort(m_seeds.begin(), m_seeds.end()); }

    string toString() const; 

protected:
    svec<SeedCandid> m_seeds;     /// Vector of seeds that represents set of seeds belonging to one sequence 
};
//======================================================

//======================================================
class AllSeedCandids
{
public:
    AllSeedCandids():m_seeds(), m_mutex() {}

    AllSeedCandids(int size):m_seeds(), m_mutex() {resize(size);}
   
    const SeedArray& operator[](int i) const { return m_seeds[i]; }
    SeedArray& operator[](int i)             { return m_seeds[i]; }

    void addSeed(int targetIndex, int queryIndex, int targetOffset, int queryOffset, int length) {
        m_seeds[targetIndex].addSeed(queryIndex, queryOffset, targetOffset, length);
    }

    void addSeedSync(int targetIndex, int queryIndex, int targetOffset, int queryOffset, int length) {
        SeedCandid sC = SeedCandid(queryIndex, queryOffset, targetOffset, length);
        m_mutex.Lock();
        m_seeds[targetIndex].addSeed(queryIndex, queryOffset, targetOffset, length);
        m_mutex.Unlock();
    }

    void addSeedArraySync(int targetIndex, const SeedArray& sArray) {
        m_mutex.Lock();
        m_seeds[targetIndex] = sArray;
        m_mutex.Unlock();
    }

    int getNumSeeds(int seqIndex) const {
        return m_seeds[seqIndex].getNumSeeds();
    }

    void resize(int sz)  { m_seeds.resize(sz);     }
    int  getSize() const { return m_seeds.isize(); }

    /**  0: binary 1: ascii 2: statisitics */
    void write(const string& overlapFile, int mode) const;  
    void writeAsc(const string& readOverlapFile) const; 
    void writeBin(const string& readOverlapFile) const; 
    /** Writes out the number of overlaps per read */
    void writeStats(const string& statFile) const; 

    void sortSeeds();

private:
    /** Internal helper class for writing/serialization */
    string getSeedString(int index) const; 
    void addSeedFromString(const string& strIn);

    svec <SeedArray> m_seeds;       /// All seed candidates, one vector per target sequence 
    ThreadMutex  m_mutex;           /// To use for locking while assigning overlaps
};
//======================================================

//======================================================
/** Class specifying set of seeds that 
    are in order and syntenic (i.e.) They can form a valid alignment 
*/
class SyntenicSeeds: public SeedArray 
{
public:
    SyntenicSeeds() : m_queryIdx(-1), m_totalSize(0), m_maxIndelSize(0), m_cumIndelSize(0) {};
    SyntenicSeeds(const SeedCandid& sC) : m_queryIdx(sC.getQueryIdx()), m_totalSize(0), m_maxIndelSize(0), m_cumIndelSize(0) { 
        initSeed(sC);
    };

    int  getQueryIdx() const           { return m_queryIdx;      }
    int  getTotalSeedLength() const    { return m_totalSize;     }
    int  getMaxIndelSize() const       { return m_maxIndelSize;  } 
     int getMaxCumIndelSize() const;  // Return the maximum between the cumulative and none cumulative indel size
    void setQueryIndex(int qIdx)       { m_queryIdx = qIdx;      }
    
    /** Returns the start and end of the seeds for query or target, whichever is smaller, 
        returns 0 for single seeds that are shorter than the given threshold */ 
    float getSeedCoverage(int singleSeedThresh); 

    /** Functions that return the offset of the first seed in the synteny */
    int getInitQueryOffset() const     { return this->m_seeds[0].getQueryOffset();  }
    int getInitTargetOffset() const    { return this->m_seeds[0].getTargetOffset(); }

    /** Functions that return the offset plus length of the last seed in the synteny */
    int getLastQueryIdx() const     { return getLatestSeed().getQueryOffset() + getLatestSeed().getSeedLength();  }
    int getLastTargetIdx() const    { return getLatestSeed().getTargetOffset() + getLatestSeed().getSeedLength(); }

    /** Functions that return the offset plus length of the last seed in the synteny */
    int getQueryCoverLength() const     { return (getLastQueryIdx()-getInitQueryOffset());   }
    int getTargetCoverLength() const    { return (getLastTargetIdx()-getInitTargetOffset()); } 

    bool operator<(SyntenicSeeds other) const { return getTotalSeedLength() < other.getTotalSeedLength(); } 

    /** Initialize with the first Seed. Note that this function must be called with first seed before adding any otheres */
    void initSeed(const SeedCandid& sc);
    /** Check to see if the given seed is syntenic with the existing seeds */
    bool checkSeedSynt(const SeedCandid& sc);
    /** Add seed that has been checked for synteny and is in increasing offset order */
    void addSeedOrdered(const SeedCandid& sC); 
    /** Add seed & check for synteny before adding, returns true if addeded & false otherwise */
    bool addSeed(const SeedCandid& sC); 
    /** Add seed & check for synteny before adding */
    void addSeed(int queryIndex, int queryOffset, int targetOffset, int length) {
        this->addSeed(SeedCandid(queryIndex, queryOffset, targetOffset, length));
    } 

private: 
    /** Must call hasSeed & make sure seeds exist before using this function */
    const SeedCandid& getLatestSeed() const { return this->m_seeds[this->m_seeds.isize()-1];   }
    bool hasSeed()                          { return !this->m_seeds.empty();                   } 

    int     m_queryIdx;              /// Index of query to which these seeds belong
    int     m_totalSize;             /// The sum of the length of all seeds 
    int     m_maxIndelSize;          /// The longest indel size that the list of seeds entail 
    int     m_cumIndelSize;          /// The accumulative indel length
};
//======================================================


#endif //_SEEDING_OBJECTS_H_


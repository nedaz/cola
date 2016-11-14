
#ifndef _SEEDING_OBJECTS_H_
#define _SEEDING_OBJECTS_H_

#include "ryggrad/src/base/ThreadHandler.h"

//======================================================
/** Class specifying coordinates for a possible seed to use in anchoring alignments */
class SeedCandid 
{
public:
    SeedCandid() : m_targetIdx(-1), m_targetOffset(-1), m_queryOffset(-1), m_length(-1) {};
    SeedCandid(int qI, int qO, int tO, int l) {
      set(qI, qO, tO, l);
    }

    void set(int qI, int qO, int tO, int l); 

    string toString() const;

    int     getTargetIdx() const           { return m_targetIdx;      }  
    int     getTargetOffset() const        { return m_targetOffset;   }  
    int     getQueryOffset() const         { return m_queryOffset;    }  
    int     getSeedLength() const          { return m_length;         }  

    bool operator < (const SeedCandid & rO) const; 

private: 
    int     m_targetIdx;       /// The index of the sequence to which this seeding match refers to
    int     m_targetOffset;    /// The position in the sequence where this overlap occurs from 
    int     m_queryOffset;     /// The position in the query sequence where this seeding match occurs from 
    int     m_length;          /// The length of the seeding match 
};
//======================================================

//======================================================
class SeedArray
{
public:
    SeedArray():m_seeds() {}

    const SeedCandid& operator[](int i) const { return m_seeds[i]; }
    SeedCandid& operator[](int i)             { return m_seeds[i]; }

    void addSeed(int targetIndex, int targetOffset, int queryOffset, int length) {
        SeedCandid sC = SeedCandid(targetIndex, targetOffset, queryOffset, length);
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

    void addSeed(int queryIndex, int targetIndex, int queryOffset, int targetOffset, int length) {
        m_seeds[queryIndex].addSeed(targetIndex, targetOffset, queryOffset, length);
    }

    void addSeedSync(int queryIndex, int targetIndex, int queryOffset, int targetOffset, int length) {
        SeedCandid sC = SeedCandid(targetIndex, targetOffset, queryOffset, length);
        m_mutex.Lock();
        m_seeds[queryIndex].addSeed(targetIndex, targetOffset, queryOffset, length);
        m_mutex.Unlock();
    }

    void addSeedArraySync(int queryIndex, const SeedArray& sArray) {
        m_mutex.Lock();
        m_seeds[queryIndex] = sArray;
        m_mutex.Unlock();
    }

    int getNumSeeds(int seqIndex) const {
        return m_seeds[seqIndex].getNumSeeds();
    }

    void resize(int sz)  { m_seeds.resize(sz);     }
    int  getSize() const { return m_seeds.isize(); }

    /**  0: binary 1: ascii 2: statisitics */
    void write(const string& overlapFile, int mode) const;  
    void writeAsc(const string& seedFile) const; 
    void writeBin(const string& seedFile) const; 
    /** Writes out the number of overlaps per read */
    void writeStats(const string& statFile) const; 

    void sortSeeds();

private:
    /** Internal helper class for writing/serialization */
    string getSeedString(int index) const; 
    void addSeedFromString(const string& strIn);

    svec <SeedArray> m_seeds;       /// All seed candidates, one vector per query sequence 
    ThreadMutex  m_mutex;           /// To use for locking while assigning overlaps
};
//======================================================

//======================================================
class SeedsSubset
{
public:
     SeedsSubset(const SeedArray& sA, int sIdx, int eIdx) : m_seeds(sA), m_startIdx(sIdx), m_endIdx(eIdx) {}
     int getNumSeeds() const                      { return m_endIdx-m_startIdx;     }
     int isize() const                            { return getNumSeeds();           }
     const SeedCandid& operator[](int i) const    { return m_seeds[m_startIdx+i];   }
     string toString() const; 

private:
    const SeedArray& m_seeds;      /// Seeds that highest scoring syntenic subset will be chosen from 
    int m_startIdx;                /// Start index in the seed set where the relevant seeds start from
    int m_endIdx;                  /// End index in the seed set where the relevant seeds end at
};
//======================================================

#endif //_SEEDING_OBJECTS_H_


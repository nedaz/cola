#ifndef _ASSEMBLER_H_
#define _ASSEMBLER_H_

#include <sstream>
#include "ryggrad/src/base/ThreadHandler.h"
#include "AlignmentParams.h"
#include "SuffixArray.h"
#include "DNASeqs.h"
#include "SeedingObjects.h"
#include "SyntenicSeeds.h" 


//======================================================
class FastAlignTargetUnit
{
public:
    FastAlignTargetUnit(const string& inputFile, int stepSize): m_targetSeqs(inputFile), 
                                                               m_suffixes(m_targetSeqs, stepSize) {}

    const DNAVector&  getTargetSeq(int seqIdx) const { return m_targetSeqs[seqIdx];        }
    int getTargetSeqSize(int seqIdx) const           { return m_targetSeqs[seqIdx].size(); }
    const string& getTargetSeqName(int seqIdx) const { return m_targetSeqs[seqIdx].Name(); }

     /** Return a vector of SuffixArrayElement entry indexes for a given 
       string of those Substrings that share a significant subsequence 
       limit specifies the number of overlaps to limit the search to (limit=0 means set the limit to string size) 
     */
    int findSeeds(const DNAVector& querySeq, int seedSizeThresh, SeedArray& seedArray) const; 

private:
    /** Returns true to indicate that seed has been found and iterating should continue, false otherwise */ 
    template<class IterType>
    bool handleIterInstance(IterType iter, map<unsigned long, int>& stringsUsed_curr, const DNAVector& querySeq,
                            int queryIterPos, int seedSizeThresh, SeedArray& seedArray) const; 
    int checkInitMatch(const DNAVector& querySeq, int queryOffset, const SuffixArrayElement& extSeqSA, int seedSizeThresh) const; 


private:
    DNASeqs                          m_targetSeqs;      /// A list of sequences from which suffixes were constructed
    SuffixArray<DNASeqs, DNAVector>  m_suffixes;       /// Structure for creating and holding suffixes from the class DNA sequences
};

//======================================================

//Forward declaration
class FindSeedsThread;
class FindSeedsSingleThread;
class AlignmentThread;
class AlignmentSingleThread;

//======================================================
class FastAlignUnit
{
    friend class FindSeedsThread;
    friend class FindSeedsSingleThread;
    friend class AlignmentThread;
    friend class AlignmentSingleThread;

public:
    // Basic Constructor used for finding overlaps
    FastAlignUnit(const string& querySeqFile, const FastAlignTargetUnit& qUnit, const AlignmentParams& params, int numOfThreads, bool revCmp=false)
                  : m_querySeqs(querySeqFile), m_targetUnit(qUnit), 
                    m_params(params), m_seeds(m_querySeqs.getNumSeqs()), m_revCmp(revCmp) {
        if(revCmp) { m_querySeqs.reverseComplementAll(); }
        findAllSeeds(numOfThreads, 0); //TODO add correct identity threshold when none exact seeding is implemented
    }

    int getTargetSeqSize(int seqIdx) const                         { return m_targetUnit.getTargetSeqSize(seqIdx);  } 
    int getQuerySeqSize(int seqIdx) const                          { return m_querySeqs[seqIdx].size();             }
    const string& getTargetSeqName(int seqIdx) const               { return m_targetUnit.getTargetSeqName(seqIdx);  } 
    const string& getQuerySeqName(int seqIdx) const                { return m_querySeqs[seqIdx].Name();             }
    int getNumQuerySeqs() const                                    { return m_querySeqs.getNumSeqs();               } 
    void writeSeeds(const string& overlapFile, int mode) const     { m_seeds.write(overlapFile, mode);              } 


   void alignAllSeqs(int numThreads, ostream& sOut);

    void alignSequence(int querySeqIdx, svec<AlignmentInfo>& cAlignmentInfos) const; 
    void alignSequence(int querySeqIdx, ostream& sOut , ThreadMutex& mtx) const; 

protected:
    const SeedArray& getSeeds(int i) const              { return m_seeds[i]; }
    const AllSeedCandids& getAllSeeds() const           { return m_seeds;    }
 
    void findAllSeeds(int numOfThreads, double identThresh); 
 
    void findSeeds(int querySeqIdx);  
    void findSyntenicBlocks(int querySeqIdx, svec<SyntenicSeeds>& syntBlocks) const;   
    SyntenicSeeds searchDPSynteny(const SeedArray& seeds, int startTIdx, int endTIdx) const; 
    const DNAVector& getTargetSeq(int seqIdx) const { return m_targetUnit.getTargetSeq(seqIdx); }

    void alignSequence(int querySeqIdx, svec<AlignmentInfo>& cAlignmentInfos, int printResults, int storeAlignmentInfo,
                       ostream& sOut, ThreadMutex& mtx) const; 
    void writeAlignment(const Alignment& algn, ostream& sOut, ThreadMutex& mtx) const; 

private:
    DNASeqs                      m_querySeqs;       /// The list of sequences for aligning 
    const FastAlignTargetUnit&   m_targetUnit;      /// An object that handles the target file and creating suffixes from it
    AlignmentParams              m_params;         /// Object containing the various parameters required for assembly
    AllSeedCandids               m_seeds;          /// All candidate seeds among the query/target sequences
    bool                         m_revCmp;         /// Flag is true if this is running for reverse complemented sequences
};

//======================================================



#endif // _ASSEMBLER_H_

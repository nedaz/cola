#ifndef _ASSEMBLER_H_
#define _ASSEMBLER_H_

#include <sstream>
#include "base/ThreadHandler.h"
#include "src/FastAlign/AlignmentParams.h"
#include "src/FastAlign/SuffixArray.h"
#include "src/FastAlign/DNASeqs.h"
#include "src/FastAlign/SeedingObjects.h"
#include "src/FastAlign/SyntenicSeeds.h" 


//======================================================
class FastAlignQueryUnit
{
public:
    FastAlignQueryUnit(const string& inputFile, int stepSize): m_querySeqs(inputFile), 
                                                               m_suffixes(m_querySeqs, stepSize) {}

    const DNAVector&  getQuerySeq(int seqIdx) const { return m_querySeqs[seqIdx];        }
    int getQuerySeqSize(int seqIdx) const           { return m_querySeqs[seqIdx].size(); }

     /** Return a vector of SuffixArrayElement entry indexes for a given 
       string of those Substrings that share a significant subsequence 
       limit specifies the number of overlaps to limit the search to (limit=0 means set the limit to string size) 
     */
    int findSeeds(const DNAVector& targetSeq, int seedSizeThresh, SeedArray& seedArray) const; 

private:
    /** Returns true to indicate that seed has been found and iterating should continue, false otherwise */ 
    template<class IterType>
    bool handleIterInstance(IterType iter, map<unsigned long, int>& stringsUsed_curr, const DNAVector& targetSeq,
                            int targetIterPos, int seedSizeThresh, SeedArray& seedArray) const; 
    int checkInitMatch(const DNAVector& targetSeq, int targetOffset, const SuffixArrayElement& extSeqSA, int seedSizeThresh) const; 


private:
    DNASeqs                          m_querySeqs;      /// A list of sequences from which suffixes were constructed
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
    FastAlignUnit(const string& targetSeqFile, const FastAlignQueryUnit& qUnit, const AlignmentParams& params, int numOfThreads, bool revCmp=false)
                  : m_targetSeqs(targetSeqFile), m_queryUnit(qUnit), 
                    m_params(params), m_seeds(m_targetSeqs.getNumSeqs()), m_revCmp(revCmp) {
        if(revCmp) { m_targetSeqs.reverseComplementAll(); }
        findAllSeeds(numOfThreads, 0); //TODO add correct identity threshold when none exact seeding is implemented
    }

    int getQuerySeqSize(int seqIdx) const                           { return m_queryUnit.getQuerySeqSize(seqIdx);  } 
    int getTargetSeqSize(int seqIdx) const                          { return m_targetSeqs[seqIdx].size();          }
    int getNumTargetSeqs() const                                    { return m_targetSeqs.getNumSeqs();            } 
    void writeSeeds(const string& overlapFile, int mode) const      { m_seeds.write(overlapFile, mode);            } 


   void alignAllSeqs(int numThreads, ostream& sOut);


protected:
    const SeedArray& getSeeds(int i) const              { return m_seeds[i]; }
    const AllSeedCandids& getAllSeeds() const           { return m_seeds;    }
 
    void findAllSeeds(int numOfThreads, double identThresh); 
    void alignSequence(int targetSeqIdx, ostream& sOut, ThreadMutex& mtx) const;
 
    void findSeeds(int targetSeqIdx);  
    void findSyntenicBlocks(int targetSeqIdx, svec<SyntenicSeeds>& syntBlocks) const;   
    SyntenicSeeds searchDPSynteny(const SeedArray& seeds, int startQIdx, int endQIdx) const; 
    const DNAVector& getQuerySeq(int seqIdx) const { return m_queryUnit.getQuerySeq(seqIdx); }

private:
    DNASeqs                     m_targetSeqs;     /// The list of sequences for aligning 
    const FastAlignQueryUnit&   m_queryUnit;      /// An object that handles the query file and creating suffixes from it
    AlignmentParams             m_params;         /// Object containing the various parameters required for assembly
    AllSeedCandids              m_seeds;          /// All candidate seeds among the target/query sequences
    bool                        m_revCmp;         /// Flag is true if this is running for reverse complemented sequences
};

//======================================================



#endif // _ASSEMBLER_H_

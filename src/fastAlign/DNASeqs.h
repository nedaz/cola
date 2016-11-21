#ifndef _DNA_SEQ_H_
#define _DNA_SEQ_H_

#include <string>
#include <stdint.h>
#include "ryggrad/src/general/DNAVector.h"
#include "AlignmentParams.h"

//======================================================
class DNASeqs {

public:
  // Default Ctor:
  DNASeqs(): m_seqs(), m_sizeInfo() {}

  // Ctor 2:
  DNASeqs(const string& fileName)
          : m_seqs(), m_sizeInfo() { 
    load(fileName);           
  }

  const DNAVector& operator[](int i) const              { return m_seqs[i];            }
  const vecDNAVector& getSeqs() const                   { return m_seqs;               }
  const DNAVector& getSeqByIndex(int idx) const         { return m_seqs[idx];          } 
  bool hasSeq(int idx) const                            { return (m_seqs.isize()>idx); }
  int getNumSeqs() const                                { return (m_seqs.size())!=0? m_seqs.size(): m_sizeInfo.size();  } 
  int getSize(int idx) const                            { return ((hasSeq(idx))? m_seqs[idx].size(): m_sizeInfo[idx]);  }

   string getSeqByIndex(int idx, int startIdx, int len) const; 
   string getSeqRCByIndex(int idx, int startIdx, int len) const;

  void reverseComplementAll() {
      m_seqs.ReverseComplement();
  }

  void write(const string& outFile) const; 
  void load(const string& inFile); 

private:
  vecDNAVector m_seqs;            /// Vector containing all sequences 
  svec<int>    m_sizeInfo;        /// The size of each sequence, this is used only when the sequence seqs are not aquired

};

//======================================================

#endif //_DNA_SEQ_H_

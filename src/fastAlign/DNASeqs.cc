#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "DNASeqs.h"

//======================================================

void DNASeqs::write(const string& sequenceFile) const {
    m_seqs.Write(sequenceFile);
}

void DNASeqs::load(const string& sequenceFile) {
    m_seqs.Read(sequenceFile);
}

string DNASeqs::getSeqByIndex(int idx, int startIdx, int len) const {
    return m_seqs[idx].Substring(startIdx, len);
} 

string DNASeqs::getSeqRCByIndex(int idx, int startIdx, int len) const {
    DNAVector rc = m_seqs[idx];
    rc.ReverseComplement();
    return rc.Substring(startIdx, len);
}


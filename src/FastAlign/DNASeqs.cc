#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "src/FastAlign/DNASeqs.h"

//======================================================

void DNASeqs::write(const string& sequenceFile) const {
    m_seqs.Write(sequenceFile);
}

void DNASeqs::load(const string& sequenceFile) {
    m_seqs.Read(sequenceFile);
}


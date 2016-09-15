#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "src/FastAlign/SeedingThreads.h"
#include "src/FastAlign/FastAlignUnit.h"


//======================================================
bool FindSeedsThread::OnDo(const string & msg) {
    int progCount = 0;
    int totSize   = (m_toIdx - m_fromIdx);
    int inc       = totSize/1000;
    if (inc < 1)
        inc = 1;

    for(int i=m_fromIdx; i<m_toIdx; i++) { 
       this->m_alignerUnit.findSeeds(i); 
       progCount++;
       if (progCount % inc == 0) 
           cout << "\r===================== " << 100.0*progCount/totSize 
                << "%  " << flush; 
    }
    return true;
}
//======================================================


//======================================================
bool FindSeedsSingleThread::OnDo(const string & msg) {
    int totSize = this->m_threadQueue.getSize();
    int inc     = totSize/1000;
    if (inc < 1)
        inc = 1;
    int currIdx = this->m_threadQueue.getNext();
    while(currIdx>=0) {
        FILE_LOG(logDEBUG3) << "Finding seeds for read idx: " << currIdx; 
        this->m_alignerUnit.findSeeds(currIdx); 
        if (currIdx  % inc == 0) 
            cout << "\r===================== " << 100.0*currIdx/totSize 
                 << "%  " << flush; 
        currIdx = this->m_threadQueue.getNext();
    }
    return true;
}
//======================================================


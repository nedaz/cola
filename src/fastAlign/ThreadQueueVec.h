#ifndef THREADQUEUEVEC_H
#define THREADQUEUEVEC_H

#include "ryggrad/src/base/SVector.h"
#include "ryggrad/src/base/ThreadHandler.h"

class ThreadQueueVec 
{
public:
    ThreadQueueVec(): m_mutex(), m_totalInstances(0), m_currIndex(0) {}
    ThreadQueueVec(int totSize): m_mutex(), m_totalInstances(totSize), m_currIndex(0) {}

    void setSize(int n) { m_totalInstances = n;    }
    int  getSize()      { return m_totalInstances; }
    void reset()        { m_currIndex = 0;         }
     
    int getNext() {
        int retVal = -1;
        m_mutex.Lock();
        if(m_currIndex < m_totalInstances) {
          retVal = m_currIndex;
          m_currIndex++; 
        }
        m_mutex.Unlock();
        return retVal;
    }

private:
    ThreadMutex  m_mutex;           /// To use for locking while updating current index
    int          m_totalInstances;  /// The total number of instances which is the maximum that currIndex can count to
    int          m_currIndex;       /// The current index that is available for processing
};

#endif //THREADQUEUEVEC_H


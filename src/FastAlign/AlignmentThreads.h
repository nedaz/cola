#ifndef _ALIGNMENT_THREADS_H_
#define _ALIGNMENT_THREADS_H_

#include "src/FastAlign/SeedingThreads.h"


class AlignmentThread : public FindSeedsThread
{
public:
  AlignmentThread(FastAlignUnit& alignerUnit,
                  int from, int to, 
                  int tn, ostream& sOut): FindSeedsThread(alignerUnit, from, to, tn),  
                                          m_sOut(sOut) {}

protected:
  virtual bool OnDo(const string & msg); 

  ostream& m_sOut;                      /// stream for alignment result output
};



class AlignmentSingleThread : public AlignmentThread
{
public:
  AlignmentSingleThread(ThreadQueueVec& tQueue,
                        FastAlignUnit& alignerUnit,
                        int tn, ostream& sOut): AlignmentThread(alignerUnit,  0, 0, tn, sOut),
                                                m_threadQueue(tQueue) {}
protected:
  virtual bool OnDo(const string & msg); 

  ThreadQueueVec& m_threadQueue;        /// Queue for picking up jobs from the queue
};





#endif // _ALIGNMENT_THREADS_H_

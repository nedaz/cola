#ifndef _SEEDING_THREADS_H_
#define _SEEDING_THREADS_H_

#include "base/ThreadHandler.h"
#include "src/FastAlign/ThreadQueueVec.h"
#include "src/FastAlign/SeedingObjects.h"


// Class forward declaration
class FastAlignUnit; 

class FindSeedsThread : public IOneThread
{
public:
  FindSeedsThread(FastAlignUnit& alignerUnit,
                  int from, int to, 
                  int tn): m_alignerUnit(alignerUnit),
                           m_fromIdx(from), m_toIdx(to), 
                           m_threadIdx(tn) {}

protected:
  virtual bool OnDie() { return true; }
  virtual bool OnDo(const string & msg); 
  virtual bool OnInitialize(const string & msg) { return true; }

  FastAlignUnit& m_alignerUnit;       /// Unit that contains the suffixes, and target sequences to be aligned 
  int m_fromIdx;                      /// Index of target sequence to run this thread from
  int m_toIdx;                        /// Index of target sequences to run this thread up to
  int m_threadIdx;                    /// Index of this thread 
};



class FindSeedsSingleThread : public FindSeedsThread
{
public:
  FindSeedsSingleThread(ThreadQueueVec& tQueue,
                        FastAlignUnit& alignerUnit,
                        int tn): FindSeedsThread(alignerUnit,  0, 0, tn),
                                 m_threadQueue(tQueue) {}
protected:
  virtual bool OnDo(const string & msg); 

  ThreadQueueVec& m_threadQueue;      /// Queue for picking up jobs from the queue
};





#endif // _SEEDING_THREADS_H_

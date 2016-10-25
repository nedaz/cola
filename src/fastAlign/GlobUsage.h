#ifndef GLOBUSAGE_H
#define GLOBUSAGE_H

#include "CachedVec.h"
#include "ryggrad/src/base/SVector.h"
#include "ryggrad/src/base/ThreadHandler.h"

class GlobalUsageHandler
{
 public:
  GlobalUsageHandler() {
    m_counter = 0;
  }
  void SetSize(int n) {
    m_list.resize(n, 0);
  }

  void Sync(VecIntInc & s);

 private:
  VecIntInc m_list;
  ThreadMutex m_mutex;
  int m_counter;
};

#endif //GLOBUSAGE_H


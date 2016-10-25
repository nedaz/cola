#include "GlobUsage.h"
#include "ryggrad/src/base/ThreadHandler.h"

void GlobalUsageHandler::Sync(VecIntInc & s)
{
    m_mutex.Lock();
    //cout << "Reporting." << endl;
    if (m_list.isize() == 0) {
        m_list.resize(s.isize());
    }
    if (s.isize() != m_list.isize()) {
        cout << "GLOBAL USAGE/Mutlithreading ERROR. " << endl;
        m_mutex.Unlock();
        throw;
    }
    int i;

    int curr = m_list.NumTotal();
    int from = s.From();
  
    //cout << "Absorbing: " << s.NumMod() << endl;
    for (i=0; i<s.NumMod(); i++) {
        int m = s.GetMod(i);
        m_list.Set(m, s[m]);
    }
    // cout << "Successful: " << m_list.NumTotal() - curr  << endl;
  
    //cout << "Propagating: " << curr - from << endl;
    int pre = s.NumMod();
    for (i=from; i<curr; i++) {
        int n = m_list.GetMod(i);
        s.Set(n, m_list[n]);
    }
    //cout << "Successful: " << s.NumMod() - pre  << endl;
    s.Reset();
    //double r = (double)m_list.NumTotal() / (double)m_list.isize()); 
  
    int now = m_list.NumTotal();
    int inc = m_list.isize() / 1000;
    if (now > m_counter + inc) {
        m_counter = now;
        cout << "\rReads processed: " << m_counter << " " << 100.*(double)m_counter/(double)m_list.isize() << " %" << flush;
    }

    //CHECK HERE
    //for (i=0; i<m_list.isize(); i++) {
    //  if (m_list[i] != s[i]) {
    //    cout << "ERROR: " << i << " " << m_list[i] << " - " << s[i] << endl;
    //  }
    //}

    m_mutex.Unlock();

}


/*void GlobalUsageHandler::Sync(VecIntInc & s)
  {
  m_mutex.Lock();
  
  if (m_list.isize() == 0) {
  m_list.resize(s.isize());
  for (int i=0; i<s.isize(); i++)
  m_list[i] = s[i];
  m_mutex.Unlock();
  return;
  }
  if (s.isize() != m_list.isize()) {
  cout << "GLOBAL USAGE/Mutlithreading ERROR. " << endl;
  m_mutex.Unlock();
  throw;
  }
  int i;
  int occ = 0;
  int rep = 0;
  int back = 0;
  for (i=0; i<m_list.isize(); i++) {
  int a = m_list[i];
  int b = s[i];

  if (b > 0 && a == 0)
  rep++;
  if (a > 0 && b == 0)
  back++;
  if (a > 0)
  occ++;

  if (a > b)
  b = a;
  m_list[i] = b;
  if (s[i] != b)
  s.Set(i, b);
  }
  cout << "Sync'ed. Total: " << occ << ", new: " << rep << ", comm: " << back << endl;
  m_mutex.Unlock();

  }
*/

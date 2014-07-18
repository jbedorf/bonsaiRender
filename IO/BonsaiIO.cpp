#include "BonsaiIO.h"

namespace BonsaiIO
{

  template<typename T>
  void Core::addType(const std::vector<T> &data)
  {
    const uint64_t n = data.size();
    assert(n > 0);
    DataStructBase *ptr = new T[n];
#pragma omp parallel for
    for (int i = 0; i < n; i++)
      ptr[i] = data[i];

    list.push_back(std::make_pair(n, ptr);
  }

  void Core::writeData(const std::string &fileName)
  {
    const int nTypes = 
    if (isMaster())
    {
    }

  }


  void IDType_t :: write(
}

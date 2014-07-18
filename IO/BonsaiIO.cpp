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

    data.push_back(std::make_pair(n, ptr));
  }

  void Core::writeData(const std::string &fileName)
  {
    const int nTypes = data.size();
    std::vector<uint64_t> nBytesLoc(nTypes), nBytesGlb(nTypes);
    for (int i = 0; i < nTypes; i++)
      nBytesLoc[i] = data[i].first*data[i]->getSize();

    MPI_Allreduce(&nBytesLoc[0], &nBytesGlb[0], nTypes, MPI_LONG_LONG, MPI_SUM, MPIComm);

    std::vector<std::pair<uint64_t,uint64_t>> begEnd(nTypes);

    for (int i = 0; i < nTypes; i++)
    {
      const uint64_t nBytesPerRank = (nBytesGlb[i] - 1 + nRank)/nRank;
      const uint64_t beg = myRank  * nBytesPerRank;
      const uint64_t end = min(beg + nBytesPerRank, nBytesGlb[i]);
      begEnd[i] = std::make_pair(beg,end);
    }

    MPI_File fh;
    MPI_Status status;
  MPI_Offset my_offset, my_current_offset;
    MPI_File_open(
        MPIComm, 
        fileName.c_str(),
        MPI_MODE_CREATE | MPI_MODE_RDWR,
        MPI_INFO_NULL,
        &fn);
    MPI_File_see(fh, beg
      

  }


  void IDType_t :: write(
}

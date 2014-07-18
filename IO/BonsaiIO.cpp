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
    for (inti = 0; i < nTypes; i++)
      nBytesLoc[i] = data[i].first*data[i]->getSize();

    MPI_Allreduce(&nBytesLoc[0], &nBytesGlb[0], nTypes, MPI_LONG_LONG, MPI_SUM, MPIComm);

    std::vector<std::pair<uint64_t,uint64_t>> begEnd(nTypes);


    uint64_t offset = 0;
    for (int i = 0; i < nTypes; i++)
    {
      const uint64_t nBytesPerRank = (nBytesGlb[i] - 1 + nRank)/nRank;
      const uint64_t beg = myRank  * nBytesPerRank;
      const uint64_t end = min(beg + nBytesPerRank, nBytesGlb[i]);
      begEnd[i] = std::make_pair(offset + beg,offset + end);
      offset += nBytesGlb[i];
    }
    const nBytes2Write = offset;

    MPI_File fh;
    MPI_Status status;

    const double tWrite = MPI_Wtime();
    const double tOpen  = MPI_Wtime();
    MPI_File_open(
        MPIComm, 
        fileName.c_str(),
        MPI_MODE_CREATE | MPI_MODE_RDWR,
        MPI_INFO_NULL,
        &fn);
    /* write header */
//    if (isMaster())
    
    const double dtOpen = MPI_Wtime() - tOpen;
    for (int type = 0; type < nTypes; type++)
    {
      MPI_Offset myOffset = begEnd[type].first;
      uint64_t nBytes = begEnd[type].second - myOffset;
      const uint64_t nMaxPerWrite = (1U << 31) - 1;
      while (nBytes > 0)
      {
        const int count = min(nBytes, nMaxPerWrite);
        MPI_File_seek(fh, myOffset, MPI_SEEK_SET);
        MPI_File_write(
            fh,
            data[type].second, 
            nMaxPerWrite, 
            MPI_BYTE, 
            &status);
        int written;
        MPI_Get_count(&status, MPI_INT, &written);
        assert(count == written);
        nBytes -= nMaxPerWrite;
      }
    }
    MPI_File_close(&fh);
    const double dtWrite = MPI_Wtime() - tWrite;

    double dtWriteMax, dtOpenMax;
    MPI_Allreduce(&dtWrite, &dtWriteMax, 1, MPI_DOUBLE, MPI_MAX, MPIComm);
    MPI_Allreduce(&dtOpen,  &dtOpenMax,  1, MPI_DOUBLE, MPI_MAX, MPIComm);

    if (isMaster())
    {
      fprintf(stderr, " wrote %g MB of data in %g sec [BW= %g MB/s] open= %g sec \n",
          static_cast<double>(nBytes2Write)/1e6, dtWriteMax,
          static_cast<double>(nBytes2Write)/1e6 / dtWriteMax,
          dtOpenMax);
    }
  }


  void IDType_t :: write(
}

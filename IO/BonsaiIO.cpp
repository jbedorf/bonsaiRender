#include "BonsaiIO.h"
  
class IDType;
{
  public:
    uint64_t ID;
  private:
    uint64_t getID() const
    {
      return _IDTypePacked & 0xFFFF000000000000ULL;
    }
    uint32_t getType() const
    {
      return static_cast<uint32_t>(_IDTypePacked >> 48);
    }
    void setID(const int64_t ID)
    {
      const uint32_t type = getType();
      _IDTypePacked = (ID & 0xFFFF000000000000ULL) | (static_cast<uint64_t>(type) << 48);
    }
    void setType(const int type)
    {
      const uint64_t ID = getID();
      _IDTypePacked  = ID | (static_cast<uint64_t>(type) << 48);
    }
};

void writeSnapshot(
    real4 *bodyPositions,
    real4 *bodyVelocities,
    uulong* bodiesIDs,
    const int n,
    const int nDomains,
    const std::string &fileName,
    const float time,
    const MPIComm &comm,
    const int nRank, const int myRank)
{
  BonsaiIO::Core out(myRank, nRank, comm, BonsaiIO::WRITEMPI, fileName);

  out.setNDomains(nDomains);

  /* write IDs */
  {
    BonsaiIO::DataType<IDType> ID("IDType", n);
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
      ID[i].setID(bodiesIDs[i]);
      int type = 0;
      if(bodyIds[i] >= DISKID  && bodyIds[i] < BULGEID)       type = 2;
      if(bodyIds[i] >= BULGEID && bodyIds[i] < DARKMATTERID)  type = 1;
      if(bodyIds[i] >= DARKMATTERID)                          type = 0;
      ID[i].setType(type);
    }
    out.write(ID);
  }

  /* write pos */
  {
    BonsaiIO::DataType<real4> pos("POS",n);
#pragma omp parallel for
    for (int i = 0; i < n; i++)
      pos[i] = bodyPositions[i];
    out.write(pos);
  }

  /* write velocities */
  {
    typedef float[3] vec3;
    BonsaiIO::Vector<vec3> vel("VEL",n);
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
      vel[i][0] = bodyVelocities[i].x;
      vel[i][1] = bodyVelocities[i].y;
      vel[i][2] = bodyVelocities[i].z;
    }
    out.write(vel);
  }
  
  /* write rhoh */
  {
    typedef float[2] vec2;
    BonsaiIO::DataType<vec2> rhoh("RHOH",n);
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
      rhoh[i][0] = 0; /* rho */
      rhoh[i][1] = 0; /*  h  */
    }
    out.write(rhoh);
  }

  out.close();
}

void readSnapshot(
    std::vector<real4>  &bodyPositions,
    std::vector<real4>  &bodyVelocities,
    std::vector<uulong> &bodiesID,
    std::vector<real2>  &rhohList,
    const std::string   &fileName,
    const MPIComm       &comm,
    const int nRank, 
    const int myRank,
    int &NTotal2,
    int &NFirst, int &NSecond, int &NThird,
    std::vector<real4> &dustPositions, std::vector<real4> &dustVelocities,
    std::vector<ullong> &dustIDs, 
    const int reduce_bodies_factor,
    const int reduce_dust_factor,
    const bool restart)
{
  NFirst = NSecond = NThird = 0;
  BonsaiIO::Core out(fileName, BonsaiIO::READMPI);

  {
    BonsaiIO::DataType<IDType> IDList("IDType");
    if (!out.read(IDList))
    {
      if (myRank == 0)
        fprintf(stderr, " FATAL: No particle ID data is found. Please make sure you passed the right file \n");
      exit(-1);
    }
    const int n = IDList.size();
    bodiesID.resize(n);
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
      bodiesID[i] = IDList[i].getID();
      switch (IDList[i].getType())
      {
        case 0:
          NFirst++;
          break;
        case 1:
          NSecond++;
          break;
        case 2:
          NThird++;
          break;
      }
    }
    NTotal2 = NFirst+NSecond+NThird;
    assert(NTotal2 == n);
  }

  {
    BonsaiIO::DataType<real4> pos("pos");
    if (!out.read(pos))
    {
      if (myRank == 0)
        fprintf(stderr, " FATAL: No particle positions data is found. Please make sure you passed the right file \n");
      exit(-1);
    }
    const int n = pos.size();
    bodyPositions.resize(n);
    assert(bodyPositions.size() == bodiesID.size());
#pragma omp parallel for
    for (int i = 0; i < n; i++)
      bodyPositions[i] = pos[i];
  }

  {
    typedef float[3] vec3;
    Bonsai::DataType<vec3> vel("VEL");
    if (!out.read(vel))
    {
      if (myRank == 0)
        fprintf(stderr, " FATAL: No particle velocity data is found. Please make sure you passed the right file \n");
      exit(-1);
    }
    const int n = vel.size();
    bodyVelocities.resize(n);
    assert(bodyVelocities.size() == bodiesID.size());
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
      bodyVelocities[i].x = vel[i][0];
      bodyVelocities[i].y = vel[i][1];
      bodyVelocities[i].z = vel[i][2];
    }
  }

  {
    typedef float[2] vec2;
    Bonsai::DataType<vec2> rhoh("RHOH");
    if (out.read.(rhoh))
    {
      if (myRank == 0)
      {
        fprintf(stderr , " -- No RHOH data is found \n");
      }
      const int n = rhoh.size();
      rhohList.resize(n);
      assert(rhohList.size() == bodiesID.size());
#pragma omp parallel for
      for (int i = 0; i < n; i++)
      {
        rhohList[i].x = rhoh[i][0];
        rhohList[i].y = rhoh[i][1];
      }
    }
    else if (myRank == 0)
    {
      fprintf(stderr , " -- No RHOH data is found \n");
    }
  }

}
                                                                                                                                                                                                                        const bool restart)



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
}

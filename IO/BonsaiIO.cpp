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
    typedef float vec3[3];
    BonsaiIO::DataType<vec3> vel("VEL",n);
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
    typedef float vec2[2];
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
    typedef float vec3[3];
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
    typedef float vec2[2];
    Bonsai::DataType<vec2> rhoh("RHOH");
    if (out.read.(rhoh))
    {
      if (myRank == 0)
        fprintf(stderr , " -- RHOH data is found \n");
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
      fprintf(stderr , " -- No RHOH data is found \n");
  }

}




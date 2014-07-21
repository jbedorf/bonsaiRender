#include "BonsaiIO.h"
#include "IDType.h"


int main(int argc, char * argv[])
{
#if 0
  IDType d;
  d.setType(12);
  d.setID(345462);

  fprintf(stderr, " id= %lu  type= %d\n",
      d.getID(), d.getType());
  return 0 ;
#endif
  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Init(&argc, &argv);
    
  int nRank, myRank;
  MPI_Comm_size(comm, &nRank);
  MPI_Comm_rank(comm, &myRank);


  if (argc < 2)
  {
    if (myRank == 0)
    {
      fprintf(stderr, " ------------------------------------------------------------------------\n");
      fprintf(stderr, " Usage: \n");
      fprintf(stderr, " %s  fileName [reduceFactor]\n", argv[0]);
      fprintf(stderr, " ------------------------------------------------------------------------\n");
    }
    exit(-1);
  }
  
  const std::string fileName(argv[1]);

  int reduceFactor = 1;
  if (argc > 2)
  {
    reduceFactor = atoi(argv[2]);
    fprintf(stderr," reduceFactor= %d\n", reduceFactor);
  }
  
  BonsaiIO::Core out(myRank, nRank, comm, BonsaiIO::READ, fileName);

  if (myRank == 0)
    out.getHeader().printFields();

  std::vector<long long int> bodyID;
  long long int NFirst = 0, NSecond = 0, NThird = 0, NTotal2 = 0;
  
  {
    BonsaiIO::DataType<IDType> IDList("IDType");
    if (!out.read(IDList))
    {
      if (myRank == 0)
        fprintf(stderr, " FATAL: No particle ID data is found. Please make sure you passed the right file \n");
      exit(-1);
    }
    const int n = IDList.getNumElements();
    bodyID.resize(n);
    long long int NFirstL = 0, NSecondL = 0, NThirdL = 0;
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
      bodyID[i] = IDList[i].getID();
      switch (IDList[i].getType())
      {
        case 0:
          NFirstL++;
          break;
        case 1:
          NSecondL++;
          break;
        case 2:
          NThirdL++;
          break;
        default:
          fprintf(stderr, "id= %lu\n", IDList[i].getID());
          fprintf(stderr, "type= %d\n", IDList[i].getType());
          assert(0);
      }
    }
    MPI_Allreduce(&NFirstL,  &NFirst, 1, MPI_LONG_LONG, MPI_SUM, comm);
    MPI_Allreduce(&NSecondL, &NSecond, 1, MPI_LONG_LONG, MPI_SUM, comm);
    MPI_Allreduce(&NThirdL,  &NThird, 1, MPI_LONG_LONG, MPI_SUM, comm);
    long long NTotal2L = NFirstL+NSecondL+NThirdL;
    NTotal2  = NFirst +NSecond +NThird ;
    if (myRank == 0)
      fprintf(stderr, "rank= %d  n= %d  nTotal2L= %lld  nTotal2= %lld\n",
          myRank, n, NTotal2L, NTotal2);
    assert(NTotal2L == n);
  }
 
  struct real4 {float x,y,z,w;} ;
  std::vector<real4> bodyPositions;
  {
    BonsaiIO::DataType<real4> pos("POS:real4");
    if (!out.read(pos))
    {
      if (myRank == 0)
        fprintf(stderr, " FATAL: No particle positions data is found. Please make sure you passed the right file \n");
      exit(-1);
    }
    const int n = pos.getNumElements();
    bodyPositions.resize(n);
    assert(bodyPositions.size() == bodyID.size());
#pragma omp parallel for
    for (int i = 0; i < n; i++)
      bodyPositions[i] = pos[i];
  }
 
  std::vector<real4> bodyVelocities; 
  {
    typedef float vec3[3];
    BonsaiIO::DataType<vec3> vel("VEL:float[3]");
    if (!out.read(vel))
    {
      if (myRank == 0)
        fprintf(stderr, " FATAL: No particle velocity data is found. Please make sure you passed the right file \n");
      exit(-1);
    }
    const int n = vel.getNumElements();
    bodyVelocities.resize(n);
    assert(bodyVelocities.size() == bodyID.size());
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
      bodyVelocities[i].x = vel[i][0];
      bodyVelocities[i].y = vel[i][1];
      bodyVelocities[i].z = vel[i][2];
      bodyVelocities[i].w = 0;
    }
  }

  struct real2 {float x,y;};
  std::vector<real2> rhohList;
  {
    typedef float vec2[2];
    BonsaiIO::DataType<vec2> rhoh("RHOH");
    if (out.read(rhoh))
    {
      if (myRank == 0)
        fprintf(stderr , " -- RHOH data is found \n");
      const int n = rhoh.getNumElements();
      rhohList.resize(n);
      assert(rhohList.size() == bodyID.size());
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

  if (myRank == 0)
    fprintf(stderr, " Bandwidth= %g MB/s\n", out.computeBandwidth()/1e6);

  out.close();




#if 0
  {
    BonsaiIO::Core out(rank, nranks, comm, BonsaiIO::WRITE, outputName);

    /* write IDs */
    {
      BonsaiIO::DataType<IDType> ID("IDType", nTotalLoc);
#pragma omp parallel for
      for (int i = 0; i< nFirstLocal; i++)
      {
        ID[i].setID(data.firstID[i]);
        ID[i].setType(0);
      }
#pragma omp parallel for
      for (int i = 0; i< nSecondLocal; i++)
      {
        ID[i+nFirstLocal].setID(data.secondID[i]);
        ID[i+nFirstLocal].setType(1);
      }
      out.write(ID);
    }
  
    /* write pos */
    {
      BonsaiIO::DataType<ReadTipsy::real4> pos("POS:real4",nTotalLoc);
#pragma omp parallel for
      for (int i = 0; i< nFirstLocal; i++)
        pos[i] = data.firstPos[i];
#pragma omp parallel for
      for (int i = 0; i< nSecondLocal; i++)
        pos[i+nFirstLocal] = data.secondPos[i];
      out.write(pos);
    }
    
    /* write vel */
    {
      typedef float vec3[3];
      BonsaiIO::DataType<vec3> vel("VEL:float[3]",nTotalLoc);
#pragma omp parallel for
      for (int i = 0; i< nFirstLocal; i++)
      {
        vel[i][0] = data.firstVel[i].x;
        vel[i][1] = data.firstVel[i].y;
        vel[i][2] = data.firstVel[i].z;
      }
#pragma omp parallel for
      for (int i = 0; i< nSecondLocal; i++)
      {
        vel[i+nFirstLocal][0] = data.secondVel[i].x;
        vel[i+nFirstLocal][1] = data.secondVel[i].y;
        vel[i+nFirstLocal][2] = data.secondVel[i].z;
      }
      out.write(vel);
    }

    out.close();

    if (rank == 0)
      fprintf(stderr, " Bandwidth= %g MB/s\n", out.computeBandwidth()/1e6);
  }
#endif




  MPI_Finalize();




  return 0;
}



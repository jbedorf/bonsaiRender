#include "BonsaiIO.h"
#include "IDType.h"
#include "read_tipsy.h"


int main(int argc, char * argv[])
{
  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Init(&argc, &argv);
    
  int nranks, rank;
  MPI_Comm_size(comm, &nranks);
  MPI_Comm_rank(comm, &rank);


  if (argc < 4)
  {
    if (rank == 0)
    {
      fprintf(stderr, " ------------------------------------------------------------------------\n");
      fprintf(stderr, " Usage: \n");
      fprintf(stderr, " %s  baseName nDomains outputName \n", argv[0]);
      fprintf(stderr, " ------------------------------------------------------------------------\n");
    }
    exit(-1);
  }
  
  const std::string baseName(argv[1]);
  const int nDomains = atoi(argv[2]);
  const std::string outputName(argv[3]);

  int reduceFactorFirst  = 1;
  int reduceFactorSecond = 1;


  ReadTipsy data(
      baseName, 
      rank, nranks,
      nDomains, 
      reduceFactorFirst,
      reduceFactorSecond);

  long long nFirstLocal = data.firstID.size();
  long long nSecondLocal = data.secondID.size();
  long long nTotalLoc = nFirstLocal + nSecondLocal;

  long long nFirst, nSecond;
  MPI_Allreduce(&nFirstLocal, &nFirst, 1, MPI_LONG, MPI_SUM, comm);
  MPI_Allreduce(&nSecondLocal, &nSecond, 1, MPI_LONG, MPI_SUM, comm);

  if (rank == 0)
  {
    fprintf(stderr, " nFirst = %lld \n", nFirst);
    fprintf(stderr, " nSecond= %lld \n", nSecond);
    fprintf(stderr, " nTotal= %lld \n", nFirst + nSecond);
  }

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
    {
      out.getHeader().printFields();
      fprintf(stderr, " Bandwidth= %g MB/s\n", out.computeBandwidth()/1e6);
    }
  }




  MPI_Finalize();




  return 0;
}



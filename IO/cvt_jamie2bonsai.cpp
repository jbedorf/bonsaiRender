#include "BonsaiIO.h"

#if 0
template<typename IO>
static double writeDM(ReadTipsy &data, IO &out)
{
  double dtWrite = 0;
  const int pCount  = data.firstID.size();
  /* write IDs */
  {
    BonsaiIO::DataType<IDType> ID("DM:IDType", pCount);
    for (int i = 0; i< pCount; i++)
    {
      ID[i].setID(data.firstID[i]);
      ID[i].setType(0);
    }
    double t0 = MPI_Wtime();
    out.write(ID);
    dtWrite += MPI_Wtime() - t0;
  }
  
  /* write pos */
  {
    BonsaiIO::DataType<ReadTipsy::real4> pos("DM:POS:real4",pCount);
    for (int i = 0; i< pCount; i++)
      pos[i] = data.firstPos[i];
    double t0 = MPI_Wtime();
    out.write(pos);
    dtWrite += MPI_Wtime() - t0;
  }
    
  /* write vel */
  {
    typedef float vec3[3];
    BonsaiIO::DataType<vec3> vel("DM:VEL:float[3]",pCount);
    for (int i = 0; i< pCount; i++)
    {
      vel[i][0] = data.firstVel[i].x;
      vel[i][1] = data.firstVel[i].y;
      vel[i][2] = data.firstVel[i].z;
    }
    double t0 = MPI_Wtime();
    out.write(vel);
    dtWrite += MPI_Wtime() - t0;
  }

  return dtWrite;
}

template<typename IO>
static double writeStars(ReadTipsy &data, IO &out)
{
  double dtWrite = 0;

  const int pCount  = data.secondID.size();

  /* write IDs */
  {
    BonsaiIO::DataType<IDType> ID("Stars:IDType", pCount);
    for (int i = 0; i< pCount; i++)
    {
      ID[i].setID(data.secondID[i]);
      ID[i].setType(1);
    }
    double t0 = MPI_Wtime();
    out.write(ID);
    dtWrite += MPI_Wtime() - t0;
  }
    
  /* write pos */
  {
    BonsaiIO::DataType<ReadTipsy::real4> pos("Stars:POS:real4",pCount);
    for (int i = 0; i< pCount; i++)
      pos[i] = data.secondPos[i];
    double t0 = MPI_Wtime();
    out.write(pos);
    dtWrite += MPI_Wtime() - t0;
  }
    
  /* write vel */
  {
    typedef float vec3[3];
    BonsaiIO::DataType<vec3> vel("Stars:VEL:float[3]",pCount);
    for (int i = 0; i< pCount; i++)
    {
      vel[i][0] = data.secondVel[i].x;
      vel[i][1] = data.secondVel[i].y;
      vel[i][2] = data.secondVel[i].z;
    }
    double t0 = MPI_Wtime();
    out.write(vel);
    dtWrite += MPI_Wtime() - t0;
  }

  return dtWrite;
}
#endif


int main(int argc, char * argv[])
{
  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Init(&argc, &argv);
    
  int nranks, rank;
  MPI_Comm_size(comm, &nranks);
  MPI_Comm_rank(comm, &rank);


  if (argc < 3)
  {
    if (rank == 0)
    {
      fprintf(stderr, " ------------------------------------------------------------------------\n");
      fprintf(stderr, " Usage: \n");
      fprintf(stderr, " %s  inputFileName outputFileName \n", argv[0]);
      fprintf(stderr, " ------------------------------------------------------------------------\n");
    }
    exit(-1);
  }
  
  const std::string inputFn(argv[1]);
  const std::string outputFn(argv[2]);

  {
    FILE *fin = fopen(inputFn.c_str(), "rb");
    assert(fin);

    int ntot, idum;
    double t, rdum;
    fread(&idum /* recsz */,  sizeof(int), 1, fin);
    fread(&ntot /* ntot  */,  sizeof(int), 1, fin);
    fread(&idum /* nnopt */, sizeof(int), 1, fin);
    fread(&rdum /* hmin  */, sizeof(double), 1, fin);
    fread(&rdum /* hmax  */, sizeof(double), 1, fin);
    fread(&rdum /* sep0  */, sizeof(double), 1, fin);
    fread(&rdum /* tf    */, sizeof(double), 1, fin);
    fread(&rdum /* dtout */, sizeof(double), 1, fin);
    fread(&idum /* nout  */, sizeof(int), 1, fin);
    fread(&idum /* nit   */, sizeof(int), 1, fin);
    fread(&t     /* t     */, sizeof(double), 1, fin);
    fread(&idum /* nav   */, sizeof(int), 1, fin);
    fread(&rdum /* alpha */, sizeof(double), 1, fin);
    fread(&rdum /* beta  */, sizeof(double), 1, fin);
    fread(&rdum /* tskip */, sizeof(double), 1, fin);
    fread(&idum /* ngr   */, sizeof(int), 1, fin);
    fread(&idum /* nrelax*/, sizeof(int), 1, fin);
    fread(&rdum /* trelax*/, sizeof(double), 1, fin);
    fread(&rdum /* dt    */, sizeof(double), 1, fin);
    fread(&rdum /* omega2*/, sizeof(double), 1, fin);

    fprintf(stderr, "ntot =%d  t= %g \n", ntot, t);
  }


#if 0
  ReadTipsy data(
      baseName, 
      rank, nranks,
      nDomains, 
      reduceFactorFirst,
      reduceFactorSecond);

  long long nFirstLocal = data.firstID.size();
  long long nSecondLocal = data.secondID.size();

  long long nFirst, nSecond;
  MPI_Allreduce(&nFirstLocal, &nFirst, 1, MPI_LONG, MPI_SUM, comm);
  MPI_Allreduce(&nSecondLocal, &nSecond, 1, MPI_LONG, MPI_SUM, comm);

  if (rank == 0)
  {
    fprintf(stderr, " nFirst = %lld \n", nFirst);
    fprintf(stderr, " nSecond= %lld \n", nSecond);
    fprintf(stderr, " nTotal= %lld \n", nFirst + nSecond);
  }

  const double tAll = MPI_Wtime();
  {
    const double tOpen = MPI_Wtime();
    BonsaiIO::Core out(rank, nranks, comm, BonsaiIO::WRITE, outputName);
    double dtOpenLoc = MPI_Wtime() - tOpen;
    double dtOpenGlb;
    MPI_Allreduce(&dtOpenLoc, &dtOpenGlb, 1, MPI_DOUBLE, MPI_MAX,comm);
    if (rank == 0)
      fprintf(stderr, "open file in %g sec \n", dtOpenGlb);

    double dtWrite = 0;

    
    if (rank == 0)
      fprintf(stderr, " write DM  \n");
    MPI_Barrier(comm);
    dtWrite += writeDM(data,out);

    if (rank == 0)
      fprintf(stderr, " write Stars\n");
    MPI_Barrier(comm);
    dtWrite += writeStars(data,out);


    double dtWriteGlb;
    MPI_Allreduce(&dtWrite, &dtWriteGlb, 1, MPI_DOUBLE, MPI_MAX,comm);
    if (rank == 0)
      fprintf(stderr, "write file in %g sec \n", dtWriteGlb);

    const double tClose = MPI_Wtime();
    out.close();
    double dtCloseLoc = MPI_Wtime() - tClose;
    double dtCloseGlb;
    MPI_Allreduce(&dtCloseLoc, &dtCloseGlb, 1, MPI_DOUBLE, MPI_MAX,comm);
    if (rank == 0)
      fprintf(stderr, "close time in %g sec \n", dtCloseGlb);

    if (rank == 0)
    {
      out.getHeader().printFields();
      fprintf(stderr, " Bandwidth= %g MB/s\n", out.computeBandwidth()/1e6);
    }
  }
  double dtAllLoc = MPI_Wtime() - tAll;
  double dtAllGlb;
  MPI_Allreduce(&dtAllLoc, &dtAllGlb, 1, MPI_DOUBLE, MPI_MAX,comm);
  if (rank == 0)
    fprintf(stderr, "All operations done in   %g sec \n", dtAllGlb);

#endif



  MPI_Finalize();




  return 0;
}



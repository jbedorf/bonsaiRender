#include "BonsaiIO.h"

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

  FILE *fin = fopen(inputFn.c_str(), "rb");
  assert(fin);

  struct __attribute__((__packed__)) header_t
  {
    int ntot;
    int nnopt;
    double hmin;
    double hmax;
    double sep0;
    double tf;
    double dtout;
    int nout;
    int nit;
    double t;
    int anv;
    double alpha;
    double beta;
    double tskip;
    int ngr;
    int nrelax;
    double trelax;
    double dt;
    double omega2;
  };

  int idum; 
  header_t h;
  fread(&idum,  sizeof(int), 1, fin); 
  assert(idum == (int)sizeof(header_t));

  fread(&h, sizeof(header_t), 1, fin);

  fread(&idum,  sizeof(int), 1, fin); 
  assert(idum == (int)sizeof(header_t));

  fprintf(stderr, "ntot= %d  t= %g \n", h.ntot, h.t);

  struct __attribute__((__packed__)) sph_t
  {
    double x,y,z;
    double am,hp,rho;
    double vx,vy,vz;
    double vxdot,vydot,vzdot;
    double u,udot;
    double grpot, mmu;
    int cc;
    double divv;
  };

  std::vector<sph_t> data(h.ntot);

  for (int i = 0; i < h.ntot; i++)
  {
    fread(&idum,  sizeof(int), 1, fin); 
    assert(idum == (int)sizeof(sph_t));

    fread(&data[i], sizeof(sph_t), 1, fin);

    fread(&idum,  sizeof(int), 1, fin); 
    assert(idum == (int)sizeof(sph_t));
  }

  fread(&idum,  sizeof(int), 1, fin); 
  assert(idum == (int)sizeof(int));

  fread(&idum, sizeof(int), 1, fin);
  assert(idum == h.ntot);

  fread(&idum,  sizeof(int), 1, fin); 
  assert(idum == (int)sizeof(int));
   

  const int ntot = h.ntot;
  BonsaiIO::Core out(rank, nranks, comm, BonsaiIO::WRITE, outputFn);

  /* write metadata */
  {
    BonsaiIO::DataType<header_t> el("SPH:header:jamieHeader_t", 1);
    el[0] = h;
    out.write(el);
  }

  /* write pos */
  {
    BonsaiIO::DataType<sph_t> el("SPH:data:jamieData_t", ntot);
    for (int i = 0; i < ntot; i++)
      el[i] = data[i];
    out.write(el);
  }

  out.close();


  MPI_Finalize();




  return 0;
}



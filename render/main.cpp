#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <unistd.h>
#include <mpi.h>
#include <sstream>
#include "read_tipsy.h"

#include "renderloop.h"
#include "anyoption.h"


int main(int argc, char * argv[])
{
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);
    
  int nranks, rank;
  MPI_Comm_size(comm, &nranks);
  MPI_Comm_rank(comm, &rank);

  assert(nranks == 1);
  assert(rank   == 0);

  std::string fileName;
  int nDomains     = -1;
  int reduceFactor =  1;
  std::string fullScreenMode    = "";
  bool stereo     = false;


  if (rank == 0)  
  {
		AnyOption opt;

#define ADDUSAGE(line) {{std::stringstream oss; oss << line; opt.addUsage(oss.str());}}

		ADDUSAGE(" ");
		ADDUSAGE("Usage:");
		ADDUSAGE(" ");
		ADDUSAGE(" -h  --help             Prints this help ");
		ADDUSAGE(" -i  --infile #         Input snapshot filename ");
    ADDUSAGE(" -n  --ndomains #       Number of domains ");
		ADDUSAGE(" -r  --reducefactor #   cut down bodies dataset by # factor [1])")
		ADDUSAGE("     --fullscreen #     set fullscreen mode string");
		ADDUSAGE("     --stereo           enable stereo rendering");
		ADDUSAGE(" ");


		opt.setFlag  ( "help" ,        'h');
		opt.setOption( "infile",       'i');
		opt.setOption( "ndomains",     'n');
		opt.setOption( "reducefactor", 'r');
    opt.setOption( "fullscreen");
    opt.setFlag("stereo");

    opt.processCommandArgs( argc, argv );


    if( ! opt.hasOptions() ||  opt.getFlag( "help" ) || opt.getFlag( 'h' ) )
    {
      /* print usage if no options or requested help */
      opt.printUsage();
      ::exit(0);
    }

    char *optarg = NULL;
    if ((optarg = opt.getValue("infile")))       fileName           = std::string(optarg);
    if ((optarg = opt.getValue("ndomains")))     nDomains           = atoi(optarg);
    if ((optarg = opt.getValue("reducefactor"))) reduceFactor       = atoi(optarg);
    if ((optarg = opt.getValue("fullscreen")))	 fullScreenMode     = std::string(optarg);
    if (opt.getFlag("stereo"))  stereo = true;

    if (fileName.empty() || nDomains == -1 || reduceFactor < 1)
    {
      opt.printUsage();
      ::exit(0);
    }

#undef ADDUSAGE
  }


  ReadTipsy data(
      fileName, 
      rank, nranks,
      nDomains, 
      reduceFactor,
      reduceFactor);

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

  initGL(argc, argv, fullScreenMode.c_str(), stereo);  
  fprintf(stderr, " -- Done -- \n");
  while(1) {}
  return 0;
}



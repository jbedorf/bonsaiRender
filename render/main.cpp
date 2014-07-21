#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <unistd.h>
#include <mpi.h>
#include <sstream>
#include <cmath>
#if 0
#include "read_tipsy.h"
#define OLDIO
#else
#include "IDType.h"
#include "BonsaiIO.h"
#endif

#include "renderloop.h"
#include "anyoption.h"
#include "RendererData.h"



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
#ifdef OLDIO
  int nDomains     = -1;
#endif
  int reduceFactor =  1;
#ifndef PARTICLESRENDERER
  std::string fullScreenMode    = "";
  bool stereo     = false;
#endif


  if (rank == 0)  
  {
		AnyOption opt;

#define ADDUSAGE(line) {{std::stringstream oss; oss << line; opt.addUsage(oss.str());}}

		ADDUSAGE(" ");
		ADDUSAGE("Usage:");
		ADDUSAGE(" ");
		ADDUSAGE(" -h  --help             Prints this help ");
		ADDUSAGE(" -i  --infile #         Input snapshot filename ");
#ifdef OLDIO
    ADDUSAGE(" -n  --ndomains #       Number of domains ");
#endif
		ADDUSAGE(" -r  --reducefactor #   cut down bodies dataset by # factor [1])")
#ifndef PARTICLESRENDERER
		ADDUSAGE("     --fullscreen #     set fullscreen mode string");
		ADDUSAGE("     --stereo           enable stereo rendering");
#endif
		ADDUSAGE(" ");


		opt.setFlag  ( "help" ,        'h');
		opt.setOption( "infile",       'i');
#ifdef OLDIO
		opt.setOption( "ndomains",     'n');
#endif
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
#ifdef OLDIO
    if ((optarg = opt.getValue("ndomains")))     nDomains           = atoi(optarg);
#endif
    if ((optarg = opt.getValue("reducefactor"))) reduceFactor       = atoi(optarg);
#ifndef PARTICLESRENDERER
    if ((optarg = opt.getValue("fullscreen")))	 fullScreenMode     = std::string(optarg);
    if (opt.getFlag("stereo"))  stereo = true;
#endif

    if (fileName.empty() ||
#ifdef OLDIO
        nDomains == -1 || 
#endif
        reduceFactor < 1)
    {
      opt.printUsage();
      ::exit(0);
    }

#undef ADDUSAGE
  }


#ifdef OLDIO
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

  const int nPtcl = nFirstLocal+ nSecondLocal;
  RendererData rData(nPtcl);
  for (int i = 0; i < nFirstLocal; i++)
  {
    const int ip = i;
    rData.posx(ip) = data.firstPos[i].x;
    rData.posy(ip) = data.firstPos[i].y;
    rData.posz(ip) = data.firstPos[i].z;
    rData.ID  (ip) = data.firstID[i];
    rData.type(ip) = 0;
    rData.attribute(RendererData::MASS, ip) = data.firstPos[i].w;
    rData.attribute(RendererData::VEL,  ip) =
      std::sqrt(
          data.firstVel[i].x*data.firstVel[i].x +
          data.firstVel[i].y*data.firstVel[i].y +
          data.firstVel[i].z*data.firstVel[i].z);
  }
  for (int i = 0; i < nSecondLocal; i++)
  {
    const int ip = i + nFirstLocal;
    rData.posx(ip) = data.secondPos[i].x;
    rData.posy(ip) = data.secondPos[i].y;
    rData.posz(ip) = data.secondPos[i].z;
    rData.ID  (ip) = data.secondID[i];
    rData.type(ip) = 1;
    rData.attribute(RendererData::MASS, ip) = data.secondPos[i].w;
    rData.attribute(RendererData::VEL,  ip) =
      std::sqrt(
          data.secondVel[i].x*data.secondVel[i].x +
          data.secondVel[i].y*data.secondVel[i].y +
          data.secondVel[i].z*data.secondVel[i].z);
  }
#else
  BonsaiIO::Core out(rank, nranks, comm, BonsaiIO::READ, fileName);
  if (rank == 0)
    out.getHeader().printFields();
  typedef float float4[4];
  typedef float float3[3];
  typedef float float2[2];
  BonsaiIO::DataType<IDType> IDList("IDType");
  BonsaiIO::DataType<float4> pos("POS:real4");
  BonsaiIO::DataType<float3> vel("VEL:float[3]");
  BonsaiIO::DataType<float2> rhoh("RHOH:float[2]");

  assert(out.read(IDList, true, reduceFactor));
  assert(out.read(pos,    true, reduceFactor));
  assert(out.read(vel,    true, reduceFactor));
  bool renderDensity = true;
  if (!out.read(rhoh,  true, reduceFactor))
  {
    if (rank == 0)
    {
      fprintf(stderr , " -- RHOH data is found \n");
      fprintf(stderr , " -- rendering w/o density info \n");
    }
    renderDensity = false;
  }
  const int nPtcl = IDList.getNumElements();
  assert(IDList.getNumElements() == pos.getNumElements());
  assert(IDList.getNumElements() == vel.getNumElements());
  if (renderDensity)
    assert(IDList.getNumElements() == pos.getNumElements());
  RendererData rData(nPtcl);
  for (int i = 0; i < nPtcl; i++)
  {
    const int ip = i;
    rData.posx(ip) = pos[i][0];
    rData.posy(ip) = pos[i][1];
    rData.posz(ip) = pos[i][2];
    rData.ID  (ip) = IDList[i].getID();
    rData.type(ip) = IDList[i].getType();
    rData.attribute(RendererData::MASS, ip) = pos[i][3];
    rData.attribute(RendererData::VEL,  ip) =
      std::sqrt(
          vel[i][0]*vel[i][0] +
          vel[i][1]*vel[i][1] +
          vel[i][2]*vel[i][2]);
    if (renderDensity)
    {
    rData.attribute(RendererData::RHO, ip) = rhoh[i][0];
    rData.attribute(RendererData::H,  ip)  = rhoh[i][1];
    }
  }
#endif
  rData.computeMinMax();

  initAppRenderer(argc, argv, rData
#ifndef PARTICLESRENDERER
      ,fullScreenMode.c_str(), stereo
#endif
      );

  fprintf(stderr, " -- Done -- \n");
  while(1) {}
  return 0;
}



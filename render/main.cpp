#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <unistd.h>
#include <sstream>
#include <cmath>
#include "IDType.h"
#include "BonsaiIO.h"

#include "renderloop.h"
#include "anyoption.h"
#include "RendererData.h"

#include <IceT.h>
#include <IceTGL.h>
#include <IceTMPI.h>

template<typename T>
static T* readBonsaiReduced(
    const int rank, const int nranks, const MPI_Comm &comm,
    const std::string &fileName,
    const int reduceDM,
    const int reduceS)
{
  if (rank == 0)
    fprintf(stderr, " ----------- \n");
  BonsaiIO::Core in(rank, nranks, comm, BonsaiIO::READ, fileName);
  if (rank == 0)
    in.getHeader().printFields();
  typedef float float5[5];
  typedef float float3[3];

  BonsaiIO::DataType<float3> posStars("Stars:XYZ:float[3]");
  BonsaiIO::DataType<float5> attrStars("Stars:VxVyVz,DENS,H:float[5]");

  if (reduceS > 0)
  {
    if (rank  == 0)
      fprintf(stderr, " Reading star data \n");
    if(!in.read(posStars,  true, reduceS)) return NULL;
    assert(in.read(attrStars, true, reduceS));
    assert(posStars.getNumElements() == attrStars.getNumElements());
  }

  BonsaiIO::DataType<float3> posDM("DM:XYZ:float[3]");
  BonsaiIO::DataType<float5> attrDM("DM:VxVyVz,DENS,H:float[5]");

  if (reduceDM > 0)
  {
    if (rank  == 0)
      fprintf(stderr, " Reading DM data \n");
    if (!in.read(posDM,  true, reduceDM)) return NULL;
    assert(in.read(attrDM, true, reduceDM));
    assert(posDM.getNumElements() == attrDM.getNumElements());
  }



  const int nS  = posStars.getNumElements();
  const int nDM = posDM.getNumElements();
  long long int nSloc = nS, nSglb;
  long long int nDMloc = nDM, nDMglb;

  MPI_Allreduce(&nSloc, &nSglb, 1, MPI_LONG, MPI_SUM, comm);
  MPI_Allreduce(&nDMloc, &nDMglb, 1, MPI_LONG, MPI_SUM, comm);
  if (rank == 0)
  {
    fprintf(stderr, "nStars = %lld\n", nSglb);
    fprintf(stderr, "nDM    = %lld\n", nDMglb);
  }


  T *rDataPtr = new T(rank,nranks,comm);
  rDataPtr->resize(nS+nDM);
  auto &rData = *rDataPtr;
  for (int i = 0; i < nS; i++)
  {
    const int ip = i;
    rData.posx(ip) = posStars[i][0];
    rData.posy(ip) = posStars[i][1];
    rData.posz(ip) = posStars[i][2];
    rData.ID  (ip) = i;
    rData.type(ip) = 1;
//    rData.attribute(RendererData::MASS, ip) = 1.0/nS;
    rData.attribute(RendererData::VEL,  ip) =
      std::sqrt(
          attrStars[i][0]*attrStars[i][0] +
          attrStars[i][1]*attrStars[i][1] +
          attrStars[i][2]*attrStars[i][2]);
    rData.attribute(RendererData::RHO, ip) = attrStars[i][3];
    rData.attribute(RendererData::H,   ip) = attrStars[i][4];
  }

  for (int i = 0; i < nDM; i++)
  {
    const int ip = i + nS;
    rData.posx(ip) = posDM[i][0];
    rData.posy(ip) = posDM[i][1];
    rData.posz(ip) = posDM[i][2];
    rData.ID  (ip) = i+1000000000;
    rData.type(ip) = 0;
//    rData.attribute(RendererData::MASS, ip) = 1.0/nDM;
    rData.attribute(RendererData::VEL,  ip) =
      std::sqrt(
          attrDM[i][0]*attrDM[i][0] +
          attrDM[i][1]*attrDM[i][1] +
          attrDM[i][2]*attrDM[i][2]);
    rData.attribute(RendererData::RHO, ip) = attrDM[i][3];
    rData.attribute(RendererData::H,   ip) = attrDM[i][4];
  }
  return rDataPtr;
}

template<typename T>
static T* readBonsaiFull(
    const int rank, const int nranks, const MPI_Comm &comm,
    const std::string &fileName,
    const int reduceDM,
    const int reduceS)
{
  if (rank == 0)
    fprintf(stderr, " ----------- \n");
  BonsaiIO::Core out(rank, nranks, comm, BonsaiIO::READ, fileName);
  if (rank == 0)
    out.getHeader().printFields();
  typedef float float4[4];
  typedef float float3[3];
  typedef float float2[2];

  BonsaiIO::DataType<IDType> IDListS("Stars:IDType");
  BonsaiIO::DataType<float4> posS("Stars:POS:real4");
  BonsaiIO::DataType<float3> velS("Stars:VEL:float[3]");
  BonsaiIO::DataType<float2> rhohS("Stars:RHOH:float[2]");

  if (reduceS > 0)
  {
    if (rank  == 0)
      fprintf(stderr, " Reading star data \n");
    if (!out.read(IDListS, true, reduceS)) return NULL;
    assert(out.read(posS,    true, reduceS));
    assert(out.read(velS,    true, reduceS));
    bool renderDensity = true;
    if (!out.read(rhohS,  true, reduceS))
    {
      if (rank == 0)
      {
        fprintf(stderr , " -- Stars RHOH data is found \n");
        fprintf(stderr , " -- rendering stars w/o density info \n");
      }
      renderDensity = false;
    }
    assert(IDListS.getNumElements() == posS.getNumElements());
    assert(IDListS.getNumElements() == velS.getNumElements());
    if (renderDensity)
      assert(IDListS.getNumElements() == posS.getNumElements());
  }

  BonsaiIO::DataType<IDType> IDListDM("DM:IDType");
  BonsaiIO::DataType<float4> posDM("DM:POS:real4");
  BonsaiIO::DataType<float3> velDM("DM:VEL:float[3]");
  BonsaiIO::DataType<float2> rhohDM("DM:RHOH:float[2]");
  if (reduceDM > 0)
  {
    if (rank  == 0)
      fprintf(stderr, " Reading DM data \n");
    if(!out.read(IDListDM, true, reduceDM)) return NULL;
    assert(out.read(posDM,    true, reduceDM));
    assert(out.read(velDM,    true, reduceDM));
    bool renderDensity = true;
    if (!out.read(rhohDM,  true, reduceDM))
    {
      if (rank == 0)
      {
        fprintf(stderr , " -- DM RHOH data is found \n");
        fprintf(stderr , " -- rendering stars w/o density info \n");
      }
      renderDensity = false;
    }
    assert(IDListS.getNumElements() == posS.getNumElements());
    assert(IDListS.getNumElements() == velS.getNumElements());
    if (renderDensity)
      assert(IDListS.getNumElements() == posS.getNumElements());
  }


  const int nS  = IDListS.getNumElements();
  const int nDM = IDListDM.getNumElements();
  long long int nSloc = nS, nSglb;
  long long int nDMloc = nDM, nDMglb;

  MPI_Allreduce(&nSloc, &nSglb, 1, MPI_LONG, MPI_SUM, comm);
  MPI_Allreduce(&nDMloc, &nDMglb, 1, MPI_LONG, MPI_SUM, comm);
  if (rank == 0)
  {
    fprintf(stderr, "nStars = %lld\n", nSglb);
    fprintf(stderr, "nDM    = %lld\n", nDMglb);
  }


  T *rDataPtr = new T(rank,nranks,comm);
  rDataPtr->resize(nS+nDM);
  auto &rData = *rDataPtr;
  for (int i = 0; i < nS; i++)
  {
    const int ip = i;
    rData.posx(ip) = posS[i][0];
    rData.posy(ip) = posS[i][1];
    rData.posz(ip) = posS[i][2];
    rData.ID  (ip) = IDListS[i].getID();
    rData.type(ip) = IDListS[i].getType();
    assert(rData.type(ip) == 1); /* sanity check */
//    rData.attribute(RendererData::MASS, ip) = posS[i][3];
    rData.attribute(RendererData::VEL,  ip) =
      std::sqrt(
          velS[i][0]*velS[i][0] +
          velS[i][1]*velS[i][1] +
          velS[i][2]*velS[i][2]);
    if (rhohS.size() > 0)
    {
      rData.attribute(RendererData::RHO, ip) = rhohS[i][0];
      rData.attribute(RendererData::H,  ip)  = rhohS[i][1];
    }
    else
    {
      rData.attribute(RendererData::RHO, ip) = 0.0;
      rData.attribute(RendererData::H,   ip) = 0.0;
    }
  }
  for (int i = 0; i < nDM; i++)
  {
    const int ip = i + nS;
    rData.posx(ip) = posDM[i][0];
    rData.posy(ip) = posDM[i][1];
    rData.posz(ip) = posDM[i][2];
    rData.ID  (ip) = IDListDM[i].getID();
    rData.type(ip) = IDListDM[i].getType();
    assert(rData.type(ip) == 0); /* sanity check */
//    rData.attribute(RendererData::MASS, ip) = posDM[i][3];
    rData.attribute(RendererData::VEL,  ip) =
      std::sqrt(
          velDM[i][0]*velDM[i][0] +
          velDM[i][1]*velDM[i][1] +
          velDM[i][2]*velDM[i][2]);
    if (rhohDM.size() > 0)
    {
      rData.attribute(RendererData::RHO, ip) = rhohDM[i][0];
      rData.attribute(RendererData::H,   ip) = rhohDM[i][1];
    }
    else
    {
      rData.attribute(RendererData::RHO, ip) = 0.0;
      rData.attribute(RendererData::H,   ip) = 0.0;
    }
  }

  return rDataPtr;
}



int main(int argc, char * argv[])
{
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);
    
  int nranks, rank;
  MPI_Comm_size(comm, &nranks);
  MPI_Comm_rank(comm, &rank);


  std::string fileName;
  int reduceDM    =  10;
  int reduceS=  1;
#ifndef PARTICLESRENDERER
  std::string fullScreenMode    = "";
  bool stereo     = false;
#endif
  int nmaxsample = 200000;
  bool doDD = false;

  {
		AnyOption opt;

#define ADDUSAGE(line) {{std::stringstream oss; oss << line; opt.addUsage(oss.str());}}

		ADDUSAGE(" ");
		ADDUSAGE("Usage:");
		ADDUSAGE(" ");
		ADDUSAGE(" -h  --help             Prints this help ");
		ADDUSAGE(" -i  --infile #         Input snapshot filename ");
		ADDUSAGE("     --reduceDM    #    cut down DM dataset by # factor [10]. 0-disable DM");
		ADDUSAGE("     --reduceS     #    cut down stars dataset by # factor [1]. 0-disable S");
#ifndef PARTICLESRENDERER
		ADDUSAGE("     --fullscreen  #    set fullscreen mode string");
		ADDUSAGE("     --stereo           enable stereo rendering");
#endif
		ADDUSAGE(" -d  --doDD             enable domain decomposition  [disabled]");
    ADDUSAGE(" -s  --nmaxsample   #   set max number of samples for DD [" << nmaxsample << "]");


		opt.setFlag  ( "help" ,        'h');
		opt.setOption( "infile",       'i');
		opt.setOption( "reduceDM");
		opt.setOption( "reduceS");
    opt.setOption( "fullscreen");
    opt.setFlag("stereo");
    opt.setFlag("doDD", 'd');
    opt.setOption("nmaxsample", 's');

    opt.processCommandArgs( argc, argv );


    if( ! opt.hasOptions() ||  opt.getFlag( "help" ) || opt.getFlag( 'h' ) )
    {
      /* print usage if no options or requested help */
      opt.printUsage();
      ::exit(0);
    }

    char *optarg = NULL;
    if ((optarg = opt.getValue("infile")))       fileName           = std::string(optarg);
    if ((optarg = opt.getValue("reduceDM"))) reduceDM       = atoi(optarg);
    if ((optarg = opt.getValue("reduceS"))) reduceS       = atoi(optarg);
#ifndef PARTICLESRENDERER
    if ((optarg = opt.getValue("fullscreen")))	 fullScreenMode     = std::string(optarg);
    if (opt.getFlag("stereo"))  stereo = true;
#endif
    if ((optarg = opt.getValue("nmaxsample"))) nmaxsample = atoi(optarg);
    if (opt.getFlag("doDD"))  doDD = true;

    if (fileName.empty() ||
        reduceDM < 0 || reduceS < 0)
    {
      opt.printUsage();
      ::exit(0);
    }

#undef ADDUSAGE
  }


  using RendererDataT = RendererDataDistribute;
  RendererDataT *rDataPtr;
  if ((rDataPtr = readBonsaiFull<RendererDataT>(rank, nranks, comm, fileName, reduceDM, reduceS))) {}
  else if ((rDataPtr = readBonsaiReduced<RendererDataT>(rank, nranks, comm, fileName, reduceDM, reduceS))) {}
  else
  {
    if (rank == 0)
      fprintf(stderr, " I don't recognize the format ... please try again , or recompile to use with old tipsy if that is what you use ..\n");
    MPI_Finalize();
    exit(-1);
  }


  assert(rDataPtr != 0);
  rDataPtr->randomShuffle();
  rDataPtr->computeMinMax();

#if 0
  fprintf(stderr, "rank= %d: min= %g %g %g  max= %g %g %g \n",
      rank, 
      rDataPtr->xminLoc(), rDataPtr->yminLoc(), rDataPtr->zminLoc(),
      rDataPtr->xmaxLoc(), rDataPtr->ymaxLoc(), rDataPtr->zmaxLoc());
#endif
  fprintf(stderr, " rank= %d: n= %d\n", rank, rDataPtr->n());
  if (doDD)
  {
    MPI_Barrier(comm);
    const double t0 = MPI_Wtime();
    rDataPtr->setNMAXSAMPLE(nmaxsample);
    rDataPtr->distribute();
//    rDataPtr->distribute();
    MPI_Barrier(comm);
    const double t1 = MPI_Wtime();
    fprintf(stderr, " rank= %d: n= %d\n", rank, rDataPtr->n());
    if (rank == 0)
      fprintf(stderr, " DD= %g sec \n", t1-t0);
  }
#if 0
  fprintf(stderr, "rank= %d: min= %g %g %g  max= %g %g %g \n",
      rank, 
      rDataPtr->xminLoc(), rDataPtr->yminLoc(), rDataPtr->zminLoc(),
      rDataPtr->xmaxLoc(), rDataPtr->ymaxLoc(), rDataPtr->zmaxLoc());
#endif

  if (rDataPtr->attributeMin(RendererData::RHO) > 0.0)
  {
    rDataPtr->rescaleLinear(RendererData::RHO, 0, 60000.0);
    rDataPtr->scaleLog(RendererData::RHO);
  }
  rDataPtr->rescaleLinear(RendererData::VEL, 0, 3000.0);
//  rDataPtr->scaleLog(RendererData::VEL);

//  rDataPtr->clamp(RendererData::VEL, 0.25, 0.25);
//  rDataPtr->scaleExp(RendererData::VEL);
  

// #if 1
//   initAppRenderer(argc, argv, 
//       rank, nranks, comm,
//       *rDataPtr,
//       fullScreenMode.c_str(), stereo);
// #else
// 
//   sleep(1);
//   MPI_Barrier(comm);
//   if (rank == 0)
//     fprintf(stderr, " -- Done -- \n");
//   MPI_Finalize();
// #endif
  
  
    //Setup the IceT context and communicators
  IceTCommunicator icetComm    = icetCreateMPICommunicator(MPI_COMM_WORLD);
  IceTContext      icetContext = icetCreateContext(icetComm);
  icetDestroyMPICommunicator(icetComm); //Save since the comm is copied to the icetContext
  
  

  initAppRenderer(argc, argv, *rDataPtr
#ifndef PARTICLESRENDERER
      ,fullScreenMode.c_str(), stereo
#endif
      );
  
  icetGLInitialize();
  

  //Start the visualization
  initAppRenderer_start();
  
  
  
  
  
  while(1) {}
  return 0;
}



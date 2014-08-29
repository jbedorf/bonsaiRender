/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

/*
   This class renders particles using OpenGL and GLSL shaders
   */

#if 0
#define _SPLOTCHSPRITES
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cassert>
#include "renderer.h"
#include "SmokeShaders.h"
//#include <nvImage.h>
#include "depthSort.h"
#include "Cubemap.h"
#include <array>
#include <algorithm>

#if defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#endif

#define COLOR_ATTENUATION 0
#define USE_MRT 0
#define USE_HALF_ANGLE 0
#define MOTION_BLUR 0

//#define NOCOPY

//extern int renderDevID;
//extern int devID;

using namespace nv;
  
  template<typename T>
static inline double4 lMatVec(const T _m[16], const double4 pos)
{
  const T (*m)[4] = (T (*)[4])_m;
  return make_double4(
      m[0][0]*pos.x + m[1][0]*pos.y + m[2][0]*pos.z + m[3][0]*pos.w,
      m[0][1]*pos.x + m[1][1]*pos.y + m[2][1]*pos.z + m[3][1]*pos.w,
      m[0][2]*pos.x + m[1][2]*pos.y + m[2][2]*pos.z + m[3][2]*pos.w,
      m[0][3]*pos.x + m[1][3]*pos.y + m[2][3]*pos.z + m[3][3]*pos.w);
}

SmokeRendererParams::SmokeRendererParams() :
  mParticleRadius(0.1f),
  mParticleScaleLog(0.0f),
  mDisplayMode(SPLOTCH),
  mWindowW(800),
  mWindowH(600),
  mFov(40.0f),
  m_downSample(1),
  m_blurDownSample(2),
  m_numSlices(64),
  m_numDisplayedSlices(m_numSlices),
  m_sliceNo(0),
  m_shadowAlpha(0.1f),
  m_volumeAlpha(0.2f),
  m_dustAlpha(1.0f),
  m_volumeColor(0.5f, 0.0f, 0.0f),
  m_doBlur(true),
  m_blurRadius(1.0f),
  m_blurPasses(2),
  m_displayLightBuffer(false),
  m_lightPos(5.0f, 5.0f, -5.0f),
  m_lightTarget(0.0f, 0.5f, 0.0f),
  m_lightColor(0.0f, 0.0f, 0.0f),
  m_colorOpacity(0.1f, 0.2f, 0.3f),
  /***************/

  m_starScaleLog(0.0f),
  m_starAlpha(1.0f),
  m_dmScaleLog(-0.4f),
  m_dmAlpha(0.1f),
  m_spriteAlpha(0.02f),
  m_transmission(0.001f),
  m_imageBrightnessPre(0.08f),
  m_gammaPre(0.4f),
  m_imageBrightnessPost(1.0f),
  m_gammaPost(1.0f),

  /***************/

  m_overBright(1.0f),
  m_overBrightThreshold(0.005f),
  m_imageBrightness(1.0f),
  m_invertedView(false),
  m_minDepth(0.0f),
  m_maxDepth(1.0f),
  m_enableAA(false),
  m_starBlurRadius(40.0f),
  m_starThreshold(1.0f),
  m_starPower(1.0f),
  m_starIntensity(0.5f),
  m_glowRadius(10.0f),
  m_glowIntensity(0.5f),
  m_ageScale(10.0f),
  m_enableVolume(false),
//  m_enableFilters(true),
  m_enableFilters(true),
  m_noiseFreq(0.05f),
  m_noiseAmp(1.0f),
  m_indirectAmount(0.5f),
  m_volumeIndirect(0.5f),
  m_volumeStart(0.5f),
  m_volumeWidth(0.1f),
  m_gamma(1.0f / 2.2f),
  m_fog(0.001f),
  //    m_cubemapTex(0),
  m_flareThreshold(0.5f),
  m_flareIntensity(0.0f),
  m_sourceIntensity(0.5f),
  m_flareRadius(50.0f),
  m_skyboxBrightness(0.5f),
  m_cullDarkMatter(true),
  m_doClipping(false),
  m_domainView(false)
{}


SmokeRenderer::SmokeRenderer(int numParticles, int maxParticles, const int _rank, const int _nrank, const MPI_Comm &_comm) :
  SmokeRendererParams(),
  rank(_rank), nrank(_nrank), comm(_comm),
  mMaxParticles(maxParticles),
  mNumParticles(numParticles),
  mPosVbo(0),
  m_pbo(0),
  mVelVbo(0),
  mColorVbo(0),
  mSizeVbo(0),
  mSizeVao(0),
  mIndexBuffer(0),
  m_lightBufferSize(512),
  m_srcLightTexture(0),
  m_lightDepthTexture(0),
  m_fbo(0),
  m_depthTex(0),
  m_rampTex(0)

#if 0
  /***************/

  m_starScaleLog(0.0f),
  m_starAlpha(1.0f),
  m_dmScaleLog(-0.4f),
  m_dmAlpha(0.1f),
  m_spriteAlpha(0.02f),
  m_transmission(0.001f),
  m_imageBrightnessPre(0.08f),
  m_gammaPre(0.4f),
  m_imageBrightnessPost(1.0f),
  m_gammaPost(1.0f),

  /***************/

  m_overBright(1.0f),
  m_overBrightThreshold(0.005f),
  m_imageBrightness(1.0f),
  m_invertedView(false),
  m_minDepth(0.0f),
  m_maxDepth(1.0f),
  m_enableAA(false),
  m_starBlurRadius(40.0f),
  m_starThreshold(1.0f),
  m_starPower(1.0f),
  m_starIntensity(0.5f),
  m_glowRadius(10.0f),
  m_glowIntensity(0.5f),
  m_ageScale(10.0f),
  m_enableVolume(false),
//  m_enableFilters(true),
  m_enableFilters(true),
  m_noiseFreq(0.05f),
  m_noiseAmp(1.0f),
  m_indirectAmount(0.5f),
  m_volumeIndirect(0.5f),
  m_volumeStart(0.5f),
  m_volumeWidth(0.1f),
  m_gamma(1.0f / 2.2f),
  m_fog(0.001f),
  //    m_cubemapTex(0),
  m_flareThreshold(0.5f),
  m_flareIntensity(0.0f),
  m_sourceIntensity(0.5f),
  m_flareRadius(50.0f),
  m_skyboxBrightness(0.5f),
  m_cullDarkMatter(true)
#endif
{
  assert(rank < nrank);
  // load shader programs
  m_simpleProg = new GLSLProgram(simpleVS, simplePS);
#if MOTION_BLUR
  m_particleProg = new GLSLProgram(mblurVS, mblurGS, particlePS);
  m_particleAAProg = new GLSLProgram(mblurVS, mblurGS, particleAAPS);
  m_particleShadowProg = new GLSLProgram(mblurVS, mblurGS, particleShadowPS);
#else
  m_particleProg = new GLSLProgram(particleVS, particlePS);
  m_particleAAProg = new GLSLProgram(particleVS, particleAAPS);
  m_particleShadowProg = new GLSLProgram(particleVS, particleShadowPS);
#endif

  //m_blurProg = new GLSLProgram(passThruVS, blur2PS);
  m_blurProg = new GLSLProgram(passThruVS, blur3x3PS);
  m_displayTexProg = new GLSLProgram(passThruVS, texture2DPS);
  m_compositeProg = new GLSLProgram(passThruVS, compositePS);

  m_starFilterProg = new GLSLProgram(passThruVS, starFilterPS);
  m_volumeProg = new GLSLProgram(volumeVS, volumePS);

  //m_downSampleProg = new GLSLProgram(passThruVS, downSample4PS);
  m_downSampleProg = new GLSLProgram(passThruVS, downSample2PS);
  m_gaussianBlurProg = new GLSLProgram(passThruVS, gaussianBlurPS);
  m_thresholdProg = new GLSLProgram(passThruVS, thresholdPS);

  m_skyboxProg = new GLSLProgram(skyboxVS, skyboxPS);

#ifdef _SPLOTCHSPRITES  /* set to 1 if you want to use point sprites */
  m_splotchProg = new GLSLProgram(splotchVS, splotchPS);
#else
  m_splotchProg = new GLSLProgram(splotchVS, splotchGS, splotchPS,
//      GL_POINTS, GL_POINTS
      GL_POINTS, GL_TRIANGLE_STRIP
      );
#endif
  m_splotch2texProg = new GLSLProgram(passThruVS, splotch2texPS);

  glClampColorARB(GL_CLAMP_VERTEX_COLOR_ARB, GL_FALSE);
  glClampColorARB(GL_CLAMP_FRAGMENT_COLOR_ARB, GL_FALSE);

  m_fbo = new FramebufferObject();
  // create buffer for light shadows
  //createLightBuffer();

  // textures
  //    loadSmokeTextures(32,0,"perlinNoiseTextures2/noise2_");
  //    loadSmokeTextures(8, 0, "noiseTextures3/noise");
  //    m_rainbowTex = loadTexture("data/rainbow.png");
  m_rainbowTex = createRainbowTexture();

  glGenTextures(1, &mPosBufferTexture);
  m_noiseTex = createNoiseTexture(64, 64, 64);

#if 0
  m_cubemapTex = loadCubemapCross("../images/Carina_cross.ppm");
  if (!m_cubemapTex) {
    //m_cubemapTex = loadCubemap("../images/deepfield%d.ppm");
    m_cubemapTex = loadCubemap("../images/deepfield%d_1k.ppm");
  }
#endif

  m_spriteTex = createSpriteTexture(256);
  m_sphTex    = createSphTexture(256);

  initParams();

  //initCUDA();

  //	cudaGLSetGLDevice(renderDevID);

  //	mParticlePos.alloc(mMaxParticles, true, false, false);
  //	mParticleDepths.alloc(mMaxParticles, false, false, false);
  //	mParticleIndices.alloc(mMaxParticles, true, false, true);
  //	for(uint i=0; i<mMaxParticles; i++) {
  //		mParticleIndices.getHostPtr()[i] = i;
  //	}
  //	mParticleIndices.copy(GpuArray<uint>::HOST_TO_DEVICE);

  //	cudaStreamCreate(&m_copyStreamPos);
  //  cudaStreamCreate(&m_copyStreamColor);

  //  cudaStreamCreate(&m_copyStreamSortPos);
  //  cudaStreamCreate(&m_copyStreamSortDepth);
  //  cudaStreamCreate(&m_copyStreamSortIndices);

  //  cudaDeviceEnablePeerAccess(devID, 0);
  //	cudaSetDevice(devID);

  //	cudaDeviceEnablePeerAccess( renderDevID, 0 );

  //Allocate additional arrays
  mParticlePos  = (float4*)malloc(mMaxParticles*sizeof(float4));
  mParticleColors  = (float4*)malloc(mMaxParticles*sizeof(float4));
  mParticleDepths  = (float*)malloc(mMaxParticles*sizeof(float));
  mParticleIndices = (uint*)malloc(mMaxParticles*sizeof(uint));
  //  cudaMalloc( &mParticleDepths_devID, mMaxParticles*sizeof(float));
  //  cudaMalloc( &mParticleIndices_devID, mMaxParticles*sizeof(uint));

  if (isMaster())
  {
    printf("Vendor: %s\n", glGetString(GL_VENDOR));
    printf("Renderer: %s\n", glGetString(GL_RENDERER));
  }

  glutReportErrors();
}

SmokeRenderer::~SmokeRenderer()
{
  delete m_particleProg;
  delete m_particleAAProg;
  delete m_particleShadowProg;
  delete m_blurProg;
  delete m_displayTexProg;
  delete m_simpleProg;
  delete m_starFilterProg;
  delete m_compositeProg;
  delete m_volumeProg;
  delete m_skyboxProg;

  delete m_splotchProg;
  delete m_splotch2texProg;

  delete m_fbo;
  glDeleteTextures(2, m_lightTexture);
  glDeleteTextures(1, &m_lightDepthTexture);
  glDeleteTextures(1, &mPosBufferTexture);

  glDeleteTextures(4, m_imageTex);
  glDeleteTextures(1, &m_depthTex);
  glDeleteTextures(3, m_downSampledTex);

  glDeleteTextures(1, &m_noiseTex);
  //	glDeleteTextures(1, &m_cubemapTex);

  //	cudaSetDevice(renderDevID);

  //	mParticlePos.free();
  //	mParticleDepths.free();
  //	mParticleIndices.free();
  free(mParticleIndices);
  free(mParticleDepths);
  free(mParticlePos);
  free(mParticleColors);
}

GLuint SmokeRenderer::createRainbowTexture()
{
  vec4f colors[] = {
    vec4f(1.0f, 0.0f, 0.0f, 1.0f),	// r
    vec4f(1.0f, 1.0f, 0.0f, 1.0f),	// y
    vec4f(0.0f, 1.0f, 0.0f, 1.0f),	// g
    vec4f(0.0f, 0.0f, 1.0f, 1.0f),	// b
    vec4f(1.0f, 0.0f, 1.0f, 1.0f),	// m
  };

  GLuint tex;
  glGenTextures(1, &tex);

  GLenum target = GL_TEXTURE_2D;
  glBindTexture( target, tex);

  glTexParameteri( target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri( target, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri( target, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri( target, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri( target, GL_GENERATE_MIPMAP_SGIS, GL_TRUE);

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glTexImage2D(target, 0, GL_RGBA, sizeof(colors) / sizeof(vec4f), 1, 0, GL_RGBA, GL_FLOAT, &colors[0]);
  return tex;
}

#if 0
GLuint SmokeRenderer::loadTexture(char *filename)
{
  nv::Image image;
  if (!image.loadImageFromFile(filename))  {
    printf( "Failed to load image '%s'\n", filename);
    return 0;
  }

  GLuint tex;
  glGenTextures(1, &tex);

  GLenum target = GL_TEXTURE_2D;
  glBindTexture( target, tex);

  glTexParameteri( target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri( target, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri( target, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri( target, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri( target, GL_GENERATE_MIPMAP_SGIS, GL_TRUE);

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glTexImage2D(target, 0, image.getInternalFormat(), image.getWidth(), image.getHeight(), 0, image.getFormat(), image.getType(), image.getLevel(0));
  return tex;
}

void SmokeRenderer::loadSmokeTextures(int nImages, int offset, char* sTexturePrefix)
{
  nv::Image* images  = new nv::Image[nImages];
  char textureName[260];
  std::string resolved_path;
  for(int i=0; i<nImages; i++) 
  {
    sprintf_s(textureName, "%s%.3d.png", sTexturePrefix, i+offset);
    //if ( pathHelper.getFilePath( textureName, resolved_path)) 
    {
      //if (! images[i].loadImageFromFile( resolved_path.c_str())) 
      if (! images[i].loadImageFromFile(textureName)) 
      {
        printf( "Failed to load smoke texture\n");
        exit(-1);
      }
      printf("Loaded '%s'\n", textureName);
    }
  }

  // load images as 2d texture array
  glGenTextures(1, &m_textureArrayID);
  glBindTexture( GL_TEXTURE_2D_ARRAY_EXT, m_textureArrayID);

  glTexParameteri( GL_TEXTURE_2D_ARRAY_EXT, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri( GL_TEXTURE_2D_ARRAY_EXT, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  glTexParameteri( GL_TEXTURE_2D_ARRAY_EXT, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri( GL_TEXTURE_2D_ARRAY_EXT, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri( GL_TEXTURE_2D_ARRAY_EXT, GL_GENERATE_MIPMAP_SGIS, GL_TRUE);

  // 2D Texture arrays a loaded just like 3D textures
  glTexImage3D(GL_TEXTURE_2D_ARRAY_EXT, 0, GL_LUMINANCE8, images[0].getWidth(), images[0].getHeight(), nImages, 0, images[0].getFormat(), images[0].getType(), NULL);
  for (int i = 0; i < nImages; i++) 
    glTexSubImage3D( GL_TEXTURE_2D_ARRAY_EXT, 0, 0, 0, i, images[i].getWidth(), images[i].getHeight(), 1, images[i].getFormat(), images[i].getType(), images[i].getLevel(0));

}
#endif

void SmokeRenderer::setNumberOfParticles(uint n_particles)
{
  if(n_particles > this->mMaxParticles)
  {
    //Uhohhh too many particles
    fprintf(stderr, "Sorry increase the number of maxParticles \n");
    this->mNumParticles = this->mMaxParticles;
  }
  else
  {
    this->mNumParticles = n_particles;
  }
}

void SmokeRenderer::setPositions(float *pos)
{
#if 0
  for (uint i = 0; i < mNumParticles; i++)
    mParticlePos[i] = ((float4*)pos)[i];
#endif
#if 0
  //memcpy(mParticlePos.getHostPtr(), pos, mNumParticles*4*sizeof(float));
  //ParticlePos.copy(GpuArray<float4>::HOST_TO_DEVICE);
#else
  // XXX - why is this so much faster?
  //    int posVbo = mParticlePos.getVbo();
  //    glBindBuffer(GL_ARRAY_BUFFER_ARB, posVbo);
  if (!m_pbo)
  {
    glGenBuffers(1, (GLuint*)&m_pbo);
    glBindBuffer(GL_ARRAY_BUFFER_ARB, m_pbo);
    glBufferData(GL_ARRAY_BUFFER_ARB, mNumParticles * 4 * sizeof(float), pos, GL_DYNAMIC_DRAW);
  }
  assert(m_pbo);
  glBindBuffer(GL_ARRAY_BUFFER_ARB, m_pbo);
  glBufferSubData(GL_ARRAY_BUFFER_ARB, 0, mNumParticles * 4 * sizeof(float), pos);
  glBindBuffer( GL_ARRAY_BUFFER_ARB, 0);
#endif
}


#if 0
void SmokeRenderer::setPositionsDevice(float *posD)
{
  //	cudaSetDevice(renderDevID);

#ifndef NOCOPY
  mParticlePos.map();
  //  cudaMemcpy(mParticlePos.getDevicePtr(), posD, mNumParticles*4*sizeof(float), cudaMemcpyDeviceToDevice);
  //    cudaMemcpyPeerAsync(mParticlePos.getDevicePtr(), renderDevID, posD, devID, mNumParticles*4*sizeof(float), m_copyStreamPos);
  mParticlePos.unmap();
#endif

  //	cudaSetDevice(devID);
}
#endif

void SmokeRenderer::setSizes(float *sizes)
{
  if (!mSizeVbo)
  {
    // allocate
    glGenBuffers(1, &mSizeVbo);
    glBindBuffer(GL_ARRAY_BUFFER_ARB, mSizeVbo);
    glBufferData(GL_ARRAY_BUFFER_ARB, mMaxParticles * sizeof(float), sizes, GL_DYNAMIC_DRAW);                
  }

  glBindBuffer(GL_ARRAY_BUFFER_ARB, mSizeVbo);
  glBufferSubData(GL_ARRAY_BUFFER_ARB, 0, mNumParticles * sizeof(float), sizes);
  glBindBuffer( GL_ARRAY_BUFFER_ARB, 0);
}

void SmokeRenderer::setColors(float *color)
{
#if 0
  for (uint i = 0; i < mNumParticles; i++)
    mParticleColors[i] = ((float4*)color)[i];
#endif
  if (!mColorVbo)
  {
    // allocate
    glGenBuffers(1, &mColorVbo);
    glBindBuffer(GL_ARRAY_BUFFER_ARB, mColorVbo);
    // 		glBufferData(GL_ARRAY_BUFFER_ARB, mNumParticles * 4 * sizeof(float), color, GL_DYNAMIC_DRAW);
    //Jeroen, I allocate the maximum number of particles
    glBufferData(GL_ARRAY_BUFFER_ARB, mMaxParticles * 4 * sizeof(float), color, GL_DYNAMIC_DRAW);                
  }

  glBindBuffer(GL_ARRAY_BUFFER_ARB, mColorVbo);
  //glBufferData(GL_ARRAY_BUFFER_ARB, mNumParticles * 4 * sizeof(float), color, GL_DYNAMIC_DRAW);
  glBufferSubData(GL_ARRAY_BUFFER_ARB, 0, mNumParticles * 4 * sizeof(float), color);
  glBindBuffer( GL_ARRAY_BUFFER_ARB, 0);
}

#if 0
void SmokeRenderer::setColorsDevice(float *colorD)
{
  //	cudaSetDevice(renderDevID);

  if (!mColorVbo)
  {
    // allocate
    glGenBuffers(1, &mColorVbo);
    glBindBuffer(GL_ARRAY_BUFFER_ARB, mColorVbo);
    // 		glBufferData(GL_ARRAY_BUFFER_ARB, mNumParticles * 4 * sizeof(float), color, GL_DYNAMIC_DRAW);
    //Jeroen, I allocate the maximum number of particles
    glBufferData(GL_ARRAY_BUFFER_ARB, mMaxParticles * 4 * sizeof(float), NULL, GL_DYNAMIC_DRAW);    
    cutilSafeCall(cudaGLRegisterBufferObject(mColorVbo));
    cutilSafeCall(cudaGLSetBufferObjectMapFlags(mColorVbo, cudaGLMapFlagsWriteDiscard));    // CUDA writes, GL consumes
  }

#ifndef NOCOPY
  void *ptr;
  cutilSafeCall(cudaGLMapBufferObject((void **) &ptr, mColorVbo));
  //cudaMemcpy( ptr, colorD, mNumParticles * 4 * sizeof(float), cudaMemcpyDeviceToDevice );
  //  cudaMemcpyPeerAsync( ptr, renderDevID, colorD, devID, mNumParticles * 4 * sizeof(float), m_copyStreamColor );
  // cutilSafeCall(cudaGLUnmapBufferObject(mColorVbo));
#endif

  //	cudaSetDevice(devID);
}
#endif

void SmokeRenderer::depthSort(float4 *pos)
{
  calcVectors();
  float4 modelViewZ = make_float4(m_modelView._array[2], m_modelView._array[6], m_modelView._array[10], m_modelView._array[14]);
  depthSortCUDA(pos, mParticleDepths, (int *) mParticleIndices, modelViewZ, mNumParticles);

  if (!mIndexBuffer)
  {
    glGenBuffersARB(1, (GLuint*)&mIndexBuffer);
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, mIndexBuffer);
    glBufferData(GL_ARRAY_BUFFER_ARB, mNumParticles * sizeof(uint), mParticleIndices, GL_DYNAMIC_DRAW);
  }
  glBindBuffer(GL_ARRAY_BUFFER_ARB, mIndexBuffer);
  glBufferSubData(GL_ARRAY_BUFFER_ARB, 0, mNumParticles * sizeof(uint), mParticleIndices);
  glBindBuffer( GL_ARRAY_BUFFER_ARB, 0);
}


void SmokeRenderer::depthSortCopy()
{
  //  depthSort(mParticlePos);
#if 0
  cudaSetDevice(renderDevID);

#ifndef NOCOPY
  mParticleIndices.map();

  //  cudaMemcpyPeerAsync(mParticleDepths.getDevicePtr(), renderDevID, mParticleDepths_devID, devID, mNumParticles*sizeof(float), m_copyStreamSortDepth);
  // cudaMemcpyPeerAsync(mParticleIndices.getDevicePtr(), renderDevID, mParticleIndices_devID, devID, mNumParticles*sizeof(uint), m_copyStreamSortIndices);

  mParticleIndices.unmap();
#endif

  //    mParticleIndices.map();
  //    mParticlePos.map();
  //
  //    float4 modelViewZ = make_float4(m_modelView._array[2], m_modelView._array[6], m_modelView._array[10], m_modelView._array[14]);
  //    depthSortCUDA(mParticlePos.getDevicePtr(), mParticleDepths.getDevicePtr(), (int *) mParticleIndices.getDevicePtr(), modelViewZ, mNumParticles);
  //
  //    mParticlePos.unmap();
  //    mParticleIndices.unmap();

  //	cudaSetDevice(devID);
#endif
}

// draw points from vertex buffer objects
void SmokeRenderer::drawPoints(int start, int count, bool sorted)
{
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_pbo);
  //    glBindBufferARB(GL_ARRAY_BUFFER_ARB, mParticlePos.getVbo());
  glVertexPointer(4, GL_FLOAT, 0, 0);
  glEnableClientState(GL_VERTEX_ARRAY);                

  if (mColorVbo) {
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, mColorVbo);
    glColorPointer(4, GL_FLOAT, 0, 0);
    glEnableClientState(GL_COLOR_ARRAY);
  }

  if (mVelVbo) {
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, mVelVbo);
    glClientActiveTexture(GL_TEXTURE0);
    glTexCoordPointer(4, GL_FLOAT, 0, 0);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  }

#if 0
  if (mSizeVbo)
  {
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, mSizeVbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 1, GL_FLOAT, GL_FALSE, 0, 0);
  }
#endif

  if (sorted) {
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, mIndexBuffer);
    glDrawElements(GL_POINTS, count, GL_UNSIGNED_INT, (void*) (start*sizeof(unsigned int)) );
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);
  } else {
    glDrawArrays(GL_POINTS, start, count);
  }

  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);

  glClientActiveTexture(GL_TEXTURE0);
  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
}

// draw points using given shader program
void SmokeRenderer::drawPointSprites(GLSLProgram *prog, int start, int count, bool shadowed, bool sorted)
{
  glEnable(GL_DEPTH_TEST);
  glDepthMask(GL_FALSE);  // don't write depth
  glEnable(GL_BLEND);

  GLuint vertexLoc = -1;
  if (!mSizeVao && mSizeVbo)
  {
    glGenVertexArrays(1, &mSizeVao);
    glBindVertexArray(mSizeVao);
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, mSizeVbo);
    vertexLoc = prog->getAttribLoc("pointRadiusAttr");
    glEnableVertexAttribArray(vertexLoc);
    glVertexAttribPointer(vertexLoc , 1, GL_FLOAT, 0, 0, 0);
  }

  prog->enable();
  glBindVertexArray(mSizeVao);
  prog->setUniform1f("pointRadius", mParticleRadius);
  prog->setUniform1f("ageScale", m_ageScale);
  prog->setUniform1f("dustAlpha", m_dustAlpha);
  prog->setUniform1f("overBright", m_overBright);
  prog->setUniform1f("overBrightThreshold", m_overBrightThreshold);
  prog->setUniform1f("fogDist", m_fog);
  prog->setUniform1f("cullDarkMatter", (float) m_cullDarkMatter);

  //prog->bindTexture("rampTex", m_rampTex, GL_TEXTURE_2D, 0);
  //prog->bindTexture("rampTex", m_rainbowTex, GL_TEXTURE_2D, 0);
  //prog->bindTexture("spriteTex",  m_textureArrayID, GL_TEXTURE_2D_ARRAY_EXT, 1);
  prog->bindTexture("spriteTex",  m_spriteTex, GL_TEXTURE_2D, 1);
  if (shadowed) {
    prog->bindTexture("shadowTex", m_lightTexture[m_srcLightTexture], GL_TEXTURE_2D, 2);
#if USE_MRT
    //prog->setUniform2f("shadowTexScale", m_lightBufferSize / (float) m_imageW, m_lightBufferSize / (float) m_imageH);
#else
    //prog->setUniform2f("shadowTexScale", 1.0f, 1.0f);
#endif
    prog->setUniform1f("indirectAmount", m_indirectAmount);
    prog->setUniform1f("alphaScale", m_spriteAlpha);

  } else {
    prog->setUniform1f("alphaScale", m_shadowAlpha);
    prog->setUniform1f("transmission", m_transmission);
  }

#if MOTION_BLUR==0
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  //prog->setUniform1f("pointScale", mWindowH / mInvFocalLen);
  prog->setUniform1f("pointScale", viewport[3] / mInvFocalLen);

  //glClientActiveTexture(GL_TEXTURE0);
  glActiveTexture(GL_TEXTURE0);
  glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);
  glEnable(GL_POINT_SPRITE_ARB);
#endif

  // draw points
  drawPoints(start, count, sorted);

  prog->disable();

  glDisable(GL_POINT_SPRITE_ARB);
  glDisable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);
  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);
  glBindVertexArray(0);

}

// calculate vectors for half-angle slice rendering
void SmokeRenderer::calcVectors()
{
  // get model view matrix
  glGetFloatv(GL_MODELVIEW_MATRIX, (float *) m_modelView.get_value());

  // calculate eye space light vector
  m_lightVector = normalize(m_lightPos);
  m_lightPosEye = m_modelView * vec4f(m_lightPos, 1.0);

  m_viewVector = -vec3f(m_modelView.get_row(2));
#if USE_HALF_ANGLE
  // calculate half-angle vector between view and light
  if (dot(m_viewVector, m_lightVector) > 0) {
    m_halfVector = normalize(m_viewVector + m_lightVector);
    m_invertedView = false;
  } else {
    m_halfVector = normalize(-m_viewVector + m_lightVector);
    m_invertedView = true;
  }

  // calculate light view matrix
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  gluLookAt(m_lightPos[0], m_lightPos[1], m_lightPos[2], 
      m_lightTarget[0], m_lightTarget[1], m_lightTarget[2],
      0.0, 1.0, 0.0);

  // calculate light projection matrix
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluPerspective(45.0, 1.0, 1.0, 200.0);

  glGetFloatv(GL_MODELVIEW_MATRIX, (float *) m_lightView.get_value());
  glGetFloatv(GL_PROJECTION_MATRIX, (float *) m_lightProj.get_value());

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
#else
  // camera-aligned slices
  m_halfVector = m_viewVector;
  m_lightView = m_modelView;
  glGetFloatv(GL_PROJECTION_MATRIX, (float *) m_lightProj.get_value());
  m_invertedView = false;
#endif

  // construct shadow matrix
  matrix4f scale;
  scale.set_scale(vec3f(0.5, 0.5, 0.5));
  matrix4f translate;
  translate.set_translate(vec3f(0.5, 0.5, 0.5));

  m_shadowMatrix = translate * scale * m_lightProj * m_lightView * inverse(m_modelView);

  // calc object space eye position
  m_eyePos = inverse(m_modelView) * vec4f(0.0, 0.0, 0.0, 1.0);

  // calc half vector in eye space
  m_halfVectorEye = m_modelView * vec4f(m_halfVector, 0.0);
}

// draw quad for volume rendering
void SmokeRenderer::drawVolumeSlice(int i, bool shadowed)
{
  glEnable(GL_BLEND);
  glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
  //glBlendFunc(GL_ONE, GL_ONE);
  glDisable(GL_DEPTH_TEST);

  //glMatrixMode(GL_MODELVIEW);
  //glPushMatrix();
  //glLoadIdentity();

  //glColor4f(1.0, 1.0, 1.0, m_volumeAlpha);
  glColor4f(m_volumeColor[0], m_volumeColor[1], m_volumeColor[2], m_volumeAlpha);

  m_volumeProg->enable();
  m_volumeProg->bindTexture("noiseTex", m_noiseTex, GL_TEXTURE_3D, 0);
  if (shadowed) {
    m_volumeProg->bindTexture("shadowTex", m_lightTexture[m_srcLightTexture], GL_TEXTURE_2D, 1);
  }
  m_volumeProg->setUniform1f("noiseFreq", m_noiseFreq);
  m_volumeProg->setUniform1f("noiseAmp", m_noiseAmp);
  m_volumeProg->setUniform1f("indirectLighting", shadowed ? m_volumeIndirect: 0.0f);
  m_volumeProg->setUniform1f("volumeStart", m_volumeStart);
  m_volumeProg->setUniform1f("volumeWidth", m_volumeWidth);

  float t = i / (float) m_numSlices;
  float z = m_minDepth + (m_maxDepth - m_minDepth) * t;
  drawQuad(0.5f, z);

  m_volumeProg->disable();

  //glPopMatrix();
  glDisable(GL_BLEND);
}

// draw slice of particles from camera view
void SmokeRenderer::drawSlice(int i)
{
#if USE_MRT
  GLenum buffers[] = { GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT };
  glDrawBuffers(2,buffers);

  glColorMaskIndexedEXT(0, true, true, true, true);
  glColorMaskIndexedEXT(1, false, false, false, false);
#else
  m_fbo->AttachTexture(GL_TEXTURE_2D, m_imageTex[0], GL_COLOR_ATTACHMENT0_EXT);
  //m_fbo->AttachTexture(GL_TEXTURE_2D, m_depthTex, GL_DEPTH_ATTACHMENT_EXT);
  m_fbo->AttachTexture(GL_TEXTURE_2D, 0, GL_DEPTH_ATTACHMENT_EXT);
#endif
  glViewport(0, 0, m_imageW, m_imageH);

  if (m_enableVolume) {
    drawVolumeSlice(i, true);
  }

  glColor4f(1.0, 1.0, 1.0, m_spriteAlpha);
  if (m_invertedView) {
    // front-to-back
    glBlendFunc(GL_ONE_MINUS_DST_ALPHA, GL_ONE);
  } else {
    // back-to-front
    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    //glBlendFunc(GL_ONE, GL_ONE);
  }
  drawPointSprites(m_particleShadowProg, i*m_batchSize, m_batchSize, true, true);
}

// draw slice of particles from light's point of view
void SmokeRenderer::drawSliceLightView(int i)
{
#if USE_MRT
  GLenum buffers[] = { GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT };
  glDrawBuffers(2, buffers);

  glColorMaskIndexedEXT(0, false, false, false, false);
  glColorMaskIndexedEXT(1, true, true, true, true);
#else
  m_fbo->AttachTexture(GL_TEXTURE_2D, m_lightTexture[m_srcLightTexture], GL_COLOR_ATTACHMENT0_EXT);
  //m_fbo->AttachTexture(GL_TEXTURE_2D, m_lightDepthTexture, GL_DEPTH_ATTACHMENT_EXT);
  m_fbo->AttachTexture(GL_TEXTURE_2D, 0, GL_DEPTH_ATTACHMENT_EXT);
#endif

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadMatrixf((GLfloat *) m_lightView.get_value());

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadMatrixf((GLfloat *) m_lightProj.get_value());

  glViewport(0, 0, m_lightBufferSize, m_lightBufferSize);

  if (m_enableVolume) {
    drawVolumeSlice(i, false);
  }

#if COLOR_ATTENUATION
  glColor4f(m_colorOpacity[0], m_colorOpacity[1], m_colorOpacity[2], m_shadowAlpha);
  glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_COLOR);
#else
  glColor4f(1.0, 1.0, 1.0, m_shadowAlpha);
  //    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
  //    glBlendFunc(GL_ONE, GL_ONE);	// additive
#endif

  drawPointSprites(m_particleProg, i*m_batchSize, m_batchSize, false, true);

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

// draw slice of particles from light's point of view
// version with anti-aliasing
void SmokeRenderer::drawSliceLightViewAA(int i)
{
#if USE_MRT
  GLenum buffers[] = { GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT };
  glDrawBuffers(2, buffers);

  glColorMaskIndexedEXT(0, false, false, false, false);
  glColorMaskIndexedEXT(1, true, true, true, true);
#else
  m_fbo->AttachTexture(GL_TEXTURE_2D, m_lightTexture[m_srcLightTexture], GL_COLOR_ATTACHMENT0_EXT);
  //m_fbo->AttachTexture(GL_TEXTURE_2D, m_lightDepthTexture, GL_DEPTH_ATTACHMENT_EXT);
  m_fbo->AttachTexture(GL_TEXTURE_2D, 0, GL_DEPTH_ATTACHMENT_EXT);
#endif

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadMatrixf((GLfloat *) m_lightView.get_value());

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadMatrixf((GLfloat *) m_lightProj.get_value());

  glViewport(0, 0, m_lightBufferSize, m_lightBufferSize);

#if 0
  glColor4f(m_colorOpacity[0] * m_shadowAlpha, m_colorOpacity[1] * m_shadowAlpha, m_colorOpacity[2] * m_shadowAlpha, 1.0);
  glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_COLOR);
#else
  glColor4f(1.0, 1.0, 1.0, m_shadowAlpha);
  //    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
#endif

  /*
     m_particleProg->enable();
     m_particleProg->setUniform1f("sliceNo", i);
     m_particleProg->setUniformfv("sortVector", &m_halfVector[0], 3, 1);
     m_particleProg->setUniform1f("numParticles", mNumParticles);
     m_particleProg->setUniform1f("numSlices", m_numSlices);
     m_particleProg->bindTexture("positionSampler", mPosBufferTexture, GL_TEXTURE_BUFFER_EXT, 3);
     */

  float sliceWidth = (m_maxDepth - m_minDepth) / (float) m_numSlices;
  float sliceZ = m_minDepth + (sliceWidth * i);
  //printf("%d: z = %f\n", i, sliceZ);
  m_particleAAProg->enable();
  m_particleAAProg->setUniform1f("sliceZ", sliceZ);
  m_particleAAProg->setUniform1f("invSliceWidth", 1.0f / sliceWidth);

  //drawPointSprites(m_particleProg, i*m_batchSize, m_batchSize, false);

  // render previous and current batch
  int start = (i-1)*m_batchSize;
  int end = start + m_batchSize*2;
  start = max(start, 0);
  end = min(end, (int) mNumParticles);
  drawPointSprites(m_particleAAProg, start, end - start, false, true);

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

// draw particles as slices with shadowing
void SmokeRenderer::drawSlices()
{
  m_batchSize = mNumParticles / m_numSlices;
  m_srcLightTexture = 0;

  setLightColor(m_lightColor);

  // clear light buffer
  m_fbo->Bind();
  glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
  m_fbo->AttachTexture(GL_TEXTURE_2D, m_lightTexture[m_srcLightTexture], GL_COLOR_ATTACHMENT0_EXT);
  m_fbo->AttachTexture(GL_TEXTURE_2D, 0, GL_COLOR_ATTACHMENT1_EXT);
  m_fbo->AttachTexture(GL_TEXTURE_2D, 0, GL_DEPTH_ATTACHMENT_EXT);
  //glClearColor(1.0 - m_lightColor[0], 1.0 - m_lightColor[1], 1.0 - m_lightColor[2], 0.0);
  glClearColor(m_lightColor[0], m_lightColor[1], m_lightColor[2], 0.0);
  //glClearColor(0.0f, 0.0f, 0.0f, 0.0f);	// clear to transparent
  glClear(GL_COLOR_BUFFER_BIT);

  // clear volume image
  m_fbo->AttachTexture(GL_TEXTURE_2D, m_imageTex[0], GL_COLOR_ATTACHMENT0_EXT);
  glClearColor(0.0, 0.0, 0.0, 0.0); 
  glClear(GL_COLOR_BUFFER_BIT);

#if 0
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  static float rotate = 0;
  glRotatef(rotate, 1,0,0);
  rotate += 0.1f;
  drawSkybox(m_cubemapTex);
  glPopMatrix();
#else
  //        drawSkybox(m_cubemapTex);
#endif

  /*
  // bind vbo as buffer texture
  glBindTexture(GL_TEXTURE_BUFFER_EXT, mPosBufferTexture);
  glTexBufferEXT(GL_TEXTURE_BUFFER_EXT, GL_RGBA32F_ARB, mPosVbo);
  */

#if USE_MRT
  // write to both color attachments
  m_fbo->AttachTexture(GL_TEXTURE_2D, m_imageTex[0], GL_COLOR_ATTACHMENT0_EXT);
  m_fbo->AttachTexture(GL_TEXTURE_2D, m_lightTexture[m_srcLightTexture], GL_COLOR_ATTACHMENT1_EXT);
  m_fbo->AttachTexture(GL_TEXTURE_2D, m_depthTex, GL_DEPTH_ATTACHMENT_EXT);

  GLenum buffers[] = { GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT };
  glDrawBuffers(2, buffers);
#endif

  glActiveTexture(GL_TEXTURE0);
  glMatrixMode(GL_TEXTURE);
  glLoadMatrixf((GLfloat *) m_shadowMatrix.get_value());

  // render slices
  if (m_numDisplayedSlices > m_numSlices) m_numDisplayedSlices = m_numSlices;

  for(int i=0; i<m_numDisplayedSlices; i++) {
#if 0
    // draw slice from camera view, sampling light buffer
    drawSlice(i);
    // draw slice from light view to light buffer, accumulating shadows
    drawSliceLightView(i);
    if (m_doBlur) {
      blurLightBuffer();
    }
#else
    // opposite order
    if (m_enableAA) {
      drawSliceLightViewAA(i);
    } else {
      drawSliceLightView(i);
    }
    if (m_doBlur) {
      blurLightBuffer();
    }

    drawSlice(i);
#endif
  }

#if USE_MRT
  glColorMaskIndexedEXT(0, true, true, true, true);
  glColorMaskIndexedEXT(1, false, false, false, false);
#endif
  m_fbo->Disable();

  glActiveTexture(GL_TEXTURE0);
  glMatrixMode(GL_TEXTURE);
  glLoadIdentity();
}


// blur light buffer to simulate scattering effects
void SmokeRenderer::blurLightBuffer()
{
  glViewport(0, 0, m_lightBufferSize, m_lightBufferSize);

  m_blurProg->enable();
  m_blurProg->setUniform2f("texelSize", 1.0f / (float) m_lightBufferSize, 1.0f / (float) m_lightBufferSize);
  glDisable(GL_DEPTH_TEST);

  for(int i=0; i<m_blurPasses; i++) {
#if 1
    // single pass
    m_fbo->AttachTexture(GL_TEXTURE_2D, m_lightTexture[1 - m_srcLightTexture], GL_COLOR_ATTACHMENT0_EXT);

    m_blurProg->bindTexture("tex", m_lightTexture[m_srcLightTexture], GL_TEXTURE_2D, 0);
    //m_blurProg->setUniform1f("blurRadius", m_blurRadius);
    m_blurProg->setUniform1f("blurRadius", m_blurRadius*(i+1));
    //m_blurProg->setUniform1f("blurRadius", m_blurRadius*powf(2.0f, (float) i));
    drawQuad();

    m_srcLightTexture = 1 - m_srcLightTexture;
#else
    // separable
    m_fbo->AttachTexture(GL_TEXTURE_2D, m_lightTexture[1 - m_srcLightTexture], GL_COLOR_ATTACHMENT0_EXT);
    m_blurProg->bindTexture("tex", m_lightTexture[m_srcLightTexture], GL_TEXTURE_2D, 0);
    //m_blurProg->setUniform1f("blurRadius", m_blurRadius);
    m_blurProg->setUniform1f("blurRadius", m_blurRadius*powf(2.0f, (float) i));
    m_blurProg->setUniform2f("texelSize", 1.0f / (float) m_lightBufferSize, 0.0f);
    drawQuad();
    m_srcLightTexture = 1 - m_srcLightTexture;

    m_fbo->AttachTexture(GL_TEXTURE_2D, m_lightTexture[1 - m_srcLightTexture], GL_COLOR_ATTACHMENT0_EXT);
    m_blurProg->bindTexture("tex", m_lightTexture[m_srcLightTexture], GL_TEXTURE_2D, 0);
    m_blurProg->setUniform2f("texelSize", 0.0f, 1.0f / (float) m_lightBufferSize);
    drawQuad();
    m_srcLightTexture = 1 - m_srcLightTexture;
#endif
  }
  m_blurProg->disable();
}

// post-process final volume image
void SmokeRenderer::processImage(GLSLProgram *prog, GLuint src, GLuint dest)
{
  m_fbo->Bind();
  m_fbo->AttachTexture(GL_TEXTURE_2D, dest, GL_COLOR_ATTACHMENT0_EXT);
  m_fbo->AttachTexture(GL_TEXTURE_2D, 0, GL_DEPTH_ATTACHMENT_EXT);

  prog->enable();
  prog->bindTexture("tex", src, GL_TEXTURE_2D, 0);

  glDisable(GL_DEPTH_TEST);
  drawQuad();

  prog->disable();
  m_fbo->Disable();
}

// display texture to screen
void SmokeRenderer::displayTexture(GLuint tex, float scale)
{
  m_displayTexProg->enable();
  m_displayTexProg->bindTexture("tex", tex, GL_TEXTURE_2D, 0);
  m_displayTexProg->setUniform1f("scale", scale);
  m_displayTexProg->setUniform1f("gamma", m_gamma);
  drawQuad();
  m_displayTexProg->disable();
}

#define DIAGONAL_STARS 0

void SmokeRenderer::doStarFilter()
{
  glViewport(0, 0, m_imageW, m_imageH);

#if 1
  // threshold
  m_thresholdProg->enable();
  m_thresholdProg->setUniform1f("scale", m_starPower);
  m_thresholdProg->setUniform1f("threshold", m_starThreshold);
  processImage(m_thresholdProg, m_imageTex[0], m_imageTex[3]);
#endif

  // star filter
  // horizontal
  m_starFilterProg->enable();
  m_starFilterProg->setUniform1f("radius", m_starBlurRadius);
#if DIAGONAL_STARS
  m_starFilterProg->setUniform2f("texelSize", 2.0f / (float) m_imageW, 2.0f / (float) m_imageH);	// diagonal
#else
  m_starFilterProg->setUniform2f("texelSize", 2.0f / (float) m_imageW, 0.0f);	// axis aligned
#endif
  m_starFilterProg->bindTexture("kernelTex", m_rainbowTex, GL_TEXTURE_2D, 1);
  processImage(m_starFilterProg, m_imageTex[3], m_imageTex[1]);
  //processImage(m_starFilterProg, m_imageTex[0], m_imageTex[1]);

  // vertical
  m_starFilterProg->enable();
#if DIAGONAL_STARS
  m_starFilterProg->setUniform2f("texelSize", -2.0f / (float) m_imageW, 2.0f / (float) m_imageH);	// diagonal
#else
  m_starFilterProg->setUniform2f("texelSize", 0.0f, 2.0f / (float) m_imageW);	// axis aligned
#endif
  processImage(m_starFilterProg, m_imageTex[3], m_imageTex[2]);
  //processImage(m_starFilterProg, m_imageTex[0], m_imageTex[2]);
}

void SmokeRenderer::downSample()
{
  // downsample
  glViewport(0, 0, m_downSampledW, m_downSampledH);
  m_downSampleProg->enable();
  //m_downSampleProg->setUniform2f("texelSize", 1.0f / (float) m_imageW, 1.0f / (float) m_imageH);
  processImage(m_downSampleProg, m_imageTex[0], m_downSampledTex[0]);
  m_downSampleProg->disable();
}

// anamorphic flare?
void SmokeRenderer::doFlare()
{
#if 1
  // threshold
  m_thresholdProg->enable();
  m_thresholdProg->setUniform1f("scale", 1.0f);
  m_thresholdProg->setUniform1f("threshold", m_flareThreshold);
  processImage(m_thresholdProg, m_downSampledTex[0], m_downSampledTex[1]);
#endif

#if 1
  m_gaussianBlurProg->enable();
  m_gaussianBlurProg->setUniform1f("radius", m_flareRadius);
  m_gaussianBlurProg->setUniform2f("texelSize", 2.0f / (float) m_downSampledW, 0.0f);
  //m_gaussianBlurProg->setUniform2f("texelSize", 1.0f / (float) m_downSampledW, 0.0f);
  processImage(m_gaussianBlurProg, m_downSampledTex[1], m_downSampledTex[2]);
#else
  m_starFilterProg->enable();
  m_starFilterProg->setUniform1f("radius", m_flareRadius);
  m_starFilterProg->setUniform2f("texelSize", 2.0f / (float) m_downSampledW, 0.0f);	// axis aligned
  m_starFilterProg->bindTexture("kernelTex", m_rainbowTex, GL_TEXTURE_2D, 1);
  processImage(m_starFilterProg, m_downSampledTex[1], m_downSampledTex[2]);
#endif
}

void SmokeRenderer::doGlowFilter()
{
  // blur
  m_gaussianBlurProg->enable();
  m_gaussianBlurProg->setUniform1f("radius", m_glowRadius);
  //m_gaussianBlurProg->setUniform2f("texelSize", 2.0f / (float) m_downSampledW, 0.0f);
  m_gaussianBlurProg->setUniform2f("texelSize", 1.0f / (float) m_downSampledW, 0.0f);
  processImage(m_gaussianBlurProg, m_downSampledTex[0], m_downSampledTex[1]);

  m_gaussianBlurProg->enable();
  //m_gaussianBlurProg->setUniform2f("texelSize", 0.0f, 2.0f / (float) m_downSampledH);
  m_gaussianBlurProg->setUniform2f("texelSize", 0.0f, 1.0f / (float) m_downSampledH);
  processImage(m_gaussianBlurProg, m_downSampledTex[1], m_downSampledTex[0]);
  m_gaussianBlurProg->disable();
}

// composite final volume image on top of scene
void SmokeRenderer::compositeResult()
{
  if (m_enableFilters) {
    if (m_starBlurRadius > 0.0f && m_starIntensity > 0.0f) {
      doStarFilter();
    }

    if (m_glowIntensity > 0.0f || m_flareIntensity > 0.0f) {
      downSample();
    }
    if (m_flareIntensity > 0.0f) {
      doFlare();
    }
    if (m_glowRadius > 0.0f && m_glowIntensity > 0.0f) {
      doGlowFilter();
    }
  }

  glViewport(0, 0, mWindowW, mWindowH);
  glDisable(GL_DEPTH_TEST);
  glDepthMask(GL_FALSE);
  glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
  glDisable(GL_BLEND);

  if (m_enableFilters) {
    m_compositeProg->enable();
    m_compositeProg->bindTexture("tex", m_imageTex[0], GL_TEXTURE_2D, 0);
    m_compositeProg->bindTexture("blurTexH", m_imageTex[1], GL_TEXTURE_2D, 1);
    m_compositeProg->bindTexture("blurTexV", m_imageTex[2], GL_TEXTURE_2D, 2);
    m_compositeProg->bindTexture("glowTex", m_downSampledTex[0], GL_TEXTURE_2D, 3);
    m_compositeProg->bindTexture("flareTex", m_downSampledTex[2], GL_TEXTURE_2D, 4);
    m_compositeProg->setUniform1f("scale", m_imageBrightness);
    m_compositeProg->setUniform1f("sourceIntensity", m_sourceIntensity);
    m_compositeProg->setUniform1f("glowIntensity", m_glowIntensity);
    m_compositeProg->setUniform1f("starIntensity", m_starIntensity);
    m_compositeProg->setUniform1f("flareIntensity", m_flareIntensity);
    m_compositeProg->setUniform1f("gamma", m_gamma);
    drawQuad();
    m_compositeProg->disable();
  } else {
    displayTexture(m_imageTex[0], m_imageBrightness);
    //displayTexture(m_downSampledTex[0], m_imageBrightness);
  }

  glDisable(GL_BLEND);
  glDepthMask(GL_TRUE);
}

void SmokeRenderer::drawBounds()
{
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  glPushMatrix();
  glTranslatef(0.0f, 0.0f, m_minDepth);
  glColor3f(0.0f, 1.0f, 0.0f);
  drawQuad();
  glPopMatrix();

  glPushMatrix();
  glTranslatef(0.0f, 0.0f, m_maxDepth);
  glColor3f(1.0f, 0.0f, 0.0f);
  drawQuad();
  glPopMatrix();

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glPopMatrix();
}

void SmokeRenderer::renderSprites(bool sort)
{
#if 1
  // post
  m_fbo->Bind();
  m_fbo->AttachTexture(GL_TEXTURE_2D, m_imageTex[0], GL_COLOR_ATTACHMENT0_EXT);
  m_fbo->AttachTexture(GL_TEXTURE_2D, 0, GL_DEPTH_ATTACHMENT_EXT);
  //m_fbo->AttachTexture(GL_TEXTURE_2D, m_depthTex, GL_DEPTH_ATTACHMENT_EXT);
  glViewport(0, 0, m_imageW, m_imageH);
  glClearColor(0.0, 0.0, 0.0, 0.0); 
  //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glClear(GL_COLOR_BUFFER_BIT);

  glDisable(GL_BLEND);
  //    drawSkybox(m_cubemapTex);
#endif

  glColor4f(1.0, 1.0, 1.0, m_spriteAlpha);
  if (sort) {
    calcVectors();
    depthSortCopy();
    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    drawPointSprites(m_particleProg, 0, mNumParticles, false, true);	
  } else {
    glBlendFunc(GL_ONE, GL_ONE);
    drawPointSprites(m_particleProg, 0, mNumParticles, false, false);
  }

#if 1
  m_fbo->Disable();
  compositeResult();
#endif
}

void SmokeRenderer::render()
{
  MPI_Bcast(this, sizeof(SmokeRendererParams),  MPI_BYTE, 0, MPI_COMM_WORLD);
#if 1
  switch(mDisplayMode) {
    case POINTS:
      glPointSize(2.0f);
      glEnable(GL_DEPTH_TEST);
      glColor4f(1.0, 1.0, 1.0, 1.0f);
      m_simpleProg->enable();
      drawPoints(0, mNumParticles, false);
      m_simpleProg->disable();
      glPointSize(1.0f);
      break;

#if 0
    case SPRITES:
      renderSprites(false);
      break;

    case SPRITES_SORTED:
      renderSprites(true);
      break;

    case VOLUMETRIC:
      calcVectors();
      depthSortCopy();
      drawSlices();
      compositeResult();
      //drawBounds();
      break;
#endif
    
    case SPLOTCH:
      splotchDraw();
      break;
    case SPLOTCH_SORTED:
      splotchDrawSort();
      break;

    case NUM_MODES:
      break;
  }
#endif

#if 0
  if (m_displayLightBuffer) {
    // display light buffer to screen
    glViewport(0, 0, m_lightBufferSize, m_lightBufferSize);
    glDisable(GL_DEPTH_TEST);
    displayTexture(m_lightTexture[m_srcLightTexture], 1.0f);
    glViewport(0, 0, mWindowW, mWindowH);
  }
#endif

  glutReportErrors();
}
    
static void lSetClippingPlane(const GLenum planeid, const float4 &plane)
{
  double eq[] = {plane.x, plane.y, plane.z, plane.w};
  glClipPlane(planeid, eq);
}

void lCompose(
    const float4* imgSrc,
    const float*  depthSrc,
    float4* imgDst,
    const int rank, const int nrank, const MPI_Comm &comm,
    const int2 imgCrd,
    const int2 imgSize,
    const int2 viewportSize,
    const int showDomain,
    const std::vector<int> &compositeOrder)
{
  constexpr int master = 0;

  using imgData_t = std::array<float,5>;
  constexpr int mpiImgDataSize = sizeof(imgData_t)/sizeof(float);
  static std::vector<imgData_t> sendbuf;

  /* copy img pixels ot send buffer */

  const int imgNPix = imgSize.x*imgSize.y;
  if (imgNPix > 0)
  {
    sendbuf.resize(imgNPix);
    if (compositeOrder.empty() && depthSrc)
    {
      /* if there is no global composition order and depth is != NULL */
#pragma omp parallel for schedule(static)
      for (int i = 0; i < imgNPix; i++)
        sendbuf[i] = imgData_t{{imgSrc[i].x, imgSrc[i].y, imgSrc[i].z, imgSrc[i].w, depthSrc[i]}};
    }
    else if (compositeOrder.empty())
    {
      /* if there is no global composition order and no depth is passed */
#pragma omp parallel for schedule(static)
      for (int i = 0; i < imgNPix; i++)
        sendbuf[i] = imgData_t{{imgSrc[i].x, imgSrc[i].y, imgSrc[i].z, imgSrc[i].w, static_cast<float>(rank)}};
    }
    else
    {
      /* if there is global composition order */
#pragma omp parallel for schedule(static)
      for (int i = 0; i < imgNPix; i++)
        sendbuf[i] = imgData_t{{imgSrc[i].x, imgSrc[i].y, imgSrc[i].z, imgSrc[i].w, static_cast<float>(compositeOrder[rank])}};
    }
  }

  /* compute which parts of img are sent to which rank */
  
  const int nPixels = viewportSize.x*viewportSize.y;
  const int nPixelsPerRank = (nPixels+nrank-1)/nrank; 

  const int x0 = imgCrd.x;
  const int y0 = imgCrd.y;
  const int x1 = imgCrd.x + imgSize.x;
  const int y1 = imgCrd.y + imgSize.y;

  const int w = viewportSize.x;

  const int imgBeg   =  y0   *w + x0;
  const int imgEnd   = (y1-1)*w + x1-1;
  const int imgWidth = imgSize.x;

  using imgMetaData_t = std::array<int,6>;
  constexpr  int mpiImgMetaDataSize = sizeof(imgMetaData_t)/sizeof(int);

  static std::vector<imgMetaData_t> srcMetaData(nrank);
  int totalSendCount = 0;
  for (int p = 0; p < nrank; p++)
  {
    /* domain scanline beginning & end */
    const int pbeg =    p * nPixelsPerRank;
    const int pend = pbeg + nPixelsPerRank-1;

    /* clip image with the domain scanline */
    const int clipBeg = std::min(pend, std::max(pbeg,imgBeg));
    const int clipEnd = std::max(pbeg, std::min(pend,imgEnd));

    int sendcount = 0;
    if (clipBeg < clipEnd)
    {
      const int i0 = clipBeg % w;
      const int j0 = clipBeg / w;
      const int i1 = clipEnd % w;
      const int j1 = clipEnd / w;
      assert(j0 >= y0);
      assert(j1 <  y1);

      /* compute number of pixels to send: 
       * multiply hight (j1-j0+1) by image width */
      /* but subtract top and bottom corners */

      const int dxtop = std::max(0,std::min(i0-x0,   imgWidth));
      const int dxbtm = std::max(0,std::min(x1-i1-1, imgWidth));
      sendcount = (j1-j0+1)*imgWidth - dxtop - dxbtm;
    }
    
    srcMetaData[p] = imgMetaData_t{{x0,y0,x1,y1,totalSendCount,sendcount}};
    totalSendCount += sendcount;
  }
  assert(totalSendCount == (x1-x0)*(y1-y0));

  /* exchange metadata info */

  static std::vector<imgMetaData_t> rcvMetaData(nrank);
  MPI_Alltoall(&srcMetaData[0], mpiImgMetaDataSize, MPI_INT, &rcvMetaData[0], mpiImgMetaDataSize, MPI_INT, comm);

  /* prepare counts & displacements for alltoallv */

  static std::vector<int> sendcount(nrank), senddispl(nrank+1);
  static std::vector<int> recvcount(nrank), recvdispl(nrank+1);
  senddispl[0] = recvdispl[0] = 0;
  for (int p = 0; p < nrank; p++)
  {
    sendcount[p  ] = srcMetaData[p][5] * mpiImgDataSize;
    recvcount[p  ] = rcvMetaData[p][5] * mpiImgDataSize;
    senddispl[p+1] = senddispl[p] + sendcount[p];
    recvdispl[p+1] = recvdispl[p] + recvcount[p];
  }

  static std::vector<imgData_t> recvbuf;
  {
    if (recvdispl[nrank] > 0)
      recvbuf.resize(recvdispl[nrank] / mpiImgDataSize);
    const double t0 = MPI_Wtime();
    MPI_Alltoallv(
        &sendbuf[0], &sendcount[0], &senddispl[0], MPI_FLOAT,
        &recvbuf[0], &recvcount[0], &recvdispl[0], MPI_FLOAT,
        comm);
    double nsendrecvloc = (senddispl[nrank] + recvdispl[nrank])*sizeof(float);
    double nsendrecv;
    MPI_Allreduce(&nsendrecvloc, &nsendrecv, 1, MPI_DOUBLE, MPI_SUM, comm);
    const double t1 = MPI_Wtime();
    if (rank == master)
    {
      const double dt = t1-t0;
      const double bw = nsendrecv / dt;
      fprintf(stderr, " MPI_Alltoallv: dt= %g  BW= %g MB/s  mem= %g MB\n", dt, bw/1e6, nsendrecv/1e6);
    }
  }

  /* pixel composition */


  constexpr int NRANKMAX = 1024;
  assert(nrank <= NRANKMAX);
    
  for (int p = 0; p < nrank+1; p++)
    recvdispl[p] /= mpiImgDataSize;

  static std::vector<float4> imgLoc;
  imgLoc.resize(nPixelsPerRank);
  const int pixelBeg =              rank * nPixelsPerRank;
  const int pixelEnd = std::min(pixelBeg + nPixelsPerRank, nPixels);
#pragma omp parallel for schedule(static)
  for (int idx = pixelBeg; idx < pixelEnd; idx++)
  {
    int pcount = 0;
    imgData_t imgData[NRANKMAX];

    const int i = idx % viewportSize.x;
    const int j = idx / viewportSize.x;

    for (int p = 0; p < nrank; p++)
      if (showDomain == -1 || showDomain == p)
      {
        const int x0   = rcvMetaData[p][0];
        const int y0   = rcvMetaData[p][1];
        const int x1   = rcvMetaData[p][2];
        const int y1   = rcvMetaData[p][3];
        const int offs = rcvMetaData[p][4];
        const int cnt  = rcvMetaData[p][5];
        const int idx =  (j-y0)*(x1-x0) + (i-x0) - offs;
        if (x0  <= i && i   < x1 &&
            y0  <= j && j   < y1 && 
            idx >= 0 && idx < cnt)
          imgData[pcount++] = recvbuf[recvdispl[p] + idx];
      }

    std::sort(imgData, imgData+pcount, 
        [](const imgData_t &a, const imgData_t &b) { return a[4] < b[4]; });

    float4 dst = make_float4(0.0f);
    for (int p = 0; p < pcount; p++)
    {
      auto &src = imgData[p];
      src[0] *= 1.0f - dst.w;
      src[1] *= 1.0f - dst.w;
      src[2] *= 1.0f - dst.w;
      src[3] *= 1.0f - dst.w;

      dst.x += src[0];
      dst.y += src[1];
      dst.z += src[2];
      dst.w += src[3];

      dst.w = std::min(dst.w, 1.0f);
    }
    imgLoc[idx - pixelBeg] = dst;
  }

  /* gather composited part of images into a single image on the master rank */
  {
    const double t0 = MPI_Wtime();
    MPI_Gather(&imgLoc[0], nPixelsPerRank*4, MPI_FLOAT, imgDst, 4*nPixelsPerRank, MPI_FLOAT, master, comm);
    const double t1 = MPI_Wtime();
    if (master == rank)
    {
      const double dt        = t1 - t0;
      const double nsendrecv = nPixelsPerRank*4*nrank*sizeof(float);
      const double bw        = nsendrecv / dt;
      fprintf(stderr, " MPI_Gather: dt= %g  BW= %g MB/s  mem= %g MB\n", dt, bw/1e6, nsendrecv/1e6);
    }
  }
}

#if 0
static void lCompose(
    float4 *src, float4 *dst, float *depth,
    const int rank, const int nrank, const MPI_Comm &comm,
    const int2 wCrd, const int2 wSize,
    const int2 viewPort, const bool resize = true)
{
  constexpr int master = 0;

#if 0  /* could be NULL if nothing is drawn */
  assert(src   != NULL);
  assert(dst   != NULL);
  assert(depth != NULL);
#endif

  assert(wCrd.x >= 0);
  assert(wCrd.y >= 0);
  assert(wCrd.x + wSize.x <= viewPort.x);
  assert(wCrd.y + wSize.y <= viewPort.y);

  static int4 localTile;           /* (xmin,ymin,xmax,ymax) local compositing tile */
  static std::vector<float4> imgLoc, imgGlb;
  static int npx = -1, npy;

  if (resize)
    npx = -1;

  if (npx == -1)
  {
    /* factorize nrank = npx*npy */
    int n_tmp = static_cast<int>(std::sqrt(nrank+0.1));
    while(nrank%n_tmp)
      n_tmp--;
    npy = n_tmp;  
    npx = nrank/npy;
    assert(npx*npy == nrank);

    const int irank = rank % npx;
    const int jrank = rank / npx;
    assert(irank < npx);
    assert(jrank < npy);

    const int winxloc = (viewPort.x + npx - 1) / npx;
    const int winyloc = (viewPort.y + npy - 1) / npy;
    
    localTile.x = irank * winxloc;
    localTile.y = jrank * winyloc;
    localTile.z = winxloc;
    localTile.w = winyloc;

    /* allocate buffers */
    imgLoc.resize(winxloc*winyloc);
    imgGlb.resize(winxloc*winyloc*nrank);
    assert((int)imgGlb.size() >= viewPort.x*viewPort.y);
  }
  
  using vec5 = std::array<float,5>;  /* rgb, alpha, depth */
  constexpr int mpiDataSize = sizeof(vec5)/sizeof(float);

  /************************************/
  /* sort pixels by destination ranks */
  /************************************/

  static std::vector<int4> tilesBnd(nrank); 
  static std::vector<int> sendcount(nrank), senddispl(nrank+1);
  senddispl[0] = 0;

#pragma omp parallel for schedule(static)
  for (int p = 0; p < nrank; p++)
  {
    const int irank = p % npx;
    const int jrank = p / npx;

    const int winxloc = localTile.z;
    const int winyloc = localTile.w;

    const int x0 = irank*winxloc;
    const int x1 = x0  + winxloc;

    const int y0 = jrank*winyloc;
    const int y1 = y0  + winyloc;

    using vec4 = std::array<int,4> ;

    /* clip image with the tile */
    auto clip = [](const vec4 &tile, const vec4 &image)
    {
      vec4 t;
      t[0] = std::min(tile[2], std::max(tile[0],image[0]));  /* clip xmin */
      t[1] = std::min(tile[3], std::max(tile[1],image[1]));  /* clip ymin */
      t[2] = std::max(tile[0], std::min(tile[2],image[2]));  /* clip xmax */
      t[3] = std::max(tile[1], std::min(tile[3],image[3]));  /* clip ymax */
      if (t[0] >= t[2] || t[1] >= t[3])  /* image is fully outside the tile */
        t = vec4{{0,0,0,0}};
      return t;
    };

    const auto &tile = clip(
        vec4{{x0,y0,x1,y1}},
        vec4{{wCrd.x,wCrd.y, wCrd.x+wSize.x, wCrd.y+wSize.y}});

    tilesBnd [p]   = make_int4(tile[0], tile[1], tile[2]-tile[0], tile[3]-tile[1]);
    sendcount[p]   = tilesBnd[p].z*tilesBnd[p].w*mpiDataSize;
  }

  for (int p = 0; p < nrank; p++)
    senddispl[p+1] = senddispl[p] + sendcount[p];

  static std::vector<vec5> sendbuf;
  if (senddispl[nrank] > 0)
    sendbuf.resize(senddispl[nrank] / mpiDataSize);

  /************************/
  /* populate send buffer */
  /************************/

#pragma omp parallel
  {
#pragma omp for schedule(dynamic)
    for (int p = 0; p < nrank; p++)
      if (sendcount[p] > 0)
      {
        const int xmin  = tilesBnd[p].x;
        const int ymin  = tilesBnd[p].y;
        const int xmax  = tilesBnd[p].z + xmin;
        const int ymax  = tilesBnd[p].w + ymin;
        const int displ = senddispl[p] / mpiDataSize;
//#pragma omp for schedule(static) collapse(2) nowait
        for (int j = ymin; j < ymax; j++)
          for (int i = xmin; i < xmax; i++)
          {
            const int iloc = i - wCrd.x;
            const int jloc = j - wCrd.y;
            assert(iloc >= 0); 
            assert(jloc >= 0);
            assert(iloc < wSize.x);
            assert(jloc < wSize.y);

            const int idx = jloc*wSize.x + iloc;
            sendbuf[displ + (j-ymin)*(xmax-xmin)+(i-xmin)] = vec5{{
              src[idx].x,src[idx].y,src[idx].z, src[idx].w, 
                depth[idx]}};
          }
      }
  }
  
  /***************************/
  /* exchange tiles metadata */
  /***************************/

  static std::vector<int4> tilesBndRecv(nrank);
  MPI_Alltoall(&tilesBnd[0], 4, MPI_INT, &tilesBndRecv[0], 4, MPI_INT, comm);

  static std::vector<int> recvcount(nrank), recvdispl(nrank+1);
  for (int p = 0; p < nrank; p++)
    recvcount[p] = tilesBndRecv[p].z*tilesBndRecv[p].w * mpiDataSize;
     
  recvdispl[0] = 0; 
  for (int p = 0; p < nrank; p++)
    recvdispl[p+1] = recvdispl[p] + recvcount[p];

  static std::vector<vec5> recvbuf;
  if (recvdispl[nrank] > 0)
    recvbuf.resize(recvdispl[nrank] / mpiDataSize);

  /***********************/
  /* exchange pixel data */
  /***********************/
  
  {
//    MPI_Barrier(comm);
    const double t0 = MPI_Wtime();
    MPI_Alltoallv(
        &sendbuf[0], &sendcount[0], &senddispl[0], MPI_FLOAT,
        &recvbuf[0], &recvcount[0], &recvdispl[0], MPI_FLOAT,
        comm);
    double nsendrecvloc = recvdispl[nrank] + senddispl[rank];
    double nsendrecv;
    MPI_Allreduce(&nsendrecvloc, &nsendrecv, 1, MPI_DOUBLE, MPI_SUM, comm);
    const double t1 = MPI_Wtime();
    if (rank == master)
    {
      const double dt = t1-t0;
      const double bw = sizeof(float)*nsendrecv / dt;
      fprintf(stderr, " MPI_Alltoallv: dt= %g  BW= %g MB/s  mem= %g MB\n", dt, bw/1e6, sizeof(float)*nsendrecv/1e6);
    }
  }
  
  for (int p = 0; p < nrank; p++)
  {
    recvcount[p] /= mpiDataSize;
    recvdispl[p] /= mpiDataSize;
  }
  recvdispl[nrank] /= mpiDataSize;
  

  /********************/
  /* compositing step */
  /********************/

#pragma omp parallel
  {
    const int xmin = localTile.x;
    const int ymin = localTile.y;
    const int xmax = localTile.z + xmin;
    const int ymax = localTile.w + ymin;

    std::vector<vec5> colors;
    colors.reserve(1024);
#pragma omp for schedule(static) collapse(2)
    for (int j = ymin; j < ymax; j++)
      for (int i = xmin; i < xmax; i++)
      {
        colors.clear();

        /* extract non-zero pixels */
        for (int p = 0; p < nrank; p++)
        {
          const int xminr = tilesBndRecv[p].x;
          const int yminr = tilesBndRecv[p].y;
          const int xmaxr = tilesBndRecv[p].z + xminr;
          const int ymaxr = tilesBndRecv[p].w + yminr;
          if (i >= xminr && i < xmaxr && j >= yminr && j < ymaxr)
          {
            const int iloc = i - xminr;
            const int jloc = j - yminr;
            assert(iloc >= 0);
            assert(jloc >= 0);
            assert(iloc <  localTile.z);
            assert(jloc <  localTile.w);
            const int idx = jloc*(xmaxr - xminr) + iloc;
            assert(recvdispl[p] + idx < recvdispl[p+1]);
            colors.push_back(recvbuf[recvdispl[p] + idx]);
          }
        }

        /* sort by depth */
        std::sort(
            colors.begin(), colors.end(),
            [](const vec5 &a, const vec5 &b){ return a[4] < b[4]; }
            );

        /* compose */
        float4 dst = make_float4(0.0f);
        for (auto &src : colors)
        {
          src[0] *= 1.0f - dst.w;
          src[1] *= 1.0f - dst.w;
          src[2] *= 1.0f - dst.w;
          src[3] *= 1.0f - dst.w;

          dst.x += src[0];
          dst.y += src[1];
          dst.z += src[2];
          dst.w += src[3];

          dst.w = std::min(dst.w, 1.0f);
        }
        imgLoc[(j-ymin)*(xmax-xmin)+(i-xmin)] = dst;
      }
  }
  
  /******************************************/
  /* gather local images on the master rank */
  /******************************************/

  MPI_Gather(
      &imgLoc[0], imgLoc.size()*4, MPI_FLOAT, 
      &imgGlb[0], imgLoc.size()*4, MPI_FLOAT, 
      master, comm);
  
  /* combine times together */
  const int winxloc = localTile.z;
  const int winyloc = localTile.w;
  assert((int)imgLoc.size() == winxloc*winyloc);

#pragma omp parallel for schedule(static) collapse(2)
  for (int j = 0; j < viewPort.y; j++)
    for (int i = 0; i < viewPort.x; i++)
    {
      const int iloc  = i % winxloc;
      const int jloc  = j % winyloc;

      const int irank = i / winxloc;
      const int jrank = j / winyloc;

      const int rank  = jrank * npx + irank;
      dst[j*viewPort.x+i] = imgGlb[winxloc*winyloc*rank + jloc*winxloc + iloc];
    }
  
}
#endif

static void lCompose(
    float4 *src, float4 *dst, float *depth,
    const int n, const int rank, const int nrank, const MPI_Comm &comm,
    const int showDomain,
    std::vector<int> compositingOrder)
{
  const int master = 0;
#if 0
  MPI_Reduce(src, dst, 4*n, MPI_FLOAT, MPI_SUM, master, comm);
#else
  const int nsend = (n+nrank-1)/nrank;
  static std::vector<float4> colorArray;
  colorArray.resize(nsend*nrank);

  if (!depth)
  {
#pragma omp parallel for schedule(static)
    for (int i = n; i < nsend*nrank; i++)
      src[i] = make_float4(0.0f);

    MPI_Alltoall(src, nsend*4, MPI_FLOAT, &colorArray[0], nsend*4, MPI_FLOAT, comm);

    assert(showDomain >= -1 && showDomain < nrank);
    if (showDomain == -1)
    {
#pragma omp parallel for schedule(static)
      for (int i = 0; i < nsend; i++)
        for (int p = 1; p < nrank; p++)
        {
          float4 dst = colorArray[i];
          float4 src = colorArray[i + p*nsend];
          dst.x += src.x;
          dst.y += src.y;
          dst.z += src.z;
          dst.w += src.w;
          colorArray[i] = dst;
        }
    }
    else
      std::copy(
          colorArray.begin() + showDomain*nsend,
          colorArray.begin() + showDomain*nsend + nsend,
          colorArray.begin());
  }
  else
  {
    static std::vector<float > depthArray;
    depthArray.resize(nsend*nrank);
#pragma omp parallel for schedule(static)
    for (int i = n; i < nsend*nrank; i++)
    {
      src[i] = make_float4(0.0f);
      depth[i] = 1.0f;
    }

    MPI_Alltoall(src, nsend*4, MPI_FLOAT, &colorArray[0], nsend*4, MPI_FLOAT, comm);
    MPI_Alltoall(depth, nsend, MPI_FLOAT, &depthArray[0], nsend, MPI_FLOAT, comm);

    using vec5 = std::array<float,5>;
    std::vector<vec5> colorArrayDepth(nsend*nrank);

    const bool doPixelLevel = compositingOrder.empty();

    if (doPixelLevel)
    {
      compositingOrder.resize(nrank);
      std::iota(compositingOrder.begin(), compositingOrder.end(), 0);
    }


#pragma omp parallel for schedule(static)
    for (int i = 0; i < nsend; i++)
      for (int p = 0; p < nrank; p++)
      {
        const int dst = i*nrank+p;
        const int src = p*nsend+i;
        colorArrayDepth[dst] = 
          vec5{{
            colorArray[src].x,colorArray[src].y,colorArray[src].z,colorArray[src].w,
            depthArray[src]
          }};
      }

   
   if (showDomain == -1)
   { 
#pragma omp parallel for schedule(static)
     for (int i = 0; i < nsend; i++)
     {
       const int stride = i*nrank;
       if (doPixelLevel)
       {
         std::sort(
             colorArrayDepth.begin() + stride, 
             colorArrayDepth.begin() + stride + nrank,
             [](const vec5 &a, const vec5 &b){ return a[4] < b[4]; }
             );
       }

       float4 _dst = make_float4(0.0);
       for (auto p : compositingOrder)
       {
         float4 _src = make_float4(
             colorArrayDepth[stride+p][0],
             colorArrayDepth[stride+p][1],
             colorArrayDepth[stride+p][2],
             colorArrayDepth[stride+p][3]
             );

         _src.x *= 1.0f - _dst.w;
         _src.y *= 1.0f - _dst.w;
         _src.z *= 1.0f - _dst.w;
         _src.w *= 1.0f - _dst.w;

         _dst.x += _src.x;
         _dst.y += _src.y;
         _dst.z += _src.z;
         _dst.w += _src.w;
         _dst.w = std::min(1.0f, _dst.w);
         assert(_dst.w >= 0.0f);
       }
       colorArray[i] = _dst;
     }
   }
   else
     std::copy(
         colorArray.begin() + showDomain*nsend,
         colorArray.begin() + showDomain*nsend + nsend,
         colorArray.begin());
  }

  MPI_Gather(&colorArray[0], nsend*4, MPI_FLOAT, dst, 4*nsend, MPI_FLOAT, master, comm);
#endif
}

std::array<int,4> SmokeRenderer::getVisibleViewport() const
{
  const float3 r0 = m_xlow;
  const float3 r1 = m_xhigh;
  const double3 bBoxVtx[] = {
    make_double3(r0.x,r0.y,r0.z),
    make_double3(r1.x,r0.y,r0.z),
    make_double3(r0.x,r1.y,r0.z),
    make_double3(r1.x,r1.y,r0.z),
    make_double3(r0.x,r0.y,r1.z),
    make_double3(r1.x,r0.y,r1.z),
    make_double3(r0.x,r1.y,r1.z),
    make_double3(r1.x,r1.y,r1.z)
  };

  int viewport[4];
  glGetIntegerv( GL_VIEWPORT, viewport);

  const int w = viewport[2];
  const int h = viewport[3];

  std::array<int,4> visibleViewport{{w,h,0,0}};

  for (auto v : bBoxVtx)
  {
#if 1
    double x,y,z;
    gluProject(v.x,v.y,v.z, m_modelViewWin,m_projectionWin,viewport,&x,&y,&z);
#else
    const double4 pos0 = make_double4(v.x,v.y,v.z,1.0);
    const double4 posO = lMatVec(m_modelViewWin,pos0);
    const double4 posP = lMatVec(m_projectionWin,posO);

    const double wclip = 1.0/posP.w;
    double3 posV = make_double3(posP.x*wclip, posP.y*wclip, posP.z*wclip);

#if 0
    fprintf(stderr, " rank= %d: xyz= %g %g %g  pxyz= %g %g %g \n",
        rank, v.x,v.y,v.z, posV.x,posV.y,posV.z);
#endif

    double x = (posV.x + 1.0)*0.5*w;
    double y = (posV.y + 1.0)*0.5*h;
#endif


    x = std::max(x, static_cast<double>(0));
    y = std::max(y, static_cast<double>(0));
    x = std::min(x, static_cast<double>(w));
    y = std::min(y, static_cast<double>(h));

    visibleViewport[0] = std::min(visibleViewport[0], static_cast<int>(floor(x)));
    visibleViewport[1] = std::min(visibleViewport[1], static_cast<int>(floor(y)));
    visibleViewport[2] = std::max(visibleViewport[2], static_cast<int>(ceil(x)));
    visibleViewport[3] = std::max(visibleViewport[3], static_cast<int>(ceil(y)));
    
  }
#if 0
  fprintf(stderr, "rank= %d:  %d %d  - %d %d - %d %d \n",
      rank, 
      visibleViewport[0],
      visibleViewport[1],
      visibleViewport[2],
      visibleViewport[3],
      w,h);
#endif
#if 0
  visibleViewport[0]= std::max(visibleViewport[0],0);
  visibleViewport[1]= std::max(visibleViewport[1],0);
  visibleViewport[2]= std::min(visibleViewport[2],w);
  visibleViewport[3]= std::min(visibleViewport[3],h);
#endif


  visibleViewport[2] -= visibleViewport[0];
  visibleViewport[3] -= visibleViewport[1];

  return visibleViewport;
}

void SmokeRenderer::splotchDraw()
{
  m_fbo->Bind();
  m_fbo->AttachTexture(GL_TEXTURE_2D, m_imageTex[0], GL_COLOR_ATTACHMENT0_EXT);
  m_fbo->AttachTexture(GL_TEXTURE_2D, 0, GL_DEPTH_ATTACHMENT_EXT);
  glViewport(0, 0, m_imageW, m_imageH);
  glClearColor(0.0, 0.0, 0.0, 0.0); 
  glClear(GL_COLOR_BUFFER_BIT);
  glClear(GL_DEPTH_BUFFER_BIT);
  glDisable(GL_BLEND);


  const int start = 0;
  const int count = mNumParticles;

  calcVectors();
  glBlendFunc(GL_ONE, GL_ONE);

  glDisable(GL_DEPTH_TEST);
  glDepthMask(GL_FALSE);  // don't write depth
  glEnable(GL_BLEND);

  auto &prog = m_splotchProg;

  GLuint vertexLoc = -1;
  if (!mSizeVao && mSizeVbo)
  {
    glGenVertexArrays(1, &mSizeVao);
    glBindVertexArray(mSizeVao);
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, mSizeVbo);
    vertexLoc = prog->getAttribLoc("particleSize");
    glEnableVertexAttribArray(vertexLoc);
    glVertexAttribPointer(vertexLoc , 1, GL_FLOAT, 0, 0, 0);
  }

  if (m_doClipping)
  {
#ifdef _SPLOTCHSPRITES 
    lSetClippingPlane(GL_CLIP_PLANE0, m_clippingPlane[0]);
    lSetClippingPlane(GL_CLIP_PLANE1, m_clippingPlane[1]);
    lSetClippingPlane(GL_CLIP_PLANE2, m_clippingPlane[2]);
    lSetClippingPlane(GL_CLIP_PLANE3, m_clippingPlane[3]);
    lSetClippingPlane(GL_CLIP_PLANE4, m_clippingPlane[4]);
    lSetClippingPlane(GL_CLIP_PLANE5, m_clippingPlane[5]);
#else
    glEnable(GL_CLIP_DISTANCE0);
    glEnable(GL_CLIP_DISTANCE1);
    glEnable(GL_CLIP_DISTANCE2);
    glEnable(GL_CLIP_DISTANCE3);
    glEnable(GL_CLIP_DISTANCE4);
    glEnable(GL_CLIP_DISTANCE5);
#endif
  }

  prog->enable();
  glBindVertexArray(mSizeVao);

  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  // VS
  prog->setUniform1f("spriteScale", viewport[3] / mInvFocalLen);
  prog->setUniform1f("starScale", powf(10.0f, m_starScaleLog));
  prog->setUniform1f("starAlpha", m_starAlpha);
  prog->setUniform1f("dmScale",  powf(10.0f, m_dmScaleLog));
  prog->setUniform1f("dmAlpha",  m_dmAlpha);
  prog->setUniform1f("spriteSizeMax", 5.0*powf(10.0f, m_spriteSizeMaxLog));
  // PS
  prog->bindTexture("spriteTex",  m_sphTex, GL_TEXTURE_2D, 1);
  prog->setUniform1f("alphaScale", m_spriteAlpha);
  prog->setUniform1f("transmission", m_transmission);
  prog->setUniform1f("resx", m_imageW);
  prog->setUniform1f("resy", m_imageH);

  prog->setUniformfv("p0o", (GLfloat*)&m_clippingPlane[0], 4, 1);
  prog->setUniformfv("p1o", (GLfloat*)&m_clippingPlane[1], 4, 1);
  prog->setUniformfv("p2o", (GLfloat*)&m_clippingPlane[2], 4, 1);
  prog->setUniformfv("p3o", (GLfloat*)&m_clippingPlane[3], 4, 1);
  prog->setUniformfv("p4o", (GLfloat*)&m_clippingPlane[4], 4, 1);
  prog->setUniformfv("p5o", (GLfloat*)&m_clippingPlane[5], 4, 1);

  prog->setUniform1f("sorted", 0);

  //glClientActiveTexture(GL_TEXTURE0);
  glActiveTexture(GL_TEXTURE0);
  glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);
  glEnable(GL_POINT_SPRITE_ARB);

  drawPoints(start,count,false);

  prog->disable();
    
#ifdef _SPLOTCHSPRITES 
  glDisable(GL_CLIP_PLANE0);
  glDisable(GL_CLIP_PLANE1);
  glDisable(GL_CLIP_PLANE2);
  glDisable(GL_CLIP_PLANE3);
  glDisable(GL_CLIP_PLANE4);
  glDisable(GL_CLIP_PLANE5);
#else
  glDisable(GL_CLIP_DISTANCE0);
  glDisable(GL_CLIP_DISTANCE1);
  glDisable(GL_CLIP_DISTANCE2);
  glDisable(GL_CLIP_DISTANCE3);
  glDisable(GL_CLIP_DISTANCE4);
  glDisable(GL_CLIP_DISTANCE5);
#endif

#if 1
  glFlush();
  glFinish();
#endif

#if 1
  {
    const double t0 = MPI_Wtime();
    GLint w, h, internalformat;

    glBindTexture(GL_TEXTURE_2D, m_imageTex[0]);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH,  &w);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_INTERNAL_FORMAT, &internalformat);
    glBindTexture(GL_TEXTURE_2D,0);

    static std::vector<float4> imgLoc, imgGlb;
    imgLoc.resize(2*w*h);
    imgGlb.resize(2*w*h);
    const int imgSize = w*h*4*sizeof(float);

    const double t1 = MPI_Wtime();

    static GLuint pbo_id[2];
    if (!pbo_id[0])
    {
      const int pbo_size = 8*1920*1080*4*sizeof(float);
      glGenBuffers(2, pbo_id);
      glBindBuffer(GL_PIXEL_PACK_BUFFER, pbo_id[0]);
      glBufferData(GL_PIXEL_PACK_BUFFER, pbo_size, 0, GL_STATIC_READ);
      glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
      glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pbo_id[1]);
      glBufferData(GL_PIXEL_UNPACK_BUFFER, pbo_size, 0, GL_STATIC_DRAW);
      glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
    }
    assert(pbo_id[0] && pbo_id[1]);

    /** fetch image ***/
    glBindTexture(GL_TEXTURE_2D, m_imageTex[0]);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, pbo_id[0]);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_FLOAT, 0);
    //      glReadPixels(0, 0, w, h, GL_RGBA, GL_FLOAT, 0);
    glFinish();
    const double t2 = MPI_Wtime();

    GLvoid *rptr = glMapBufferRange(GL_PIXEL_PACK_BUFFER, 0, imgSize, GL_MAP_READ_BIT);

#pragma omp parallel for schedule(static)
    for (int i = 0; i < w*h; i++)
      imgLoc[i] = reinterpret_cast<float4*>(rptr)[i];

    //      memcpy(&imgLoc[0], rptr, imgSize);
    glUnmapBuffer(GL_PIXEL_PACK_BUFFER);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    glBindTexture(GL_TEXTURE_2D,0);


    glFinish();
    const double t3 = MPI_Wtime();

    lCompose(&imgLoc[0], &imgGlb[0], NULL, w*h, rank, nrank, comm,
        m_domainView ? m_domainViewIdx : -1,
        compositingOrder);
    glFinish();
    const double t4 = MPI_Wtime();

    if (isMaster())
    {
      glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pbo_id[1]);
      GLvoid *wptr = glMapBufferRange(GL_PIXEL_UNPACK_BUFFER, 0, imgSize, GL_MAP_WRITE_BIT);

#pragma omp parallel for schedule(static)
      for (int i = 0; i < w*h; i++)
        reinterpret_cast<float4*>(wptr)[i] = imgGlb[i];
      glFinish();
      const double t5 = MPI_Wtime();

      glUnmapBuffer(GL_PIXEL_UNPACK_BUFFER);

      glBindTexture(GL_TEXTURE_2D, m_imageTex[0]);
      glTexImage2D(GL_TEXTURE_2D, 0, internalformat, w,h,0,GL_RGBA,GL_FLOAT, 0);
      glFinish();
      glBindTexture(GL_TEXTURE_2D,0);
      glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
      const double t6 = MPI_Wtime();

      if (1)
        fprintf(stderr, 
            "total= %g: d2h= %g cpy= %g  mpi= %g  cpy= %g h2d= %g :: bwMPI= %g bwD2H= %g  bwH2D= %g\n", t6-t0,
            t2-t1,   t3-t2,       t4-t3,   t5-t4,     t6-t5,
            3.0*imgSize/(t4-t3)/1e6, imgSize/(t2-t1)/1e6, imgSize/(t6-t5)/1e6);
    }

  }
#endif

  m_fbo->Disable();



#if 1 
  glDisable(GL_BLEND);
  m_splotch2texProg->enable();
  m_splotch2texProg->bindTexture("tex", m_imageTex[0], GL_TEXTURE_2D, 0);
  m_splotch2texProg->setUniform1f("scale_pre", 0.1*m_imageBrightnessPre);
  m_splotch2texProg->setUniform1f("gamma_pre", m_gammaPre);
  m_splotch2texProg->setUniform1f("scale_post", m_imageBrightnessPost);
  m_splotch2texProg->setUniform1f("gamma_post", m_gammaPost);
  m_splotch2texProg->setUniform1f("sorted", 0);
  drawQuad();
  m_splotch2texProg->disable();
#else
#error
  {
#if 0
    m_fbo->Bind();
    m_fbo->AttachTexture(GL_TEXTURE_2D, m_imageTex[0], GL_COLOR_ATTACHMENT0_EXT);
    m_fbo->AttachTexture(GL_TEXTURE_2D, 0, GL_DEPTH_ATTACHMENT_EXT);
    glViewport(0, 0, m_imageW, m_imageH);
    glClearColor(0.0, 0.0, 0.0, 0.0); 
    glClear(GL_COLOR_BUFFER_BIT);
    /* run prog */
    m_fbo->Disable();
#endif

    if (m_enableFilters) {
      if (m_starBlurRadius > 0.0f && m_starIntensity > 0.0f) 
      {
        doStarFilter();
      }

      if (m_glowIntensity > 0.0f || m_flareIntensity > 0.0f) {
        downSample();
      }
      if (m_flareIntensity > 0.0f) {
        doFlare();
      }
      if (m_glowRadius > 0.0f && m_glowIntensity > 0.0f) {
        doGlowFilter();
      }
    }

    glViewport(0, 0, mWindowW, mWindowH);
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_BLEND);

    if (m_enableFilters) {
      m_compositeProg->enable();
      m_compositeProg->bindTexture("tex", m_imageTex[0], GL_TEXTURE_2D, 0);
      m_compositeProg->bindTexture("blurTexH", m_imageTex[1], GL_TEXTURE_2D, 1);
      m_compositeProg->bindTexture("blurTexV", m_imageTex[2], GL_TEXTURE_2D, 2);
      m_compositeProg->bindTexture("glowTex", m_downSampledTex[0], GL_TEXTURE_2D, 3);
      m_compositeProg->bindTexture("flareTex", m_downSampledTex[2], GL_TEXTURE_2D, 4);
      m_compositeProg->setUniform1f("scale", m_imageBrightness);
      m_compositeProg->setUniform1f("sourceIntensity", m_sourceIntensity);
      m_compositeProg->setUniform1f("glowIntensity", m_glowIntensity);
      m_compositeProg->setUniform1f("starIntensity", m_starIntensity);
      m_compositeProg->setUniform1f("flareIntensity", m_flareIntensity);
      m_compositeProg->setUniform1f("gamma", m_gamma);
      drawQuad();
      m_compositeProg->disable();
    } else {
      displayTexture(m_imageTex[0], m_imageBrightness);
      //displayTexture(m_downSampledTex[0], m_imageBrightness);
    }

    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);
  }
#endif
}

#if 0
float4 lPlaneEquation(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3)
{
  float4 eq;
  eq.x = (y1*(z2 – z3)) + (y2*(z3 – z1)) + (y3*(z1 – z2));
  eq.x = (z1*(x2 – x3)) + (z2*(x3 – x1)) + (z3*(x1 – x2));
  eq.z = (x1*(y2 – y3)) + (x2*(y3 – y1)) + (x3*(y1 – y2));
  eq.w = -((x1*((y2*z3) – (y3*z2))) + (x2*((y3*z1) – (y1*z3))) + (x3*((y1*z2) – (y2*z1))));
  return eq;
}
#endif

void SmokeRenderer::splotchDrawSort()
{
  m_fbo->Bind();
  m_fbo->AttachTexture(GL_TEXTURE_2D, m_imageTex[0], GL_COLOR_ATTACHMENT0_EXT);
  m_fbo->AttachTexture(GL_TEXTURE_2D, m_depthTex, GL_DEPTH_ATTACHMENT_EXT);
  glViewport(0, 0, m_imageW, m_imageH);
  glClearColor(0.0, 0.0, 0.0, 0.0); 
  glClearDepth(1.0f);
  glDepthMask(GL_TRUE);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glDisable(GL_BLEND);


  const int start = 0;
  const int count = mNumParticles;

  calcVectors();
  depthSortCopy();
  glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);   

  auto &prog = m_splotchProg;

  GLuint vertexLoc = -1;
  if (!mSizeVao && mSizeVbo)
  {
    glGenVertexArrays(1, &mSizeVao);
    glBindVertexArray(mSizeVao);
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, mSizeVbo);
    vertexLoc = prog->getAttribLoc("particleSize");
    glEnableVertexAttribArray(vertexLoc);
    glVertexAttribPointer(vertexLoc , 1, GL_FLOAT, 0, 0, 0);
  }

  if (m_doClipping)
  {
#ifdef _SPLOTCHSPRITES 
    lSetClippingPlane(GL_CLIP_PLANE0, m_clippingPlane[0]);
    lSetClippingPlane(GL_CLIP_PLANE1, m_clippingPlane[1]);
    lSetClippingPlane(GL_CLIP_PLANE2, m_clippingPlane[2]);
    lSetClippingPlane(GL_CLIP_PLANE3, m_clippingPlane[3]);
    lSetClippingPlane(GL_CLIP_PLANE4, m_clippingPlane[4]);
    lSetClippingPlane(GL_CLIP_PLANE5, m_clippingPlane[5]);
#else
    glEnable(GL_CLIP_DISTANCE0);
    glEnable(GL_CLIP_DISTANCE1);
    glEnable(GL_CLIP_DISTANCE2);
    glEnable(GL_CLIP_DISTANCE3);
    glEnable(GL_CLIP_DISTANCE4);
    glEnable(GL_CLIP_DISTANCE5);
#endif
  }



  /******   depth buffer pass *****/
  
  glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
  glDepthFunc(GL_LEQUAL);
  glDisable(GL_BLEND);
  glEnable(GL_DEPTH_TEST);
  glDepthMask(GL_TRUE);  // don't write depth

  prog->enable();
  glBindVertexArray(mSizeVao);

  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  // VS
  prog->setUniform1f("spriteScale", viewport[3] / mInvFocalLen);
  prog->setUniform1f("starScale", powf(10.0f, m_starScaleLog));
  prog->setUniform1f("starAlpha", m_starAlpha);
  prog->setUniform1f("dmScale",  powf(10.0f, m_dmScaleLog));
  prog->setUniform1f("dmAlpha",  m_dmAlpha);
  prog->setUniform1f("spriteSizeMax", powf(10.0f, m_spriteSizeMaxLog));
  // PS
  prog->bindTexture("spriteTex",  m_sphTex, GL_TEXTURE_2D, 1);
  prog->setUniform1f("alphaScale", m_spriteAlpha);
  prog->setUniform1f("transmission", m_transmission);
  prog->setUniform1f("resx", m_imageW);
  prog->setUniform1f("resy", m_imageH);

  prog->setUniformfv("p0o", (GLfloat*)&m_clippingPlane[0], 4, 1);
  prog->setUniformfv("p1o", (GLfloat*)&m_clippingPlane[1], 4, 1);
  prog->setUniformfv("p2o", (GLfloat*)&m_clippingPlane[2], 4, 1);
  prog->setUniformfv("p3o", (GLfloat*)&m_clippingPlane[3], 4, 1);
  prog->setUniformfv("p4o", (GLfloat*)&m_clippingPlane[4], 4, 1);
  prog->setUniformfv("p5o", (GLfloat*)&m_clippingPlane[5], 4, 1);

  prog->setUniform1f("sorted", 2.0);

  //glClientActiveTexture(GL_TEXTURE0);
  glActiveTexture(GL_TEXTURE0);
  glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);
  glEnable(GL_POINT_SPRITE_ARB);

  drawPoints(start,count,true);

  prog->disable();

  /***** generate image pass *****/
 
  glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
  glDisable(GL_DEPTH_TEST);
  glDepthMask(GL_FALSE);  // don't write depth
  glEnable(GL_BLEND);
  m_fbo->AttachTexture(GL_TEXTURE_2D, 0, GL_DEPTH_ATTACHMENT_EXT);

  prog->enable();

  // VS
  prog->setUniform1f("spriteScale", viewport[3] / mInvFocalLen);
  prog->setUniform1f("starScale", powf(10.0f, m_starScaleLog));
  prog->setUniform1f("starAlpha", m_starAlpha);
  prog->setUniform1f("dmScale",  powf(10.0f, m_dmScaleLog));
  prog->setUniform1f("dmAlpha",  m_dmAlpha);
  prog->setUniform1f("spriteSizeMax", powf(10.0f, m_spriteSizeMaxLog));
  // PS
  prog->bindTexture("spriteTex",  m_sphTex, GL_TEXTURE_2D, 1);
  prog->setUniform1f("alphaScale", m_spriteAlpha);
  prog->setUniform1f("transmission", m_transmission);
  prog->setUniform1f("resx", m_imageW);
  prog->setUniform1f("resy", m_imageH);

  prog->setUniform1f("sorted", 1.0);

  //glClientActiveTexture(GL_TEXTURE0);
  glActiveTexture(GL_TEXTURE0);
  glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);
  glEnable(GL_POINT_SPRITE_ARB);

  drawPoints(start,count,true);

  prog->disable();


#ifdef _SPLOTCHSPRITES 
  glDisable(GL_CLIP_PLANE0);
  glDisable(GL_CLIP_PLANE1);
  glDisable(GL_CLIP_PLANE2);
  glDisable(GL_CLIP_PLANE3);
  glDisable(GL_CLIP_PLANE4);
  glDisable(GL_CLIP_PLANE5);
#else
  glDisable(GL_CLIP_DISTANCE0);
  glDisable(GL_CLIP_DISTANCE1);
  glDisable(GL_CLIP_DISTANCE2);
  glDisable(GL_CLIP_DISTANCE3);
  glDisable(GL_CLIP_DISTANCE4);
  glDisable(GL_CLIP_DISTANCE5);
#endif
  /********* compose ********/

#if 1
  glFlush();
  glFinish();
#endif


#if 1
  {
    GLint w, h, internalformat;

    glBindTexture(GL_TEXTURE_2D, m_imageTex[0]);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH,  &w);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_INTERNAL_FORMAT, &internalformat);
    glBindTexture(GL_TEXTURE_2D,0);

    static std::vector<float4> imgLoc, imgGlb;
    static std::vector<float > depth, depthLoc;
    imgLoc.resize(2*w*h);
    imgGlb.resize(2*w*h);
    depth.resize(2*w*h);


    static GLuint pbo_id[2];
    if (!pbo_id[0])
    {
      const int pbo_size = 8*4096*3072*sizeof(float4);
      glGenBuffers(2, pbo_id);
      glBindBuffer(GL_PIXEL_PACK_BUFFER, pbo_id[0]);
      glBufferData(GL_PIXEL_PACK_BUFFER, pbo_size, 0, GL_STATIC_READ);
      glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
      glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pbo_id[1]);
      glBufferData(GL_PIXEL_UNPACK_BUFFER, pbo_size, 0, GL_STATIC_DRAW);
      glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
    }
    assert(pbo_id[0] && pbo_id[1]);
   
    /* determine visible viewport bounds */ 

    const auto &viewport = getVisibleViewport();
    int2 wCrd  = make_int2(viewport[0], viewport[1]);
    int2 wSize = make_int2(viewport[2], viewport[3]);
    const int2 viewPort = make_int2(mWindowW,mWindowH);

    /***** fetch depth buffer *****/
    glFinish();
    MPI_Barrier(comm);
    const double t00 = MPI_Wtime();

    glBindTexture(GL_TEXTURE_2D, m_depthTex);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, pbo_id[0]);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
    GLvoid *rptr = glMapBufferRange(GL_PIXEL_PACK_BUFFER, 0, w*h*sizeof(float), GL_MAP_READ_BIT);
    
    glFinish();
    MPI_Barrier(comm);
    const double t10 = MPI_Wtime();

    float dmin = +HUGE, dmax = -HUGE;
#pragma omp parallel for schedule(static) reduction(min:dmin) reduction(max:dmax)
    for (int i = 0; i < w*h; i++)
    {
      depth[i] = reinterpret_cast<float*>(rptr)[i];
      dmin = std::min(dmin, depth[i]);
      dmax = std::max(dmax, depth[i]);
    }
//    fprintf(stderr, "rank= %d: dmin= %g  dmax= %g\n", rank, dmin,dmax);

    glUnmapBuffer(GL_PIXEL_PACK_BUFFER);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    glBindTexture(GL_TEXTURE_2D,0);
    
    glFinish();
    MPI_Barrier(comm);
    const double t20 = MPI_Wtime();
   
    /* determine real bounds to which pixels are written */ 
    { 
      int xmin = w;
      int ymin = h;
      int xmax = 0;
      int ymax = 0;

#pragma omp parallel for schedule(static) collapse(2) reduction(min:xmin,ymin) reduction(max:xmax,ymax)
      for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++)
          if (depth[j*w+i] < 1.0f)
          {
            xmin = std::min(xmin,i);
            ymin = std::min(ymin,j);
            xmax = std::max(xmax,i+1);
            ymax = std::max(ymax,j+1);
          }
      
      if (xmin == w)
        xmin = ymin = xmax = ymax = 0.0;

      wCrd  = make_int2(xmin,ymin);
      wSize = make_int2(xmax-xmin, ymax-ymin);
    }

    glFinish();
    MPI_Barrier(comm);
    const double t30 = MPI_Wtime();

    /***** fetch image *****/

    glBindTexture(GL_TEXTURE_2D, m_imageTex[0]);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, pbo_id[0]);
    if (wSize.x*wSize.y > 0)
      glReadPixels(wCrd.x, wCrd.y, wSize.x, wSize.y, GL_RGBA, GL_FLOAT, 0);

    glFinish();
    MPI_Barrier(comm);
    const double t40 = MPI_Wtime();

    if (wSize.x*wSize.y > 0)
      rptr = glMapBufferRange(GL_PIXEL_PACK_BUFFER, 0, wSize.x*wSize.y*sizeof(float4), GL_MAP_READ_BIT);

    if (wSize.x*wSize.y > 0)
    {
      imgLoc.resize(wSize.x*wSize.y);
      depthLoc.resize(wSize.x*wSize.y);
    }
    else
    {
      imgLoc.clear();
      depthLoc.clear();
    }
    
    glFinish();
    MPI_Barrier(comm);
    const double t50 = MPI_Wtime();

#pragma omp parallel for schedule(static)
    for (int i = 0; i < wSize.x*wSize.y; i++)
    {
      imgLoc[i] = reinterpret_cast<float4*>(rptr)[i];
      const int ix = i % wSize.x;
      const int iy = i / wSize.x;
      depthLoc[i] = depth[(iy+wCrd.y)*w+(ix+wCrd.x)];
    }

    glUnmapBuffer(GL_PIXEL_PACK_BUFFER);
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    glBindTexture(GL_TEXTURE_2D,0);
    
    glFinish();
    MPI_Barrier(comm);
    const double t60 = MPI_Wtime();

    compositingOrder.clear();
    lCompose(
        &imgLoc[0], &depthLoc[0], &imgGlb[0], 
        rank, nrank, comm,
        wCrd, wSize, viewPort,
        m_domainView ? m_domainViewIdx : -1,
        compositingOrder);

    glFinish();
    MPI_Barrier(comm);
    const double t70 = MPI_Wtime();

    if (isMaster())
    {
      /***** place back to fbo *****/

      glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pbo_id[1]);
      GLvoid *wptr = glMapBufferRange(GL_PIXEL_UNPACK_BUFFER, 0, w*h*sizeof(float4), GL_MAP_WRITE_BIT);

#pragma omp parallel for schedule(static)
      for (int i = 0; i < w*h; i++)
        reinterpret_cast<float4*>(wptr)[i] = imgGlb[i];

      glFinish();
      const double t80 = MPI_Wtime();

      glUnmapBuffer(GL_PIXEL_UNPACK_BUFFER);

      glBindTexture(GL_TEXTURE_2D, m_imageTex[0]);
      glTexImage2D(GL_TEXTURE_2D, 0, internalformat, w,h,0,GL_RGBA,GL_FLOAT, 0);
      glFinish();
      glBindTexture(GL_TEXTURE_2D,0);
      glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);

      glFinish();
      const double t90 = MPI_Wtime();

#if 0
      if (1)
        fprintf(stderr, 
            "total= %g: d2h= %g cpy= %g  mpi= %g  cpy= %g h2d= %g :: bwMPI= %g bwD2H= %g  bwH2D= %g\n", t6-t0,
                          t2-t1, t3-t2,   t4-t3,   t5-t4,  t6-t5,
            3.0*imgSize/(t4-t3)/1e6, imgSize/(t2-t1)/1e6, imgSize/(t6-t5)/1e6);
#else
      if (1)
        fprintf(stderr, "total= %g: depth= [%g %g %g] img= [%g %g %g] mpi= %g  wb= [ %g %g ]\n", t90 - t00,
            t10-t00, t20-t10, t30-t20,
            t40-t30, t50-t40, t60-t50,
            t70-t60,
            t80-t70, t90-t80);
#endif
    }

  }
#endif

  m_fbo->Disable();



  glDisable(GL_BLEND);
  m_splotch2texProg->enable();
  m_splotch2texProg->bindTexture("tex", m_imageTex[0], GL_TEXTURE_2D, 0);
  m_splotch2texProg->setUniform1f("scale_pre", 0.1*m_imageBrightnessPre);
  m_splotch2texProg->setUniform1f("gamma_pre", m_gammaPre);
  m_splotch2texProg->setUniform1f("scale_post", m_imageBrightnessPost);
  m_splotch2texProg->setUniform1f("gamma_post", m_gammaPost);
  m_splotch2texProg->setUniform1f("sorted", 1.0f);
  drawQuad();
  m_splotch2texProg->disable();
}

// render scene depth to texture
// (this is to ensure that particles are correctly occluded in the low-resolution render buffer)
void SmokeRenderer::beginSceneRender(Target target)
{
  m_fbo->Bind();
  if (target == LIGHT_BUFFER) {
    m_fbo->AttachTexture(GL_TEXTURE_2D, m_lightTexture[m_srcLightTexture], GL_COLOR_ATTACHMENT0_EXT);
    m_fbo->AttachTexture(GL_TEXTURE_2D, m_lightDepthTexture, GL_DEPTH_ATTACHMENT_EXT);

    glViewport(0, 0, m_lightBufferSize, m_lightBufferSize);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadMatrixf((GLfloat *) m_lightView.get_value());

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadMatrixf((GLfloat *) m_lightProj.get_value());
  } else {
    m_fbo->AttachTexture(GL_TEXTURE_2D, m_imageTex[0], GL_COLOR_ATTACHMENT0_EXT);
    m_fbo->AttachTexture(GL_TEXTURE_2D, m_depthTex, GL_DEPTH_ATTACHMENT_EXT);

    glViewport(0, 0, m_imageW, m_imageH);
  }
  glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
  glDepthMask(GL_TRUE);
  glClear(GL_DEPTH_BUFFER_BIT);
}

void SmokeRenderer::endSceneRender(Target target)
{
  m_fbo->Disable();
  if (target == LIGHT_BUFFER) {
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
  }
  glViewport(0, 0, mWindowW, mWindowH);
  glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
}

// create an OpenGL texture
  GLuint
SmokeRenderer::createTexture(GLenum target, int w, int h, GLint internalformat, GLenum format, void *data)
{
  GLuint texid;
  glGenTextures(1, &texid);
  glBindTexture(target, texid);

  glTexParameteri(target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(target, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(target, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  glTexImage2D(target, 0, internalformat, w, h, 0, format, GL_FLOAT, data);
  return texid;
}

inline float frand()
{
  return rand() / (float) RAND_MAX;
}

inline float sfrand()
{
  return frand()*2.0f-1.0f;
}

GLuint SmokeRenderer::createNoiseTexture(int w, int h, int d)
{
  int size = w*h*d;
  float *data = new float [size];
  float *ptr = data;
  for(int i=0; i<size; i++) {
    *ptr++ = sfrand();
    //*ptr++ = sfrand();
    //*ptr++ = sfrand();
    //*ptr++ = sfrand();
  }

  GLuint texid;
  glGenTextures(1, &texid);
  GLenum target = GL_TEXTURE_3D;
  glBindTexture(target, texid);

  glTexParameteri(target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(target, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(target, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(target, GL_TEXTURE_WRAP_T, GL_REPEAT);

  glTexImage3D(GL_TEXTURE_3D, 0, GL_LUMINANCE16F_ARB, w, h, d, 0, GL_LUMINANCE, GL_FLOAT, data);
  delete [] data;

  return texid;
}

float * SmokeRenderer::createSplatImage(int n)
{
  float *img = new float[n*n];
  for(int y=0; y<n; y++) {
    float v = (y / (float) (n-1))*2.0f-1.0f;
    for(int x=0; x<n; x++) {
      float u = (x / (float) (n-1))*2.0f-1.0f;
      float d = sqrtf(u*u + v*v);
      if (d > 1.0f) d = 1.0f;
      float i = 1.0f - d*d*(3.0f - 2.0f*d);	// smoothstep
      img[y*n+x] = i;
    }
  }
  return img;
}

GLuint SmokeRenderer::createSpriteTexture(int size)
{
  float *img = createSplatImage(size);

  GLuint tex = createTexture(GL_TEXTURE_2D, size, size, GL_LUMINANCE8, GL_LUMINANCE, img);
  delete [] img;
  glGenerateMipmapEXT(GL_TEXTURE_2D);
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  return tex;
}
static inline float lWkernel(const float q2)
{
  const float q = sqrtf(q2);
  const float sigma = 8.0f/M_PI;

  const float qm = 1.0f - q;
  if      (q < 0.5f) return sigma * (1.0f + (-6.0f)*q*q*qm);
  else if (q < 1.0f) return sigma * 2.0f*qm*qm*qm;

  return 0.0f;
}
GLuint SmokeRenderer::createSphTexture(int size)
{
  const float scale = 1.0f/lWkernel(0.0f);
  float *img = new float[size*size];
  for (int j = 0; j < size; j++)
    for (int i = 0; i < size; i++)
    {
      const float dx = ((i+0.5f)/size - 0.5f) * 2.01f;
      const float dy = ((j+0.5f)/size - 0.5f) * 2.01f;
      const float q2 = dx*dx + dy*dy;
      img[j*size+i] = lWkernel(q2)*scale;
    }

  GLuint tex = createTexture(GL_TEXTURE_2D, size, size, GL_LUMINANCE8, GL_LUMINANCE, img);
  delete [] img;
  glGenerateMipmapEXT(GL_TEXTURE_2D);
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  return tex;
}

// create textures for off-screen rendering
void SmokeRenderer::createBuffers(int w, int h)
{
  if (m_imageTex[0]) {
    glDeleteTextures(4, m_imageTex);
    glDeleteTextures(1, &m_depthTex);

    glDeleteTextures(3, m_downSampledTex);
  }

  mWindowW = w;
  mWindowH = h;

  m_imageW = w / m_downSample;
  m_imageH = h / m_downSample;
  if (isMaster())
    printf("image size: %d %d\n", m_imageW, m_imageH);

  // create texture for image buffer
  GLint format = GL_RGBA32F;
// format = GL_RGBA16F;
  //GLint format = GL_LUMINANCE16F_ARB;
  //GLint format = GL_RGBA8;
  m_imageTex[0] = createTexture(GL_TEXTURE_2D, m_imageW, m_imageH, format, GL_RGBA);
  m_imageTex[1] = createTexture(GL_TEXTURE_2D, m_imageW, m_imageH, format, GL_RGBA);
  m_imageTex[2] = createTexture(GL_TEXTURE_2D, m_imageW, m_imageH, format, GL_RGBA);
  m_imageTex[3] = createTexture(GL_TEXTURE_2D, m_imageW, m_imageH, format, GL_RGBA);
  m_imageTex[4] = createTexture(GL_TEXTURE_2D, m_imageW, m_imageH, format, GL_RGBA);

//  m_depthTex = createTexture(GL_TEXTURE_2D, m_imageW, m_imageH, GL_DEPTH_COMPONENT24_ARB, GL_DEPTH_COMPONENT);
  m_depthTex = createTexture(GL_TEXTURE_2D, m_imageW, m_imageH, GL_DEPTH_COMPONENT32_ARB, GL_DEPTH_COMPONENT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

  m_downSampledW = m_imageW / m_blurDownSample;
  m_downSampledH = m_imageH / m_blurDownSample;
  if (isMaster())
    printf("downsampled size: %d %d\n", m_downSampledW, m_downSampledH);

  m_downSampledTex[0] = createTexture(GL_TEXTURE_2D, m_downSampledW, m_downSampledH, format, GL_RGBA);
  m_downSampledTex[1] = createTexture(GL_TEXTURE_2D, m_downSampledW, m_downSampledH, format, GL_RGBA);
  m_downSampledTex[2] = createTexture(GL_TEXTURE_2D, m_downSampledW, m_downSampledH, format, GL_RGBA);

  createLightBuffer();
}

// create textures for light buffer
  void
SmokeRenderer::createLightBuffer()
{
  if (m_lightTexture[0]) {
    glDeleteTextures(1, &m_lightTexture[0]);
    glDeleteTextures(1, &m_lightTexture[1]);
    glDeleteTextures(1, &m_lightDepthTexture);
  }

  GLint format = GL_RGBA16F_ARB;
  //GLint format = GL_RGBA8;
  //GLint format = GL_LUMINANCE16F_ARB;

#if USE_MRT
  // textures must be same size to be bound to same FBO at same time
  m_lightTextureW = std::max(m_lightBufferSize, m_imageW);
  m_lightTextureH = std::max(m_lightBufferSize, m_imageH);
#else
  m_lightTextureW = m_lightBufferSize;
  m_lightTextureH = m_lightBufferSize;
#endif

  m_lightTexture[0] = createTexture(GL_TEXTURE_2D, m_lightTextureW, m_lightTextureH, format, GL_RGBA);
  // make shadows clamp to light color at edges
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

  m_lightTexture[1] = createTexture(GL_TEXTURE_2D, m_lightTextureW, m_lightTextureH, format, GL_RGBA);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

  m_lightDepthTexture = createTexture(GL_TEXTURE_2D, m_lightTextureW, m_lightTextureH, GL_DEPTH_COMPONENT24_ARB, GL_DEPTH_COMPONENT);

  m_fbo->AttachTexture(GL_TEXTURE_2D, m_lightTexture[m_srcLightTexture], GL_COLOR_ATTACHMENT0_EXT);
  m_fbo->AttachTexture(GL_TEXTURE_2D, 0, GL_COLOR_ATTACHMENT1_EXT);
  m_fbo->AttachTexture(GL_TEXTURE_2D, m_lightDepthTexture, GL_DEPTH_ATTACHMENT_EXT);
  m_fbo->IsValid();
}

  void
SmokeRenderer::setLightColor(vec3f c)
{
  m_lightColor = c;

  // set light texture border color
  //    GLfloat borderColor[4] = { 1.0 - m_lightColor[0], 1.0 - m_lightColor[1], 1.0 - m_lightColor[2], 0.0 };
  GLfloat borderColor[4] = { m_lightColor[0], m_lightColor[1], m_lightColor[2], 0.0 };

  glBindTexture(GL_TEXTURE_2D, m_lightTexture[0]);
  glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

  glBindTexture(GL_TEXTURE_2D, m_lightTexture[1]);
  glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

  glBindTexture(GL_TEXTURE_2D, 0);
}

void SmokeRenderer::setWindowSize(int w, int h)
{
  mAspect = (float) mWindowW / (float) mWindowH;
  mInvFocalLen = (float) tan(mFov*0.5*NV_PI/180.0);

  createBuffers(w, h);
}

void SmokeRenderer::drawQuad(float s, float z)
{
  glBegin(GL_QUADS);
  glTexCoord2f(0.0, 0.0); glVertex3f(-s, -s, z);
  glTexCoord2f(1.0, 0.0); glVertex3f(s, -s, z);
  glTexCoord2f(1.0, 1.0); glVertex3f(s, s, z);
  glTexCoord2f(0.0, 1.0); glVertex3f(-s, s, z);
  glEnd();
}

void SmokeRenderer::drawVector(vec3f v)
{
  glBegin(GL_LINES);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3fv((float *) &v[0]);
  glEnd();
}

// render vectors to screen for debugging
void SmokeRenderer::debugVectors()
{
  glColor3f(1.0, 1.0, 0.0);
  drawVector(m_lightVector);

  glColor3f(0.0, 1.0, 0.0);
  drawVector(m_viewVector);

  glColor3f(0.0, 0.0, 1.0);
  drawVector(-m_viewVector);

  glColor3f(1.0, 0.0, 0.0);
  drawVector(m_halfVector);
}

void SmokeRenderer::drawSkybox(GLuint tex)
{
#if 0
  if (!m_cubemapTex)
    return;
#else
  return;
#endif

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_CUBE_MAP, tex);
  m_skyboxProg->enable();
  m_skyboxProg->bindTexture("tex", tex, GL_TEXTURE_CUBE_MAP, 0);

  glDisable(GL_DEPTH_TEST);
  glDepthMask(GL_FALSE);
  glColor3f(m_skyboxBrightness, m_skyboxBrightness, m_skyboxBrightness);
  //glColor3f(0.25f, 0.25f, 0.25f);
  //glColor3f(0.5f, 0.5f, 0.5f);
  //glColor3f(1.0f, 1.0f, 1.0f);

  glutSolidCube(2.0);

  m_skyboxProg->disable();

  glEnable(GL_DEPTH_TEST);
  glDepthMask(GL_TRUE);
}

void SmokeRenderer::initParams()
{
  // m_splotchProg
  // ---------------
  // spriteScale
  // starScale
  // starAlpha
  // dmScale
  // dmAlpha
  // spriteSizeMax
  // alphaScale
  // transmission
  // gamma pre/post
  // brightness pre/post
  //
  //  m_compositeProg
  //  ---------------------
  //  imageTex
  //  blurTexH
  //  blurTexV
  //  glowTex
  //  flareTex
  //  scale 
  //  sourceIntensity
  //  glowIntensity
  //  starInensity
  //  flareIntensity
  //  gamma
  //
  //  m_thresholdProg
  //  ----------
  //  scale (starPower)
  //  threshold (m_starThreshold)
  //
  //  m_starFilterProg
  //  ----------
  //  radius
  //  texeslSize
  //
  //  m_gaussiaBlurProg
  //  -----------
  //  radius
  m_params = new ParamListGL("render_params");

  m_params->AddParam(new Param<float>("star scale [log]", m_starScaleLog,     -1.0f, 1.0f, 0.001f, &m_starScaleLog));
  m_params->AddParam(new Param<float>("star alpha      ", m_starAlpha,         0.0f, 1.0f, 0.001f, &m_starAlpha));

  m_params->AddParam(new Param<float>("dm scale   [log]", m_dmScaleLog,       -1.0f, 1.0f, 0.001f, &m_dmScaleLog));
  m_params->AddParam(new Param<float>("dm alpha        ", m_dmAlpha,           0.0f, 1.0f, 0.001f, &m_dmAlpha));

  m_params->AddParam(new Param<float>("max size [log]  ", m_spriteSizeMaxLog, -1.0f, 1.0f, 0.001f, &m_spriteSizeMaxLog));
  m_params->AddParam(new Param<float>("alpha           ", m_spriteAlpha,       0.0f, 1.0f, 0.001f, &m_spriteAlpha));
  m_params->AddParam(new Param<float>("transmission    ", m_transmission,      0.0f, 0.1f, 0.001f, &m_transmission));

  m_params->AddParam(new Param<float>("brightness [pre]",  m_imageBrightnessPre, 0.0f, 1.0f, 0.001f, &m_imageBrightnessPre));
  m_params->AddParam(new Param<float>("gamma [pre]",       m_gammaPre,           0.0f, 2.0f, 0.001f, &m_gammaPre));
  m_params->AddParam(new Param<float>("brightness [post]", m_imageBrightnessPost, 0.0f, 1.0f, 0.001f, &m_imageBrightnessPost));
  m_params->AddParam(new Param<float>("gamma [post]",      m_gammaPost,           0.0f, 2.0f, 0.001f, &m_gammaPost));
 
#if 0 
  m_params->AddParam(new Param<int>("slices", m_numSlices, 1, 256, 1, &m_numSlices));
  m_params->AddParam(new Param<int>("displayed slices", m_numDisplayedSlices, 1, 256, 1, &m_numDisplayedSlices));
  

  m_params->AddParam(new Param<float>("sprite size", mParticleRadius, 0.0f, 0.2f, 0.001f, &mParticleRadius));
  m_params->AddParam(new Param<float>("scale [log]", mParticleScaleLog, -1.0f, 1.0f, 0.01f, &mParticleScaleLog));
   
  m_params->AddParam(new Param<float>("dust scale", m_ageScale, 0.0f, 50.0f, 0.1f, &m_ageScale));
  m_params->AddParam(new Param<float>("dust alpha", m_dustAlpha, 0.0f, 0.1f, 0.01f, &m_dustAlpha));

  m_params->AddParam(new Param<float>("light color r", m_lightColor[0], 0.0f, 1.0f, 0.01f, &m_lightColor[0]));
  m_params->AddParam(new Param<float>("light color g", m_lightColor[1], 0.0f, 1.0f, 0.01f, &m_lightColor[1]));
  m_params->AddParam(new Param<float>("light color b", m_lightColor[2], 0.0f, 1.0f, 0.01f, &m_lightColor[2]));

#if 0
  m_params->AddParam(new Param<float>("color opacity r", m_colorOpacity[0], 0.0f, 1.0f, 0.01f, &m_colorOpacity[0]));
  m_params->AddParam(new Param<float>("color opacity g", m_colorOpacity[1], 0.0f, 1.0f, 0.01f, &m_colorOpacity[1]));
  m_params->AddParam(new Param<float>("color opacity b", m_colorOpacity[2], 0.0f, 1.0f, 0.01f, &m_colorOpacity[2]));
#endif

  m_params->AddParam(new Param<float>("alpha", m_spriteAlpha, 0.0f, 1.0f, 0.001f, &m_spriteAlpha));
  m_params->AddParam(new Param<float>("shadow alpha", m_shadowAlpha, 0.0f, 1.0f, 0.001f, &m_shadowAlpha));
  m_params->AddParam(new Param<float>("transmission", m_transmission, 0.0f, 0.1f, 0.001f, &m_transmission));
  m_params->AddParam(new Param<float>("indirect lighting", m_indirectAmount, 0.0f, 1.0f, 0.001f, &m_indirectAmount));

#if 0
  // volume stuff
  m_params->AddParam(new Param<float>("volume alpha", m_volumeAlpha, 0.0f, 1.0f, 0.01f, &m_volumeAlpha));
  m_params->AddParam(new Param<float>("volume indirect", m_volumeIndirect, 0.0f, 1.0f, 0.01f, &m_volumeIndirect));

  m_params->AddParam(new Param<float>("volume color r", m_volumeColor[0], 0.0f, 1.0f, 0.01f, &m_volumeColor[0]));
  m_params->AddParam(new Param<float>("volume color g", m_volumeColor[1], 0.0f, 1.0f, 0.01f, &m_volumeColor[1]));
  m_params->AddParam(new Param<float>("volume color b", m_volumeColor[2], 0.0f, 1.0f, 0.01f, &m_volumeColor[2]));
  m_params->AddParam(new Param<float>("volume noise freq", m_noiseFreq, 0.0f, 1.0f, 0.01f, &m_noiseFreq));
  m_params->AddParam(new Param<float>("volume noise amp", m_noiseAmp, 0.0f, 2.0f, 0.01f, &m_noiseAmp));
  m_params->AddParam(new Param<float>("volume start", m_volumeStart, 0.0f, 1.0f, 0.01f, &m_volumeStart));
  m_params->AddParam(new Param<float>("volume width", m_volumeWidth, 0.0f, 1.0f, 0.01f, &m_volumeWidth));
#endif

  m_params->AddParam(new Param<float>("fog", m_fog, 0.0f, 0.1f, 0.001f, &m_fog));

  m_params->AddParam(new Param<float>("over bright multiplier", m_overBright, 0.0f, 100.0f, 1.0f, &m_overBright));
  //m_params->AddParam(new Param<float>("over bright threshold", m_overBrightThreshold, 0.0f, 1.0f, 0.001f, &m_overBrightThreshold));
  m_params->AddParam(new Param<float>("star brightness", m_overBrightThreshold, 0.0f, 10.0f, 0.001f, &m_overBrightThreshold));
  m_params->AddParam(new Param<float>("image brightness", m_imageBrightness, 0.0f, 2.0f, 0.1f, &m_imageBrightness));
  m_params->AddParam(new Param<float>("image gamma", m_gamma, 0.0f, 2.0f, 0.0f, &m_gamma));

  m_params->AddParam(new Param<float>("blur radius", m_blurRadius, 0.0f, 10.0f, 0.1f, &m_blurRadius));
  m_params->AddParam(new Param<int>("blur passes", m_blurPasses, 0, 10, 1, &m_blurPasses));

  m_params->AddParam(new Param<float>("source intensity", m_sourceIntensity, 0.0f, 1.0f, 0.01f, &m_sourceIntensity));
  m_params->AddParam(new Param<float>("star blur radius", m_starBlurRadius, 0.0f, 100.0f, 1.0f, &m_starBlurRadius));
  m_params->AddParam(new Param<float>("star threshold", m_starThreshold, 0.0f, 100.0f, 0.1f, &m_starThreshold));
  m_params->AddParam(new Param<float>("star power", m_starPower, 0.0f, 100.0f, 0.1f, &m_starPower));
  m_params->AddParam(new Param<float>("star intensity", m_starIntensity, 0.0f, 1.0f, 0.1f, &m_starIntensity));
  m_params->AddParam(new Param<float>("glow radius", m_glowRadius, 0.0f, 100.0f, 1.0f, &m_glowRadius));
  m_params->AddParam(new Param<float>("glow intensity", m_glowIntensity, 0.0f, 1.0f, 0.01f, &m_glowIntensity));
  m_params->AddParam(new Param<float>("flare intensity", m_flareIntensity, 0.0f, 1.0f, 0.01f, &m_flareIntensity));
  m_params->AddParam(new Param<float>("flare threshold", m_flareThreshold, 0.0f, 10.0f, 0.01f, &m_flareThreshold));
  m_params->AddParam(new Param<float>("flare radius", m_flareRadius, 0.0f, 100.0f, 0.01f, &m_flareRadius));

  m_params->AddParam(new Param<float>("skybox brightness", m_skyboxBrightness, 0.0f, 1.0f, 0.01f, &m_skyboxBrightness));
#endif
}

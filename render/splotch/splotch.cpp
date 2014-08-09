#include "Splotch.h"

float4 Splotch::modelView(const float4 pos) const
{
  auto &m = modelViewMatrix;
  return make_float4(
      m[0][0]*pos.x + m[0][1]*pos.y + m[0][2]*pos.z + m[0][3]*pos.w,
      m[1][0]*pos.x + m[1][1]*pos.y + m[1][2]*pos.z + m[1][3]*pos.w,
      m[2][0]*pos.x + m[2][1]*pos.y + m[2][2]*pos.z + m[2][3]*pos.w,
      m[3][0]*pos.x + m[3][1]*pos.y + m[3][2]*pos.z + m[3][3]*pos.w);
}
float4 Splotch::projection(const float4 pos) const
{
  auto &m = projectionMatrix;
  return make_float4(
      m[0][0]*pos.x + m[0][1]*pos.y + m[0][2]*pos.z + m[0][3]*pos.w,
      m[1][0]*pos.x + m[1][1]*pos.y + m[1][2]*pos.z + m[1][3]*pos.w,
      m[2][0]*pos.x + m[2][1]*pos.y + m[2][2]*pos.z + m[2][3]*pos.w,
      m[3][0]*pos.x + m[3][1]*pos.y + m[3][2]*pos.z + m[3][3]*pos.w);
}

void Splotch::transform(const bool perspective)
{
  const int np = vtxArray.size();
  vtxArrayView.realloc(np);
  depthArray.resize(np);

  const auto &colorMapTex = *colorMapTexPtr;


  int nVisible = 0;
#pragma omp parallel for schedule(runtime) reduction(+:nVisible)
  for (int i = 0; i < np; i++)
  {
    const auto &vtx = vtxArray[i];
    const float4 posO = make_float4(vtx.pos.x,vtx.pos.y,vtx.pos.z,1.0f);
    float4 posV = modelView(posO);
    const float depth = posV.z;
    const float dist  = length(posV);
    float3 col = make_float3(-1.0f);


    if (depth >= depthMin && depth <= depthMax)
    {
      posV = projection(posV);
      posV.x *= 1.0f/posV.z;
      posV.y *= 1.0f/posV.z;

      posV.x = (posV.x + 1.0f) * 0.5f * width;
      posV.y = (1.0f - posV.y) * 0.5f * height;
        
      posV.w = vtx.pos.h * 0.5f * width / dist;
      assert(vtx.pos.h > 0.0f);
      using std::sqrt;
      using std::max;
      posV.w = sqrt(posV.w*posV.w + minHpix*minHpix);
      posV.w = min(posV.w, maxHpix);
      assert(posV.w > 0.0f);

      if (   posV.x - posV.w <= width
          && posV.x + posV.w >= 0
          && posV.y - posV.w <= height
          && posV.y + posV.w >= 0)
      {
        const float s = vtx.attr.rho;
        const float t = vtx.attr.vel;
        assert(s>=0.0f && s<=1.0f);
        assert(t>=0.0f && t<=1.0f);
        const auto &tex = colorMapTex(s,t);
        col = make_float3(tex[0],tex[1],tex[2]);
      }
      else
        posV.w = -1.0;
    }
    else
      posV.w = -1.0f;

    depthArray  [i] = depth;
    vtxArrayView[i] = 
    {
      pos2d_t(posV.x, posV.y, posV.w),
      make_float4(col, 1.0f),
      vtx.attr
    };

    nVisible += vtxArrayView[i].isVisible();
  }
  fprintf(stderr, "nParticles= %d nVisible= %d\n", np, nVisible);
}

void Splotch::depthSort()
{
  const int np = depthArray.size();

  using pair = std::pair<float,int>;
  std::vector<pair> depthMap;
  depthMap.reserve(np);

  for (int i = 0; i < np; i++)
    if (vtxArrayView[i].isVisible())
      depthMap.push_back(std::make_pair(depthArray[i],i));

  __gnu_parallel::sort(depthMap.begin(), depthMap.end(),
      [](const pair &a, const pair &b) { return a.first < b.first;} );

  const int npVis = depthMap.size();
  VertexArrayView vtxView(npVis);

#pragma omp parallel for 
  for (int i = 0; i < npVis; i++)
  {
    const auto &map = depthMap[i];
    assert(map.second >= 0);
    assert(map.second < np);
#if 0
    fprintf(stderr, "i= %d  npVis= %d  map.second= %d  np= %d\n",
        i, npVis, map.second, np);
#endif
    vtxView[i] = vtxArrayView[map.second];
  }

  swap(vtxArrayView,vtxView);
}

// assumes atomic execution
Splotch::Quad Splotch::rasterize(const VertexView &vtx, const Splotch::Quad &range, std::vector<color_t> &fb)
{
  using std::max;
  using std::min;
  using std::exp;
  using std::floor;
  using std::ceil;

  Quad q;
  q.x0  = max(range.x0, floor(vtx.pos.x - vtx.pos.h));
  q.x1  = min(range.x1, ceil (vtx.pos.x + vtx.pos.h));
  q.y0  = max(range.y0, floor(vtx.pos.y - vtx.pos.h));
  q.y1  = min(range.y1, ceil (vtx.pos.y + vtx.pos.h));

  const float invh  = 1.0f/vtx.pos.h;
  const float invh2 = invh*invh;
  const float width  = range.x1 - range.x0;
  const float height = range.y1 - range.y0;
  for (float iy = q.y0; iy < q.y1; iy++)
  {
    const float dy   = iy - vtx.pos.y;
    const float qy   = dy*dy * invh2;
    const float facy = exp(-qy);
    for(float ix = q.x0; ix < q.x1; ix++)
    {
      const float dx = ix - vtx.pos.x;
      const float qx = dx*dx * invh2;
      const float facx = exp(-qx);

      float4 color = vtx.color;
      color.w = facx*facy; /* alpha */

      const int idx = (ix - range.x0) + width*(iy - range.y0);
      assert(idx >= 0);
      assert(idx < width*height);
      fb[idx] = Blending::getColor<Blending::ONE,Blending::SRC_ALPHA>(fb[idx],color);
    }
  }
  return q;
}

void Splotch::render()
{
  const int np = vtxArrayView.size();

  image.resize(width*height);
  std::fill(image.begin(), image.end(), make_float4(0.0f));

#define NTHREADMAX 256
  std::vector<color_t> fbVec[NTHREADMAX];

#pragma omp parallel
  {
    const int nt = omp_get_num_threads();
    assert(nt <= NTHREADMAX);
    const int tid = omp_get_thread_num();
    auto &fb = fbVec[tid];
    fb.resize(width*height);
    std::fill(fb.begin(), fb.end(), make_float4(0.0f));

    Quad range;
    range.x0 = 0;
    range.x1 = width;
    range.y0 = 0;
    range.y1 = height;

    if (tid == 0)
      fprintf(stderr, "rasterize begin .. \n");

#pragma omp for schedule(runtime)
    for (int i = 0; i < np; i++)
    {
#if 0
      fprintf(stderr, "i:= %d  np= %d: x= %g  y= %g  h= %g\n",
          i, np, 
          vtxArrayView[i].pos.x, vtxArrayView[i].pos.y, 
          vtxArrayView[i].pos.h);
#endif
      rasterize(vtxArrayView[i], range, fb);
    }

    if (tid == 0)
      fprintf(stderr, "rasterize end .. \n");

#pragma omp barrier

#pragma omp for schedule(runtime) collapse(2)
    for (int j = 0; j < height; j++)
      for (int i = 0; i < width; i++)
      {
        const int idx = j*width + i;
        for (int k = 0; k < nt; k++)
        {
          image[idx] = Blending::getColor<Blending::ONE,Blending::SRC_ALPHA>(image[idx], fbVec[k][idx]);
        }
      }
  }
}


void Splotch::finalize()
{
#pragma omp for schedule(runtime) collapse(2)
  for (int j = 0; j < height; j++)
    for (int i = 0; i < width; i++)
    {
      const int idx = j*width + i;
      const float4 src = image[idx];
      float4 dst;
#if 0
      fprintf(stderr, " (%3d,%3d): %g %g %g \n",
          i,j, src.x,src.y,src.z);
#endif

      dst.x = 1.0f - exp(-src.x);
      dst.y = 1.0f - exp(-src.y);
      dst.z = 1.0f - exp(-src.z);
      dst.w = 1.0f;

      image[idx] = dst;
    }
}


void Splotch::genImage(const bool perspective)
{
  fprintf(stderr , " --- transform \n");
  transform(perspective);
  fprintf(stderr , " --- depthSort \n");
  depthSort();
  fprintf(stderr , " --- render \n");
  render();
  fprintf(stderr , " --- finalize \n");
  finalize();
}

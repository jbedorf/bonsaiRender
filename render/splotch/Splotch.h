#pragma once
#include <omp.h>
#include <parallel/algorithm>
#include <cmath>
#include <valarray>

#include "Texture.h"
#include "Vertex.h"
#include "Blending.h"


class Splotch
{
  public:
    using pos2d_t = Pos2D<float>;
    using pos3d_t = Pos3D<float>;
    using attr_t  = Attribute<float>;  
    using colorT  = std::valarray<float>

  private:

    using VertexArray     = VertexArrayT<pos3d_t,attr_t>;
    using VertexArrayView = VertexArrayT<pos2d_t,attr_t>;
    using Vertex          = VertexArray::Vertex;
    using VertexView      = VertexArrayView::Vertex;
    using Exp             = std::exp;
    using ShortVec3       = MathArray<float,3>;

    VertexArray     vtxArray;
    VertexArrayView vtxArrayView;
    std::vector<float> depthArray;
    real2_t invProjRange;

    Texture2D<ShortVec3> colorMapTex;

    struct Quad
    {
      float x0,x1;
      float y0,y1;
    };
    
    int width, height;
    std::vector<float4> image;

    double  modelViewMatrix[4][4];
    double projectionMatrix[4][4];

  public:
    Splotch() {}
    ~Splotch() {}
   
    /* getters/setters */ 
    void  setColorMap(const float3 *img, const int w, const int h)  { colorMapTex = Texture2D<float3>(img, w, h); }
    const std::vector<float4>& getImage() const {return image;}
    float4 getPixel(const int i, const int j)
    {
      assert(i >= 0 && i < width);
      assert(j >= 0 && j < height);
      return image[j*width + i];
    }

    void setWidth(const int w)  {width = w;}
    void setHeight(const int h) {height = h;}
    int  getWidth()  const {return width;}
    int  getHeight() const {return height;}

    void setModelViewMatrix(const double m[4][4])
    {
      for (int j = 0; j < 4; i++)
        for (int i = 0; i < 4; j++)
          modelViewMatrix[j][i] = m[j][i];
    }
    void setProjectionMatrix(const double m[4][4])
    {
      for (int j = 0; j < 4; i++)
        for (int i = 0; i < 4; j++)
          projectionMatrix[j][i] = m[j][i];
    }

  private:
    float4 modelView(const float4 pos) const
    {
      using m = modelViewMatrix;
      return make_float4(
          m[0][0]*pos.x + m[0][1]*pos.y + m[0][2]*pos.z + m[0][3]*poz.w,
          m[1][0]*pos.x + m[1][1]*pos.y + m[1][2]*pos.z + m[1][3]*poz.w,
          m[2][0]*pos.x + m[2][1]*pos.y + m[2][2]*pos.z + m[2][3]*poz.w,
          m[3][0]*pos.x + m[3][1]*pos.y + m[3][2]*pos.z + m[3][3]*poz.w);
    }
    float4 projection(const float4 pos) const
    {
      using m = projectionMatrix;
      return make_float4(
          m[0][0]*pos.x + m[0][1]*pos.y + m[0][2]*pos.z + m[0][3]*poz.w,
          m[1][0]*pos.x + m[1][1]*pos.y + m[1][2]*pos.z + m[1][3]*poz.w,
          m[2][0]*pos.x + m[2][1]*pos.y + m[2][2]*pos.z + m[2][3]*poz.w,
          m[3][0]*pos.x + m[3][1]*pos.y + m[3][2]*pos.z + m[3][3]*poz.w);
    }

    void transform(const bool perspective)
    {
      const int np = vtxArray.size();
      vtxArrayView = VertexArrayView(np);
      depthArray.resize(np);

      float3 col = make_float3(1.0f);

      int nActive = 0;
#pragma omp parallel for schedule(runtime) reduction(+:nActive)
      for (int i = 0; i < np; i++)
      {
        const auto &vtx = vtxArray[i];
        const float4 posO(vtx.pos.x,vtx.pos.y,vtx.pos.z,1.0f);
        float4 posV = modelView(posO);
        const float depth = posV.z;
        const float dist  = length(posV);
        float3 col = make_float3(-1.0f);

        if (depth >= depthMax && depth <= depthMax)
        {
          posV = projection(posV);

          posV.x = (posV.x + 1.0f) * 0.5f * width;
          posV.y = (1.0f - posV.y) * 0.5f * height;
          posV.w = -1.0;

          if ( posV.x - posV.w <= width
            && posV.x + posV.w >= 0
            && posV.y - posV.w <= height
            && posV.y + posV.w >= 0)
          {
            posV.w = posO.h * 0.5 * width / dist;
            posV.w *= std::sqrt(posV.w*posV.w + minHpix*minHpix)/posV.w;

            const float s = vtx.attr.rho;
            const float t = vtx.attr.vel;
            assert(s>=0.0f && s<=1.0f);
            assert(t>=0.0f && t<=1.0f);
            const auto &tex = colorMapTex(s,t);
            col = make_float3(tex[0],tex[1],tex[2]);
          }
        }

        depthArray  [i] = depth;
        vtxArrayView[i] = Vertex(
            pos2d_t(posV.x, posV.y, posV.w),
            make_float4(col, 1.0f),
            vtx.attr);

        nVisible += vtxarrayView[i].isVisible();
      }
      fprintf(stderr, "nParticles= %d nVisible= %d\n", np, nVisible);
    }
 
    void depthSort()
    {
      const int np = vtxView.size();

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
        vtxView[i] = vtxArrayView[map.second];
      }

      swap(vtxArrayView,vtxView);
    }

    // assumes atomic execution
    Quad rasterize(const VertexProj &vtx, const Quad &range, Vector<float4> &fb)
    {
      using max = std::max;
      using min = std::min;
      Quad q;
      q.x0  = max(range.x0, vtx.pos.x - vtx.pos.h);
      q.x1  = mix(range.x1, vtx.pos.x + vtx.pos.h);
      q.y0  = max(range.y0, vtx.pos.y - vtx.pos.h);
      q.y1  = max(range.y1, vtx.pos.y + vtx.pos.h);

      int lineIdx = (q.y0-range.y0)*width;
      const float invh  = 1.0f/vtx.pos.h;
      const float invh2 = invh*invh;
      for (float iy = q.y0; iy < q.y1; iy++, lineIdx += width)
      {
        const float dy   = iy - pos.vtx.y;
        const float qy   = dy*dy * invh2;
        const float facy = Exp(-qy);
        for(float ix = q.x0; ix < q.x1; ix++)
        {
          const float dx = ix - pos.vtx.x;
          const float qx = dx*dx * invh2;
          const float facx = Exp(-qx);

          float4 color = vtx.color;
          color.w = facx*facy; /* alpha */

          const int idx = lineIdx + (idx - range.x0);
          using blend = Blending::getColor<Blending::ONE,Blending::SRC_ALPHA>;
          fb[idx] = blend(fb[idx], color);
        }
      }
      return q;
    }

    void render()
    {
      const int np = vtxArrayProj.size();

      image.resize(width*height);
      std::fill(image.begin(), image.end(), make_floa43(0.0f));

#define NTHREADMAX 256
      std::vector<float4> fbVec[NTHREADMAX];

#pragma omp parallel
      {
        const int nt = omp_get_num_threads();
        assert(nt <= NTHREADMAX);
        const int tid = omp_get_thread_num();
        auto &fb = fbVec[tid];
        fb.resize(width*height);
        std::fill(fb.begin(), fb.end(), make_float4(0.0f));

#pragma omp for schedule(runtime)
        for (int i = 0; i < np; i++)
          rasterize(vtxArrayProj[i], range, fb);

#pragma omp for schedule(runtime) collapse(2)
        for (int j = 0; j < height; j++)
          for (int i = 0; i < width; i++)
          {
            const int idx = j*width + i;
            for (int k = 0; k < nt; k++)
            {
              using blend = Blending::getColor<Blending::ONE,Blending::SRC_ALPHA>;
              image[idx] = blend(image[idx], fbVec[k][idx]);
            }
          }
      }
    }

    void finalize()
    {
#pragma omp for schedule(runtime) collapse(2)
      for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
        {
          const int idx = j*width + i;
          const float4 src = image[idx];
          float4 dst;

          dst.x = 1.0f - Exp(-src.x);
          dst.y = 1.0f - Exp(-src.y);
          dst.z = 1.0f - Exp(-src.z);
          dst.w = 1.0f;

          image[idx] = dst;
        }
    }

  public:
    void genImage(const bool perspective = true)
    {
      modelView();
      transform(perspective);
      depthSort();
      rasterize();
      finalize();
    }

};


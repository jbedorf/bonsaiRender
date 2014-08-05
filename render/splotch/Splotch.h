#pragma once
#include <omp.h>
#include <parallel/algorithm>


template<typename real_t>
struct Pos2D
{
  real_t x, y, h;
  Pos2D() {}
  Pos2D(const real_t &_x, const real_t &_y, const real_t &_h) : x(_x), y(_y), h(_h) {}
  bool isVisible() const { return h > 0.0f; }
};
template<typename real_t>
struct Pos3D
{
  real_t x,y,z,h;
  Pos3D() {}
  Pos3D(const real_t &_x, const real_t &_y, const real_t &_z, const real_t &_h) : x(_x), y(_y), z(_z), h(_h) {}
};
template<typename real_t>
struct Attribute
{
  real_t rho, vel, I;
  Attribute() {}
  Attribute(const real_t &_rho, const real_t &_vel, const real_t &_I) : rho(_rho), vel(_vel), I(_I) {}
};



template<typename Tpos, typename Tattr>
class VertexArrayT
{
  private:
    Tpos  *_pos;
    Tattr *_attr;
    int    _size;
    int    _capacity;
  public:
    struct Vertex
    {
      Tpos &pos;
      Tattr &attr;
      Attr(Tpos &_pos, Tattr &_attr) :
        pos(_pos) ,attr(_attr);
    };
    DataVector() : _pos(NULL), _attr(NULL), _size(0) {}
    DataVector(const Tpos *pos, const Tattr *attr, const int size) 
    {
      assert(size > 0);
      _size = size;
      _capacity = size;
      _pos  = new Tpos[_size];
      _attr = new Tattr[_size];
#pragma omp parallel for schedule(static)
      for (int i = 0; i < _size; i++)
      {
        _pos [i] = pos [i];
        _attr[i] = attr[i];
      }
    }
    ~DataVector() 
    {
      if (_size > 0)
      {
        delete[] _pos;
        delete[] _attr;
        _pos  = NULL;
        _attr = NULL;
      }
    }
    Tpos& pos(const int i) {return _pos[i];}
    const Tpos& pos(const int i) const {return _pos[i];}
    Tattr& attr(const int i) {return _attr[i];}
    const Tattr& attr(const int i) const {return _attr[i];}

    Vertex& operator[](const int i) {return Vertex(_pos[i], _attr[i]);}
    const Vertex& operator[](const int i) const {return const Attr(_pos[i],_attr[i]);}
    int size() const {return _size;}
};

class Splotch
{
  public:
    using pos2d_t = Pos2D<float>;
    using pos3d_t = Pos3D<float>;
    using attr_t  = Attribute<float>;  

  private:

    using VertexArray     = VertexArrayT<pos3d_t,attr_t>;
    using VertexArrayProj = VertexArrayT<pos2d_t,attr_t>;
    using Vertex          = VertexArray::Vertex;
    using VertexProj      = VertexArrayProj::Vertex;
    using Exp     = std::exp;

    VertexArray     vtxArray;
    VertexArrayProj vtxArrayProj, vtxArrayView;
    real2_t invProjRange;
    
    struct Quad
    {
      float x0,x1;
      float y0,y1;
    };
    
    enum BlendType {BLEND_ONE, BLEND_ZERO, BLEND_SRC_ALPHA, BLEND_ONE_MINUS_SRC_ALPHA};

    template<BlendType BLEND>
    float4 BLendColorT(const float4 col, const float4 src, const float4 dst)
    {
      float4 res;
      switch (BLEND)
      {
        case BLEND_ONE_MINUS_SRC_ALPHA:
          res.x *= 1.0f - src.w;
          res.y *= 1.0f - src.w;
          res.z *= 1.0f - src.w;
          res.w *= 1.0f - src.w;
        case BLEND_SRC_ALPHA:
          res.x *= src.w;
          res.y *= src.w;
          res.z *= src.w;
          res.w *= src.w;
        case BLEND_ZERO:
          res.x = res.y = res.z = rez.w = 0;
          break;
        case BLEND_ONE:
          break;
        default:
          assert(0);
      }
      return res;
    }

    template<BlendType SRC, BlendType DST>
    float4 BlendT(const float4 d, const float4 s)
    {
      float4 res;


      const float4 src = BlendColorT<SRC>(s,d,s);
      const float4 dst = BlendColorT<DST>(d,d,s);

      res.x = src.x + dst.x;
      res.y = src.y + dst.y;
      res.z = src.z + dst.z;
      res.w = src.w + dst.w;

      using max = std::max;
      res.w = min(1.0f, res.w);

      return res;
    }

    int width, height;
    std::vector<float4> image;

  public:
    Splotch() {}
    ~Splotch() {}

    void modelView(const bool project = true)
    {
      const int np = vtxArray.size();
      vtxArrayView = VertexArrayProj(np);
      vtxArrayProj = VertexArrayProj(np);

      int nActive = 0;
#pragma omp parallel for schedule(runtime) reduction(+:nActive)
      for (int i = 0; i < np; i++)
      {
        const auto &vtx = vtxArray[i];
        const float4 posO(vtx.pos.x,vtx.pos.y,vtx.pos.z,1.0f);
        float4 posV = modelViewMatrix * posO;
        posV.w = -1.0f;
        float4 posP = posV;

        if (posT.z > depthMin && posT.z > depthMax)
        {
          const float dist = project ? posT.z : zDist;
          const float xfac = project ? 1.0f/(fovfct*posT.z)  : scaleFactor;
          posV.x = res2*(posV.x+fovfct*zDist)*scaleFactor;                                          
          posV.y = res2*(posV.y+fovfct*zDist)*scaleFactor + yCorr;   
          posV.w = posO.h * res2*scaleFactor;

          posP = posV;

          if (project)
          {
            const float xfac = 1.0f/(fovfct*posT.z);
            posP.x = res2*(posP.x+fovfct*posP.z)*xfac;
            posP.y = res2*(posP.y+fovfct*posP.z)*xfac + yCorr;   
            posP.w = posO.h * res2*xfac;

            const float rcorr = std::sqrt(posP.w*posP.w + minHpix*minHpix)/posP.w;
            posP.w *= rcorr;
          }
            
          if ( posP.x - posP.w > width
            || posP.x + posP.w < 0
            || posP.y - posP.w > height
            || posP.y + posP.w < 0)
            posV.w = posP.w = -1.0;
        }

        Vertex vtxProj;
        vtxProj.pos  = pos2d_t(posP.x,posP.y,posP.w);
        vtxProj.attr = vtx.attr;
        vtxArrayProj[i] = vtxProj;
       
        Vertex vtxView;
        vtxView.pos  = pos2d_t(posV.x,posV.y,posV.w);
        vtxView.attr = vtx.attr;
        vtxArrayView[i] = vtxView;

        nVisible += vtxView.pos.isVisible();
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
        if (posView.pos.isVisible())
          depthMap.push_back(std::make_pair(posView.pos.z,i));
      
      __gnu_parallel::sort(depthMap.begin(), depthMap.end(),
          [](const pair &a, const pair &b) { return a.first < b.first;} );

      const int npVis = depthMap.size();
      VertexArrayProj vtxProj(npVis);

#pragma omp parallel for 
      for (int i = 0; i < npVis; i++)
      {
        const auto &map = depthMap[i];
        vtxProj[i] = vtxArrayProj[map.second];
      }

      vtxArrayProj = vtxProj;
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

          const float3 col3 = assignColor(vtx.attr.rho, vtx.attr.vel);
          float4 color;
          color.x = col3.x;
          color.y = col3.y;
          color.z = col3.z;
          color.w = facx*facy; /* alpha */

          const int idx = lineIdx + (idx - range.x0);
          using blend = BlendT<BLEND_ONE,BLEND_SRC_ALPHA>;
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
              using blend = BlendT<BLEND_ONE,BLEND_SRC_ALPHA>;
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

    void genImage()
    {
      modelView();
      depthSort();
      rasterize();
      finalize();
    }

};


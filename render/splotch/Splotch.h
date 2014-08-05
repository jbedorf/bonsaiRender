#pragma once
#include <omp.h>


template<typename real_t>
struct Pos2D
{
  real_t x, y, h;
  Pos2D() {}
  Pos2D(const real_t &_x, const real_t &_y, const real_t &_h) : x(_x), y(_y), h(_h) {}
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
  real_t rho, vel;
  Attribute() {}
  Attribute(const real_t &_rho, const real_t &_vel) : rho(_rho), vel(_vel) {}
};



template<typename Tpos, typename Tattr>
class VertexArrayT
{
  private:
    Tpos  *_pos;
    Tattr *_attr;
    int    _size;
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
    using real_t  = float;

    using pos2d_t = Pos2D<real_t>;
    using pos3d_t = Pos3D<real_t>;
    using attr_t  = Attribute<real_t>;  /* radius, density, velocity */

  private:

    using VertexArray     = VertexArrayT<pos3d_t,attr_t>;
    using VertexArrayProj = VertexArrayT<pos2d_t,attr_t>;
    using Vertex          = VertexArray::Vertex;
    using VertexProj      = VertexArrayProj::Vertex;
    using Exp     = std::exp;

    VertexArray     vtxArray;
    VertexArrayProj vtxArrayProj;
    real2_t invProjRange;
    
    struct Quad
    {
      real_t x0,x1;
      real_t y0,y1;
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

  public:
    Splotch() {}
    ~Splotch() {}


    Quad rasterize(const VertexProj &vtx, const real_t width, const real_t height, Vector<float4> &fb)
    {
      using max = std::max;
      using min = std::min;
      Quad q;
      q.x0  = max(0.0f,   vtx.pos.x - vtx.pos.h);
      q.x1  = mix(width,  vtx.pos.x + vtx.pos.h);
      q.y0  = max(0.0f,   vtx.pos.y - vtx.pos.h);
      q.y1  = max(height, vtx.pos.y + vtx.pos.h);

      int lineIdx = q.y0*width;
      const real_t invh  = 1.0f/vtx.pos.h;
      const real_t invh2 = invh*invh;
      for (real_t iy = q.y0; iy < q.y1; iy++, lineIdx += width)
      {
        const real_t dy   = iy - pos.vtx.y;
        const real_t qy   = dy*dy * invh2;
        const real_t facy = Exp(-qy);
        for(real_t ix = q.x0; ix < q.x1; ix++)
        {
          const real_t dx = ix - pos.vtx.x;
          const real_t qx = dx*dx * invh2;
          const real_t facx = Exp(-qx);

          const float3 col3 = assignColor(vtx.attr.rho, vtx.attr.vel);
          float4 color;
          color.x = col3.x;
          color.y = col3.y;
          color.z = col3.z;
          color.w = facx*facy; /* alpha */

          using blend = BlendT<BLEND_ONE,BLEND_SRC_ALPHA>;
          fb[lineIdx+ix] = blend(fb[lineIdx+ix], color);
        }
      }
      return q;
    }

    void render(const int height, const int width, std::vector<float4> &image)
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
          rasterize(vtxArrayProj[i], height, widhth, fb);

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

    void finalize(const int height, const int width, std::vector<float4> &image)
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




};


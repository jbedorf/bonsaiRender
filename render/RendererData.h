#pragma once

struct float3 { float x, y, z; };
struct float4 { float x,y,z,w; };
struct int2   { int x, y; };

static int2 make_int2(int _x, int _y) 
{ 
  int2 v; v.x = _x; v.y = _y;; return v;
}
static float3 make_float3(float x, float y, float z)
{
  float3 v; v.x = x; v.y = y; v.z=z; return v;
}
static float4 make_float4(float x, float y, float z, float w)
{
  float4 v; v.x = x; v.y = y; v.z=z; v.w=w; return v;
}

class RendererData
{
  public:
    typedef unsigned long long long_t;
    enum Attribute_t {
      MASS,
      VEL,
      NPROP};
  private:
    const int _n;
    float *_posx, *_posy, *_posz;
    long_t *_ID;
    int   *_type;
    float *_attribute[NPROP];

  public:
    RendererData(const int __n) : _n(__n)
  {
    _posx   = new float[_n];
    _posy   = new float[_n];
    _posz   = new float[_n];
    _ID     = new long_t[_n];
    _type   = new int  [_n];

    for (int p = 0; p < NPROP; p++)
      _attribute[p] = new float[_n];
  }

    ~RendererData()
    {
      delete[] _posx;
      delete[] _posy;
      delete[] _posz;
      delete[] _ID;
      delete[] _type;
      for (int p = 0; p < NPROP; p++)
        delete[] _attribute[p];
    }

    const int n() const { return _n; }

    float  posx(const int i) const { return _posx[i]; }
    float& posx(const int i)       { return _posx[i]; }
    float  posy(const int i) const { return _posy[i]; }
    float& posy(const int i)       { return _posy[i]; }
    float  posz(const int i) const { return _posz[i]; }
    float& posz(const int i)       { return _posz[i]; }

    float  attribute(const Attribute_t p, const int i) const {return _attribute[p][i]; }
    float& attribute(const Attribute_t p, const int i)       {return _attribute[p][i]; }

    int  type(const int i) const { return _type[i]; }
    int& type(const int i)       { return _type[i]; }

    long_t  ID(const long_t i) const { return _ID[i]; }
    long_t& ID(const long_t i)       { return _ID[i]; }
};

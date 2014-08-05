#pragma once

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

template<typename Tpos, typename Tattr, typename Tcolor)
class VertexArrayT
{
  private:
    Tpos  *_pos;
    Tcolor *_color;
    Tattr *_attr;
    int    _size;
  public:
    struct Vertex
    {
      Tpos pos;
      Tattr attr;
      Tcolor color;
      Vertex(Tpos &_pos, Tattr &_attr, Tcolor &_color) :
        pos(_pos), color(_color), attr(_attr) {}
      bool isVisible() const { return pos.isVisible(); }
    };
    struct VertexRef
    {
      Tpos &pos;
      Tcolor &color;
      Tattr &attr;
      VertexRef(Tpos &_pos, Tattr &_attr, Tcolor &_color) :
        pos(_pos), color(_color), attr(_attr) {}
      VertexRef operator=(const Vertex &v) 
      {
        pos   = v.pos;
        color = v.color;
        attr  = v.attr;
      }
      bool isVisible() const {return pos.isVisible(); }
    };

  private:
    void free()
    {
      if (_size > 0)
      {
        ::free(_pos  ); 
        ::free(_color);
        ::free(_attr );
        _size = 0;
      }
      _pos   = NULL;
      _color = NULL;
      _attr  = NULL;
    }
    void malloc(const int size)
    {
      free();
      _size = size;
      if (_size > 0)
      {
        _pos   = (Tpos  *)::malloc(sizeof(Tpos  )*_size);
        _color = (Tcolor*)::malloc(sizeof(Tcolor)*_size);
        _attr  = (Tattr *)::malloc(sizeof(Tattr )*_size);
      }
    }

  public:
    VertexArrayT(const int size = 0) : _pos(NULL), _color(NULL), _attr(NULL), _size(size) { malloc(size); }
    VertexArrayT(const Tpos *pos, const Tcolor, *color, const Tattr *attr, const int size) 
    {
      malloc(size);
#pragma omp parallel for schedule(static)
      for (int i = 0; i < _size; i++)
      {
        _pos  [i] = pos  [i];
        _color[i] = color[i];
        _attr [i] = attr [i];
      }
    }
    ~VertexArrayT()  { free(); }

    VertexArrayT& operator=(const VertexArrayT &data)
    {
      free();
      *this = VertexArrayT(data._pos, data._color, data._attr, data._size);
      return *this;
    }

    friend static void swap(VertexArrayT &a, const VertexArrayT &b)
    {
      std::swap(a._pos,   b._pos);
      std::swap(a._color, b._color);
      std::swap(a._attr,  b._attr);
      std::swap(a._size,  b._size);
    }


    VertexRef operator[](const int i) {return Vertex(_pos[i], _attr[i]);}
    const VertexRef operator[](const int i) const {return const Attr(_pos[i],_attr[i]);}
    int size() const {return _size;}
};

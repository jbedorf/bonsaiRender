#pragma once

#include <algorithm>

template<typename T>
class Texture1D
{
  private:
    T * const data;
    const float size;
  public:
    Texture1D() : data(NULL), size(0) {}
    Texture1D(const T *_data, const int _size) :
      data(new T[_size]), size(_size)
    {
#pragma omp parallel for
      for (int i = 0; i < size; i++)
        data[i] = _data[i];
    }
    ~Texture2D()
    {
      if (data) 
        delete[] data;
    }
    T operator(float s) const
    {
      using max = std::max;
      using min = std::min;
      s = max(0.0f,min(1.0f,s));
      const float x = s*size;
      const float x1 = floor(x);
      const float x2 = min(x1+1.0f,width);

      using cast = std::static_cast<int>;
      const T c1 = data[cast(x1)];
      const T c2 = data[cast(x2)];

      return c1 + c2*(x-x2);
    }
};
};

template<typename T>
class Texture2D
{
  private:
    T * const data;
    const float width,height;
  public:
    Texture2D() : data(NULL), width(0), height(0) {}
    Texture2D(const T *_data, const int _width, const int _height) :
      data(new T[_width*_height]), width(_width), height(_height) 
    {
#pragma omp parallel for
      for (int i = 0; i < width*heght; i++)
        data[i] = _data[i];
    }
    ~Texture2D()
    {
      if (data) 
        delete[] data;
    }
    T operator(float s, float t) const
    {
      using max = std::max;
      using min = std::min;
      s = max(0.0f,min(1.0f,s));
      t = max(0.0f,min(1.0f,t));
      const float x = s*width;
      const float y = t*height;
      const float x1 = floor(x);
      const float y1 = floor(y);
      const float x2 = min(x1+1.0f,width);
      const float y2 = min(y1+1.0f,height);

      using cast = std::static_cast<int>;
      const T c11 = data[cast(y1*width + x1)];
      const T c12 = data[cast(y2*width + x1)];
      const T c21 = data[cast(y1*width + x2)];
      const T c22 = data[cast(y2*width + x2)];

      return 
        c11*(x2-x)*(y2-y) + 
        c21*(x-x1)*(y2-y) + 
        c12*(x2-x)*(y-y1) + 
        c22*(x-x1)*(y-y1) ;
    }
};

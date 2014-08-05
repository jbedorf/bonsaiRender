#pragma once

template<typename T, int N>
class MathArray
{
  private:
    T data[N];

  public:

    MathArray() {}
    MathArray(const T &v)
    {
      for (auto &x : data) 
        x = v;
    }
    T& operator[](const int i)       { return data[i]; }
    T  operator[](const int i) const {return data[i]; }

#define OPERATOR(OP) { \
  MathArray& operator##OP##(const MathArray& rhs)   \
  {  \
    foreach([](int i){ data[i] ##OP## rhs[i]; });  \
    return *this;  \
  } 
  OPERATOR(*=);
#undef METHOD

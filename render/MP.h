#pragma once


#include <iostream>

namespace MP
{

  struct float4
  {
    typedef float  v4sf __attribute__ ((vector_size(16)));
    typedef double v2df __attribute__ ((vector_size(16)));
    static v4sf v4sf_abs(v4sf x){
      typedef int v4si __attribute__ ((vector_size(16)));
      v4si mask = {0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff};
      return __builtin_ia32_andps(x, (v4sf)mask);
    }
    union{
      v4sf v;
      struct{
        float x, y, z, w;
      };
    };
    float4() : v((v4sf){0.f, 0.f, 0.f, 0.f}) {}
    float4(float x, float y, float z, float w) : v((v4sf){x, y, z, w}) {}
    float4(float x) : v((v4sf){x, x, x, x}) {}
    float4(v4sf _v) : v(_v) {}
    float4 abs(){
      typedef int v4si __attribute__ ((vector_size(16)));
      v4si mask = {0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff};
      return float4(__builtin_ia32_andps(v, (v4sf)mask));
    }
    void dump(){
      std::cerr << x << " "
        << y << " "
        << z << " "
        << w << std::endl;
    }
#if 1
    v4sf operator=(const float4 &rhs){
      v = rhs.v;
      return v;
    }
    float4(const float4 &rhs){
      v = rhs.v;
    }
#endif
#if 0
    const v4sf_stream operator=(const v4sf_stream &s){
      __builtin_ia32_movntps((float *)&v, s.v);
      return s;
    }
    float4(const v4sf_stream &s){
      __builtin_ia32_movntps((float *)&v, s.v);
    }
#endif
  };


  using vector3 =  std::array<float,3>;

  // void MP_initialize(int argc, char *argv[]);
  void MP_initialize(int *argc, char ***argv);
  void MP_end();
  void MP_sync();
  int MP_myprocid();
  int MP_proccount();
  void MP_copyparams(float &dt,
      float &dtsnapout,
      int &outlogstep,
      float &tend,
      float &eps,
      float &theta,
      int &ncrit,
      float &pos_scale,
      float &vel_scale);
  void MP_convert_snap_name(int& flag, char * name);

  void MP_gather_sample_coords(int&nsample, vector3 * sample_array);
  void MP_gather_sample_coords(int &nsample, std::vector<float4> &sample_array);
  void MP_int_bcast(int&i);
  void MP_int_sum(int&i);
  void MP_sum(double& r);
  void MP_double_sum(double r[], int count, bool allreduce);
  void MP_double_max(double r[], int count, bool allreduce);
  void MP_double_bcast(double*i, int nwords);
#if 0
  void MP_exchange_particle(int ibox,
      nbody_particle * pb,
      int firstloc,
      int nparticles,
      int isource,
      int &iloc);
  int MP_exchange_particle_with_overflow_check(int ibox,
      nbody_particle * pb,
      int firstloc,
      int nparticles,
      int isource,
      int &iloc,
      int &nsend);
#endif

  int MP_intmax(int localval);
  double MP_doublemax(double localval);
  void MP_exchange_bhlist(int ibox,
      int nlist,
      int nbmax,
      vector3 * plist,
      float * mlist,
      int isource,
      int & nrecvlist,
      vector3 * precvbuf,
      float * mrecvbuf);
  void MP_Abort(int);

  template <typename T>
    void MP_alltoallv(std::vector<T> sendbuf[], std::vector<T> recvbuf[]);

#if 0
  void MP_exchange_particle_alltoall(
      nbody_particle *sendbuf,
      nbody_particle *recvbuf,
      int sendoff[],
      const int nrecv_max,
      int &nrecv);
#endif

  void MP_collect_cmterms(vector3& pos,vector3& vel,float& mass);

  void MP_print_times(std::ostream &s);
  void MP_print_treestats(float total_interactions,
      int tree_walks,
      int nisum,
      std::ostream &s);
  void MP_print_float(
      const char *prefix,
      float val,
      std::ostream &ofs,
      bool print_maxmin);
#if 0
  void MP_print_string(std::string &str, std::ostream &os, const char *format);
  static inline void MP_print_string(std::stringstream &ss, std::ostream &os, const char *format){
    std::string s = ss.str();
    MP_print_string(s, os, format);
  }
#endif

  void MP_bool_bcast(bool& i);
  void MP_char_bcast(char * name, int nword);
  void MP_string_bcast(std::string &str);
  void MP_print_md5(unsigned char md5[16]);
  bool MP_root();
  const char *MP_get_hostname(int);
  const char *MP_get_hostname();
  template <typename T>
    void MP_copy_vector(std::vector<T> &vec, const int dst, const int src);

}


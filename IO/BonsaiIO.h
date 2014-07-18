#pragma once

namespace BonsaiIO
{
  typedef float  fp32_t;
  typedef double fp64_t;

  /*********** Bonsai Header *************/

  class Core
  {
    std::vector<std::pair<uint64_t,DataStructBase*>> data;
  }

  class BonsaiHeader
  {
    const char versionString[256];
    fp64_t time;
    fp64_t dummySpace[256];
    std::vector<int> npPerType;
  };
  
  /*********** Bonsai Data Types *************/

  enum TypeName 
  {
    IDTYPE, POS, VEL,
    RHOH, GASDATA,
    TREENODE, TREESTRUCTURE,
    NTYPES
  };

  /*********** Bonsai Data Type Base class *************/

  class DataStructBase
  {
    public:
      virtual TypeName getTypeName() const = 0;
  };
  
  /*********** Bonsai Data Type Derived Classes  *************/

  class IDType_t : public DataStructBase
  {
    private:
      uint64_t _IDTypePacked;
    public:
      virtual TypeName getTypeName() const { return IDTYPE; }

      uint64_t getID() const
      {
        return _IDTypePacked & 0xFFFF000000000000ULL;
      }
      uint32_t getType() const
      {
        return static_cast<uint32_t>(_IDTypePacked >> 48);
      }
      void setID(const int64_t ID)
      {
        const uint32_t type = getType();
        _IDTypePacked = (ID & 0xFFFF000000000000ULL) | (static_cast<uint64_t>(type) << 48);
      }
      void setType(const int type)
      {
        const uint64_t ID = getID();
        _IDTypePacked  = ID | (static_cast<uint64_t>(type) << 48);
      }
  };

  class Position : public DataStructBase
  {
    private:
      fp32_t _pos[3];
    public:
      virtual TypeName getTypeName() const { return POS; }
      fp32_t  pos(const int i) const { return _pos[i]; }
      fp32_t& pos(const int i)       { return _pos[i]; }
  };
  class Velocity : public DataStructBase
  {
    private:
      fp32_t _vel[3];
    public:
      virtual TypeName getTypeName() const { return VEL; }
      fp32_t  vel(const int i) const { return _vel[i]; }
      fp32_t& vel(const int i)       { return _vel[i]; }
  };
  class DensityH : public DataStructBase
  {
    private:
      fp32_t _density;
      fp32_t _h;
    public:
      virtual TypeName getTypeName() const { return VEL; }
      fp32_t  _density() const { return _density; }
      fp32_t& _density()       { return _density; }
      fp32_t  _h() const { return _h; }
      fp32_t& _h()       { return _h; }
  };
  
  class GasData : public DensityH
  {
    private:
      fp32_t _temperature, pressure;
    public:
      virtual TypeName getTypeName() const { return VEL; }
      fp32_t  _pressure() const { return _density; }
      fp32_t& _pressure()       { return _density; }
      fp32_t  _pressure() const { return _pressure; }
      fp32_t& _pressure()       { return _pressure; }
  };


  class TreeNode : public DataStructBase
  {
    private:
      fp32_t rmin[3];
      fp32_t rmax[3];
      fp32_t _mass, _com[3];
    public:
      virtual TypeName getTypeName() const { return TREENODE; }
      fp32_t  min(const int i) const { return rmin[i]; }
      fp32_t& min(const int i)       { return rmin[i]; }
      fp32_t  max(const int i) const { return rmax[i]; }
      fp32_t& max(const int i)       { return rmax[i]; }
      fp32_t  com(const int i) const { return _com[i]; }
      fp32_t& com(const int i)       { return _com[i]; }
      fp32_t  mass() const { return _mass; }
      fp32_t& mass()       { return _mass; }
  };

  class TreeStructure : public DataStructBase
  {
    private:
      uint64_t _parent;
      uint64_t _children[8];
    public:
      virtual TypeName getTypeName() const { return TREESTRUCTURE; }
      uint64_t  parent() const { return _parent; }
      uint64_t& parent()       { return _parent; }
      uint64_t  children(const int i) const { return _children[i]; }
      uint64_t& children(const int i)       { return _children[i]; }
  };
};

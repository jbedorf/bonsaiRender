#pragma once

namespace BonsaiIO
{
  typedef long long int long_t;

  /*********** Header ******************/
  class Header
  {
    private:
      const std::string versionString;
      double time;

    public:
    Header() : versionString("V1")
    {
    }

  };

  
  /*********** Data types *************/

  class DataTypeBase
  {
    private:
      std::string name;
      long_t n;
      int   elementSize;
    public:
      virtual void DataTypeBase(const std::string _name, const int n) : name(_name), n(_n) {};
      virtual void ~DataTypeBase() {};
      const std::string& typeName() const {return name;}
  };


  template<typename T>
    class DataType : public DataTypeBase
  {
    private:
      std::vector<T> data;
    public:
      virtual void DataType(std::string _name, const int n = 0) : DataTypeBase(_name, _n)
    {
      if (n > 0)
        data.resize(n);
      dataSize = sizeof(T);
    }
      const T& operator[](const size_t i) const { return data[i]; }
      T operator[](const size_t i) const { return data[i]; }
  };

  /*********** Core reader/writer *************/

  enum IOTYPE 
  {
    FREAD, FWRITE,         /* one file per MPI proc, uses default fread/fwrite */
    READMPI, WRITEMPI      /* use MPI-IO */
  }

  class Core
  {
    private:
      const int myRank;
      const int nRank;
      const MPIComm &comm;
      const IOTYPE iotype;
      const std::string fileName;

      MPI_File fh;
      MPI_Status mpi;

      Header header;
      

    public:
        Core(const int _myRank,
            const int _nRank,
            const MPI_Comm &_comm,
            const IOTYPE _iotype,
            const std::string _fileName) :
          write(false)
          myRank(_myRank),
          nRank(_nRank),
          comm(_comm),
          iotype(_iotype),
          fileName(_fileName)
    {
      int file_open_error;
      switch (iotype)
      {
        case WRITEMPI:
          MPI_File_open(
              comm,
              fileName,c_str(),
              MPI_MODE_CREATE | MPI_MODE_RDWR,
              MPI_INFO_NULL,
              &fn);
          break;
        case READMPI:
          file_open_error = MPI_File_open(
              comm,
              fileName,c_str(),
              MPI_MODE_RDONLY,
              MPI_INFO_NUL,
              &fh);
          if (file_open_error != MPI_SUCCESS)
          {
            char error_string[256];
            int length_of_error_string, error_class;

            MPI_Error_class(file_open_error, &error_class);
            MPI_Error_string(error_class, error_string, &length_of_error_string);
            printf("%3d: %s\n", my_rank, error_string);

            MPI_Error_string(file_open_error, error_string, &length_of_error_string);
            printf("%3d: %s\n", my_rank, error_string);

            MPI_Abort(MPI_COMM_WORLD, file_open_error);
          }
          break;
      }
    }


        bool readMPIIO(DataTypeBase &data)
        {
          assert(isRead());
          if (!header.seek(data.getName()))
            return false;

          const size_t nBytesGbl = header.getNBytes();


          MPI_Offset dataOffset = header.getOffset();

          const uint64_t nBytesPerRank = (nBytesGlb - 1 + nRank) / nRank;
          const uint64_t beg = myRank * nBytesPerRank;
          const uint64_t end = std::min(beg + nBytesPerRank, nBytesGlb);
          
          data.allocate(nBytesPerRank);
          
          uint64_t nBytes = end - beg;
          const uint64_t nMaxPerRead = (1U << 31) - 1;
          while (nBytes > 0)
          {
            const int count = static_cast<int>(std::min(nBytes, nMaxPerWrite));
            assert(count > 0);
            MPI_File_read(
                fh,
                data.getDataPtr(),
                count, 
                MPI_BYTE, 
                &status);
            int read;
            MPI_Get_count(&status, MPI_INT, &read);
            assert(count == read);
            nBytes -= nMaxPerRead;
          }
          return true;
        }

        bool writeMPIIO(const DataTypeBase &data)
        {
          assert(isWrite());
          const double tWrite = MPI_Wtime();


          uint64_t nBytesLoc = data.getN() * data.unitSize();
          uint64_t nBytesGlb;
          

          MPI_Allreduce(&nBytesLoc, &nBytesGlb, 1, MPI_LONG_LONG, MPI_SUM, comm);
          header.add(data.getName(), nBytesGlb);

          const uint64_t nBytesPerRank = (nBytesGlb  - 1 + nRank) / nRank;
          const uint64_t beg = myRank * nBytesPerRank;
          const uint64_t end = std::min(beg + nBytesPerRank, nBytesGlb);

          MPI_Offset myOffset = headerOffset + dataOffset + beg;
          MPI_File_seek(fh, myOffset, MPI_SEEK_SET);

          uint64_t nBytes = end - beg;
          const uint64_t nMaxPerWrite = (1U << 31) - 1;
          while (nBytes > 0)
          {
            const int count = static_cast<int>(std::min(nBytes, nMaxPerWrite));
            assert(count > 0);
            MPI_File_write(
                fh,
                data.getDataPtr(),
                count, 
                MPI_BYTE, 
                &status);
            int written;
            MPI_Get_count(&status, MPI_INT, &written);
            assert(count == written);
            nBytes -= nMaxPerWrite;
          }

          dataOffset += nBytesGlb;

          dtWrite += MPI_Wtime() - tWrite;
          return true;
        }


        void close()
        {
          /* write header */
          MPI_File_close(&fh);
        }

  };

  class Header
  {
    const char versionString[256];
    fp64_t time;
    fp64_t dummySpace[256];
    std::vector<std::pair<uint64_t,uint64_t>> typeInfo
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

  class IDType : public DataStructBase
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

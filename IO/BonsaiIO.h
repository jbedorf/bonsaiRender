#pragma once

namespace BonsaiIO
{
  struct Exception : public std::exception
  {
    std::string s;
    Exception(std::string ss) : s(ss) {}
    ~Exception() throw () {} // Updated
    const char* what() const throw() { return s.c_str(); }
  };

  typedef long long int long_t;
  
  /*********** Data types *************/

  class DataTypeBase
  {
    private:
      std::string name;
    public:
      void DataTypeBase(const std::string _name) : name(_name) {};
      const std::string& getName() const {return name;}

      virtual size_t getUnitSize   () const = 0;
      virtual size_t getNumElements() const = 0;
      virtual size_t getNumBytes   () const = 0;
  };

  template<typename T>
    class DataType : public DataTypeBase
  {
    private:
      std::vector<T> data;
    public:
      DataType(std::string name, const int _n = 0) : DataTypeBase(name)
      {
        if (n > 0)
          data.resize(n);
      }
      void allocate(const size_t n)
      {
        data.resize(n);
      }
      const T& operator[](const size_t i) const { return data[i]; }
            T& operator[](const size_t i)       { return data[i]; }

      size_t getUnitSize   () const {return sizeof(T);}
      size_t getNumElements() const {return data.size();}
      size_t getNumBytes   () const {return data.size()*sizeof(T);}
  };

  /*********** Header ******************/

  class Header
  {
    private:
      struct dataType
      {
        std::string name;
        size_t unitSize;
        size_t numElements;
      };

      const std::string versionString;
      double time;
      std::vector<dataType> data;

    public:
      Header() : versionString("V1") {}

      int find(const std::string &name) const
      {
        const int n = data.size();
        for (int i = 0; i < n; i++)
          if (data[i].name == name)
            return i;
        return -1;
      }

      MPI_Offset getOffset(const int idx) const
      {
        assert(idx >= 0);
        MPI_Offset myOffset = headerSize();
        for (int i = 0; i < idx-1; i++)
          myOffset += data[i].unitSize * data[i].numElements;
        return myOffset;
      }

      bool add(const DataTypeBase &dataVec)
      {
        if (find(name) == -1)
        {
          dataType d;
          d.name        = dataVec.getName       ();
          d.unitSize    = dataVec.getUnitSize   ();
          d.numElements = dataVec.getNumElements();
          data.push_back(d);
          return true;
        }

        return false;
      }

      bool read(MPI_File &fh)
      {
      }

      bool write(MPI_File &fh)
      {
      }

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
        /* make sure we are in the reading phase */
        if (!isRead())
          throw Exception("Trying to read a file, while using WRITE mode.");

        /* find data set */
        const int idx = header.find(data.getName());
        if (idx == -1)
          return false;

        /* confirm that unit size has the same length */
        if (header.getUnitSize(idx) != data.getUnitSize())
          return false

        const size_t nElementsGlb = header.getNumBytes(idx);

        const uint64_t nElementsPerRank = (nElementsGlb - 1 + nRank) / nRank;
        const uint64_t beg = myRank * nElementsPerRank;
        const uint64_t end = std::min(beg + nElementsPerRank, nElementsGlb);

        data.allocate(nElementsPerRank);

        uint64_t nBytes = (end - beg) * data.getUnitSize();
        const uint64_t nMaxPerRead = (1U << 31) - 1;

        const MPI_Offset dataOffset = header.getDataOffset(dx);
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
          if (count != read)
          
            throw Exception("Count != read while reading the file.");
          nBytes -= nMaxPerRead;
        }
        return true;
      }

      bool writeMPIIO(const DataTypeBase &data)
      {
        /* make sure we are in the writing phase */
        if (!isWrite())
          throw Exception("Trying to write a file, while using READ mode.");
        const double tWrite = MPI_Wtime();

        uint64_t nBytesLoc = data.getN() * data.unitSize();

        MPI_Allreduce(&nBytesLoc, &nBytesGlb, 1, MPI_LONG_LONG, MPI_SUM, comm);
        if (isMaster())
          header.add(data);

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
        if (isMaster())
        {
          MPI_FIle_seek(fh, 0, MPI_SEEK_SET);
          header.write(&fh);
        }
        MPI_Barrier(comm);
        /* write header */
        MPI_File_close(&fh);
      }

  };

}

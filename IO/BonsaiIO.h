#pragma once

#include <mpi.h>
#include <string>
#include <vector>
#include <exception>
#include <cassert>
#include <cstdlib>


namespace BonsaiIO
{
  /*#####################################*/
  /*#####################################*/
  /*#####################################*/

  struct Exception : public std::exception
  {
    std::string s;
    Exception(std::string ss) : s(ss) {}
    ~Exception() throw () {} // Updated
    const char* what() const throw() { return s.c_str(); }
  };

  /*#####################################*/
  /*#####################################*/
  /*#####################################*/
  
  enum IOTYPE 
  {
    READ, WRITE
  };

  typedef long long int long_t;
#define MPI_LONGT MPI_LONG_LONG

  /*#####################################*/
  /*#####################################*/
  /*#####################################*/

  /* Base IO class */
  class FileIO
  {
    protected:
      bool _opened;
      IOTYPE iotype;
    public:
      FileIO() : _opened(false) {}
      virtual ~FileIO() {}
      virtual bool isOpened() const { return _opened; }
      virtual bool isRead  () const { return isOpened() && iotype == READ;  }
      virtual bool isWrite () const { return isOpened() && iotype == WRITE; }

      virtual void open(const std::string&, IOTYPE) = 0;
      virtual void close() = 0;
      virtual void seek (const size_t offset) = 0;
      virtual void read (void *data, const size_t count, const std::string &errString) = 0;
      virtual void write(const void *data, const size_t count, const std::string &errString) = 0;
  };

  /* MPI-IO implementation */
  class MPIFileIO : public FileIO
  {
    private:
      const MPI_Comm &comm;
      MPI_File fh;
      MPI_Status status;
      MPI_Offset offset;
      int checkCount;

    public:
      MPIFileIO(const MPI_Comm &_comm) : FileIO(), comm(_comm) {}
      virtual ~MPIFileIO() {}

      void open(const std::string &fileName, const IOTYPE _iotype)
      {
        assert(!isOpened());
        iotype = _iotype;
        if (iotype == READ)
        {
          const int err = 
            MPI_File_open(
                comm,
                (char*)fileName.c_str(),
                MPI_MODE_RDONLY,
                MPI_INFO_NULL,
                &fh);
          if (err != MPI_SUCCESS)
            throw Exception("Unable to open a file to read.");
        }
        else if (iotype == WRITE)
        {
          const int err = 
            MPI_File_open(
                comm,
                (char*)fileName.c_str(),
                MPI_MODE_CREATE | MPI_MODE_RDWR,
                MPI_INFO_NULL,
                &fh);

          if (err != MPI_SUCCESS)
            throw Exception("Unable open a file to write.");
        }
        else
          assert(0);
        _opened = true;
      }
      void close() 
      {
        assert(isOpened());
        MPI_File_close(&fh);
        _opened = false;
      }
      void seek(const size_t offset)
      {
        MPI_File_seek(fh, offset, MPI_SEEK_SET);
      }

      void read(void *data, const size_t count, const std::string &errString) 
      {
        assert(isRead());
        const size_t batchMax = (1U << 31) - 1;
        size_t offset = 0;
        size_t nRead = count;
        while (nRead > 0)
        {
          const int count = static_cast<int>(std::min(nRead, batchMax));
          assert(count > 0);
          MPI_File_read(
              fh,
              (void*)(reinterpret_cast<char*>(data) + offset),
              count, 
              MPI_BYTE, 
              &status);
          MPI_Get_count(&status, MPI_BYTE, &checkCount);
          if (count != checkCount)
            throw Exception(errString);
          nRead  -= count;
          offset += count;
        }
      }

      void write(const void *data, const size_t count, const std::string &errString)
      {
        assert(isWrite());
        const size_t batchMax = (1U << 31) - 1;
        size_t offset = 0;
        size_t nWrite = count;
        while (nWrite > 0)
        {
          const int count = static_cast<int>(std::min(nWrite, batchMax));
          assert(count > 0);
          MPI_File_write(
              fh,
              (void*)(reinterpret_cast<const char*>(data) + offset),
              count, 
              MPI_BYTE, 
              &status);
          MPI_Get_count(&status, MPI_BYTE, &checkCount);
          if (count != checkCount)
            throw Exception(errString);
          nWrite -= count;
          offset += count;
        }
      }
  };

  /*#####################################*/
  /*#####################################*/
  /*#####################################*/

  /*********** Data types *************/

  class DataTypeBase
  {
    private:
      std::string name;
    public:
      DataTypeBase(const std::string _name) : name(_name) {};
      const std::string& getName() const {return name;}

      virtual void   resize(const size_t)   = 0;
      virtual size_t getElementSize() const = 0;
      virtual size_t getNumElements() const = 0;
      virtual size_t getNumBytes   () const = 0;

      virtual const void*  getDataPtr() const = 0;
      virtual       void*  getDataPtr()       = 0;
  };

  template<typename T>
    class DataType : public DataTypeBase
  {
    private:
      size_t numElements;
      T *data;
      void free()
      {
        if (data != NULL)
          ::free(data);
        data = NULL;
      }
      void malloc(const size_t _numElements)
      {
        numElements = _numElements;
        data = (T*)::malloc(sizeof(T)*numElements);
      }
    public:
      DataType(std::string name, const size_t n = 0) : DataTypeBase(name), data(NULL)
      {
        if (n > 0)
          malloc(n);
      }
      ~DataType() { free(); }
      void   resize(const size_t n) { assert(n > 0); free(); malloc(n); }
      size_t getElementSize() const {return sizeof(T);}
      size_t getNumElements() const {return numElements;}
      size_t getNumBytes   () const {return numElements*sizeof(T);}

      const T& operator[](const size_t i) const { return data[i]; }
            T& operator[](const size_t i)       { return data[i]; }

      const void* getDataPtr() const { return reinterpret_cast<const void*>(&data[0]); }
            void* getDataPtr()       { return reinterpret_cast<      void*>(&data[0]); }
  };

  /*#####################################*/
  /*#####################################*/
  /*#####################################*/

  /*********** Header ******************/

  class Header
  {
    private:
      struct DataInfo
      {
        enum {NAMEMAX = 255};
        char name[NAMEMAX+1];
        size_t elementSize;
        long_t offset;
        int nRank;
      };

      std::vector<DataInfo> data;
      double time;


    public:
      void   setTime(const double _time) { time = _time;}
      double getTime() { return time; }

      void printFields() const
      {
        const int n = data.size();
        for (int i = 0; i < n; i++)
        {
          fprintf(stderr, "---------Field %d ----------------------\n", i);
          fprintf(stderr, "Name=        %s\n", data[i].name);
          fprintf(stderr, "elementSize= %d\n", (int)data[i].elementSize);
          fprintf(stderr, "nRank=       %d\n", (int)data[i].nRank);
          fprintf(stderr, "offset=      %lld\n", data[i].offset);
        }
        fprintf(stderr, "------------------------------------\n");
      }

      Header() {}

      int find(const std::string &name) const
      {
        const int n = data.size();
        for (int i = 0; i < n; i++)
          if (std::string(data[i].name) == name)
            return i;
        return -1;
      }

      long_t getOffset(const int idx) const
      {
        assert(idx >= 0);
        return data[idx].offset;
      }
      size_t getElementSize(const int idx) const
      {
        assert(idx >= 0);
        return data[idx].elementSize;
      }
      int getNRank(const int idx) const
      {
        assert(idx >= 0);
        return data[idx].nRank;
      }
      long_t getDataOffset(const int idx) const
      {
        assert(idx >= 0);
        return data[idx].offset;
      }

      bool add(const DataTypeBase &dataVec, const long_t offset, const int nRank)
      {
        assert(dataVec.getName().size() <= DataInfo::NAMEMAX);
        if (find(dataVec.getName()) == -1)
        {
          DataInfo d;
          sprintf(d.name, dataVec.getName().c_str());
          d.offset      = offset;
          d.elementSize = dataVec.getElementSize();
          d.nRank       = nRank;
          data.push_back(d);
          return true;
        }
        return false;
      }

      void write(FileIO &fh)
      {
        assert(fh.isWrite());
        char versionString[16] = "V1";
        fh.write(versionString, 16*sizeof(char), "Error writing versionString.");
        fh.write(&time, sizeof(float), "Error writing time.");
        int nData = data.size();
        fh.write(&nData, sizeof(int), "Error writing nData.");
        fh.write(&data[0], sizeof(DataInfo)*nData, "Error writing dataInfo.");
      }

      void read(FileIO &fh)
      {
        assert(fh.isRead());
        char versionString[16];
        fh.read(versionString, 16*sizeof(char), "Error reading versionString.");
        assert(std::string(versionString) == "V1");

        fh.read(&time, sizeof(float), "Error reading time.");
        int nData;
        fh.read(&nData, sizeof(int), "Error reading nData.");

        data.resize(nData);
        fh.read(&data[0], sizeof(DataInfo)*nData, "Error reading dataInfo.");
      }
  };



  /*#####################################*/
  /*#####################################*/
  /*#####################################*/

  /*********** Core reader/writer *************/

  class Core
  {
    private:
      const int myRank;
      const int nRank;
      const MPI_Comm &comm;

      long_t dataOffsetGlb;
      FileIO *fhPtr;
      FileIO &fh;

      size_t numBytes;
      double dtIO;

      Header header;
      bool isMaster() const { return myRank == 0; }


    public:
      Header const & getHeader() { return header; }

    public:
      Core(const int _myRank,
          const int _nRank,
          const MPI_Comm &_comm,
          const IOTYPE iotype,
          const std::string fileName) :
        myRank(_myRank),
        nRank(_nRank),
        comm(_comm),
        dataOffsetGlb(sizeof(long_t)),
        fhPtr(new MPIFileIO(comm)),
        fh(*fhPtr),
        numBytes(0), dtIO(0)
      {
        switch (iotype)
        {
          case WRITE:
            fh.open(fileName, WRITE);
            break;
          case READ:
            fh.open(fileName, READ);
            long_t headerOffset;
            fh.read(&headerOffset, sizeof(long_t), "Unable to read header offset.");
            fh.seek(headerOffset);
            header.read(fh);
            break;
          default:
            assert(0);
        }
      }
      ~Core() { delete fhPtr; }

      bool read(DataTypeBase &data, const bool restart = false, const int reduceFactor = 1)
      {
        const double tRead = MPI_Wtime();
        /* make sure we are in the reading phase */
        assert(fh.isRead());

        /* find data set */
        const int idx = header.find(data.getName());
        if (idx == -1)
          return false;

        /* confirm that element size has the same length */
        if (header.getElementSize(idx) != data.getElementSize())
          return false;

        /* how many ranks wrote this file */
        const int nRankFile = header.getNRank(idx);

        /* read the number of elements written per rank */
        std::vector<long_t> numElementsPerRank(nRankFile);

        long_t offset = header.getDataOffset(idx);
        fh.seek(offset);
        fh.read(&numElementsPerRank[0], sizeof(long_t)*nRankFile, "Error while reading numElementsPerRank.");
        offset   += nRankFile*sizeof(long_t);
        numBytes += nRankFile*sizeof(long_t);

        std::vector<long_t> beg(nRankFile+1, 0), end(nRankFile+1, 0);
        for (int i = 0; i < nRankFile; i++)
        {
          end[i  ] = beg[i] + numElementsPerRank[i];
          beg[i+1] = end[i];
        }
        const long_t numElementsGlb = end[nRankFile-1];

        if (!(restart && nRankFile == nRank))
        {
          const long_t _nElementsPerRank = (numElementsGlb - 1 + nRank) / nRank;
          const long_t _beg = myRank * _nElementsPerRank;
          const long_t _end = std::min(_beg + _nElementsPerRank, numElementsGlb);
          beg[myRank] = _beg;
          end[myRank] = _end;
        }

        const long_t numElementsLoc = end[myRank] - beg[myRank];
        { /* sanity check */
          long_t sumGlb, sumLoc = numElementsLoc;
          MPI_Allreduce(&sumLoc, &sumGlb, 1, MPI_LONGT, MPI_SUM, comm);
          assert(sumGlb == numElementsGlb);
        }

        offset += beg[myRank]*data.getElementSize();;
        fh.seek(offset);

        if (reduceFactor <= 1)
        {
          data.resize(numElementsLoc);
          const long_t nBytes = (end[myRank] - beg[myRank]) * data.getElementSize();
          fh.read(data.getDataPtr(), nBytes, "Error while reading data.");
          numBytes+= numElementsGlb*data.getElementSize();
        }
        else
        {
          assert(reduceFactor == 1);
          long_t numElements = 0;
          for (long_t i = 0; i < numElementsLoc; i += reduceFactor)
            numElements++;
          data.resize(numElements);
          long_t el = 0;
          const size_t size = data.getElementSize();
          for (long_t i = 0; i < numElementsLoc; i += reduceFactor, el += size)
          {
            fh.read(
                reinterpret_cast<char*>(data.getDataPtr()) + el,
                size, "Error while reading reduced data.");
          }
        }

        dtIO += MPI_Wtime() - tRead;
        return true;
      }

      bool write(const DataTypeBase &data)
      {
        const double tWrite = MPI_Wtime();

        /* make sure we are in the writing phase */
        assert(fh.isWrite());

        long_t numElementsLoc = data.getNumElements();

        /* gather numELementsLoc to all ranks */
        std::vector<long_t> numElementsPerRank(nRank);
        MPI_Allgather(
            &numElementsLoc,        1, MPI_LONGT,
            &numElementsPerRank[0], 1, MPI_LONGT,
            comm);

        std::vector<long_t> beg(nRank+1, 0), end(nRank+1, 0);
        for (int i = 0; i < nRank; i++)
        {
          end[i  ] = beg[i] + numElementsPerRank[i];
          beg[i+1] = end[i];
        }

        const size_t numElementsGlb = end[nRank-1];

        if (isMaster())
        {
          /* add data description to the header */ 
          if (!header.add(data, dataOffsetGlb, nRank))
            throw Exception("Data type is already added.");

          /* write descirption about #elements at each rank */
          /* this is handy for restart functionality to avoid domain decomposition */
          fh.seek(dataOffsetGlb);
          fh.write(&numElementsPerRank[0], sizeof(long_t)*nRank, "Error while writing numElementsPerRank.");
        }
        numBytes += nRank*sizeof(long_t);

        dataOffsetGlb += sizeof(long_t)*nRank;

        fh.seek(dataOffsetGlb + beg[myRank]*data.getElementSize());
        const long_t nBytes = (end[myRank] - beg[myRank]) * data.getElementSize();
        assert(nBytes > 0);
        fh.write(data.getDataPtr(), nBytes, "Error while writing data.");

        dataOffsetGlb += numElementsGlb*data.getElementSize();

        numBytes += numElementsGlb*data.getElementSize();
        dtIO += MPI_Wtime() - tWrite;
        return true;
      }


      void close()
      {
        if (isMaster() && fh.isWrite())
        {
          fh.seek(0);
          fh.write(&dataOffsetGlb, sizeof(long_t), "Error while writing headerOffset.");
          fh.seek(dataOffsetGlb);
          header.write(fh);
        }

        fh.close();
      }

      double computeBandwidth() const { return numBytes/dtIO; }
  };

}

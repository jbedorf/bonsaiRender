CXX = mpicxx
CC  = mpicc
LD  = mpicxx

OMPFLAGS  = -fopenmp 
#OMPFLAGS += -D_GLIBCXX_PARALLEL

OFLAGS = -O3 -g -Wall 
# OMPFLAGS=

CXXFLAGS =  -fPIC $(OFLAGS) $(OMPFLAGS) -std=c++11



LDFLAGS = 

SRCPATH = ./

SRC1 = cvt2grid.cpp
#SRC2 = readBonsai.cpp
#SRC3 = benchmark.cpp
OBJ1 = $(SRC1:%.cpp=%.o)
OBJ2 = $(SRC2:%.cpp=%.o)
OBJ3 = $(SRC3:%.cpp=%.o)

PROG1 = cvt2grid
#PROG2 = readBonsai
#PROG3 = benchmark

RM = /bin/rm

all:	  $(PROG1) $(PROG2) $(PROG3)


$(PROG1): $(OBJ1) 
	$(LD) $(LDFLAGS) $^ -o $@ $(OMPFLAGS)
$(PROG2): $(OBJ2) 
	$(LD) $(LDFLAGS) $^ -o $@ $(OMPFLAGS)
$(PROG3): $(OBJ3) 
	$(LD) $(LDFLAGS) $^ -o $@ $(OMPFLAGS)




%.o: $(SRCPATH)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	/bin/rm -rf *.o $(PROG1) $(PROG2) $(OBJ1) $(OBJ2)

$(OBJ1): BonsaiIO.h  

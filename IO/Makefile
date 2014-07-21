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

SRC1 = BonsaiIO.cpp cvt_tipsy2bonsai.cpp
SRC2 = BonsaiIO.cpp readBonsai.cpp
OBJ1 = $(SRC1:%.cpp=%.o)
OBJ2 = $(SRC2:%.cpp=%.o)

PROG1 = cvt_tipsy2bonsai
PROG2 = readBonsai

RM = /bin/rm

all:	  $(PROG1) $(PROG2)


$(PROG1): $(OBJ1) 
	$(LD) $(LDFLAGS) $^ -o $@ $(OMPFLAGS)

$(PROG2): $(OBJ2) 
	$(LD) $(LDFLAGS) $^ -o $@ $(OMPFLAGS)



%.o: $(SRCPATH)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	/bin/rm -rf *.o $(PROG1) $(PROG2) $(OBJ1) $(OBJ2)

$(OBJ1): BonsaiIO.h  read_tipsy.h
$(OBJ2): BonsaiIO.h

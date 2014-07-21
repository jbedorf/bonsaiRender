CXX = mpicxx
CC  = mpicc
LD  = mpicxx

OMPFLAGS  = -fopenmp
#OMPFLAGS += -D_GLIBCXX_PARALLEL

OFLAGS = -O3 -g -Wall 
# OMPFLAGS=

CXXFLAGS =  -fPIC $(OFLAGS) -Wstrict-aliasing=2 $(OMPFLAGS)



LDFLAGS = 

SRCPATH = ./

SRC = BonsaiIO.cpp main.cpp
OBJ = $(SRC:%.cpp=%.o)

PROG = read_tipsy

RM = /bin/rm

all:	  $(PROG)


$(PROG): $(OBJ) 
	$(LD) $(LDFLAGS) $^ -o $@ $(OMPFLAGS)


%.o: $(SRCPATH)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	/bin/rm -rf *.o $(PROG)  $(PROG1)

$(OBJ): BonsaiIO.h 

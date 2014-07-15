CXX = g++
CC  = gcc
LD  = g++
F90  = gfortran

CUDA_TK  = /usr/local/cuda
CXXFLAGS = -O3 -g -Wall -openmp
CFLAGS = -O3 -g -Wall

CXXFLAGS=-I ./vfx
LDFLAGS = -lGL -lGLU -lglut

OBJ = main.o renderloop.o anyoption.cpp

PROG= renderer

all: $(PROG)

$(PROG): $(OBJ)
	$(LD) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $< -o $@

clean:
	/bin/rm -rf $(OBJ) $(PROG)




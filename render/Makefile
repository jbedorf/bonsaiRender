CXX = mpic++
CC  = mpicc
LD  = mpic++
F90  = mpif90

CUDA_TK  = /usr/local/cuda
CXXFLAGS = -O3 -g -Wall -openmp -std=c++11
CFLAGS = -O3 -g -Wall

ENGINE=./particles
INCLUDE_FLAGS= -I./ -I$(ENGINE)/
LDFLAGS = -lGL -lGLU -lglut -lGLEW
#-L$(PWD)/lib/linux -lGLEW_x86_64

OBJ = main.o $(ENGINE)/renderloop.o anyoption.o $(ENGINE)/render_particles.o paramgl.o

PROG= renderer

all: $(PROG)

$(PROG): $(OBJ)
	$(LD) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $< -o $@

clean:
	/bin/rm -rf $(OBJ) $(PROG)


$(OBJ): anyoption.h  paramgl.h  param.h  read_tipsy.h  RendererData.h  $(ENGINE)/renderloop.h  $(ENGINE)/render_particles.h  tipsydefs.h  vector_math.h



CXX = mpic++
CC  = mpicc
LD  = mpic++
F90  = mpif90

CUDA_TK  = /usr/local/cuda
CUDA_TK = $(HOME)/usr/local/cuda-5.5
CXXFLAGS = -O3 -g -Wall -openmp -std=c++11
CXXFLAGS += -Wno-reorder
CXXFLAGS += -Wno-write-strings
CFLAGS = -O3 -g -Wall

NVCC=$(CUDATK)/bin/nvcc
NVCCFLAGS=-arch=sm_35

ENGINE=particles
ENGINE=smoke


INCLUDE_FLAGS= -I./ -I./$(ENGINE)/
INCLUDE_FLAGS += -I$(CUDA_TK)/include
LDFLAGS = -lGL -lGLU -lglut -lGLEW
#-L$(PWD)/lib/linux -lGLEW_x86_64

OBJ = main.o $(ENGINE)/renderloop.o anyoption.o $(ENGINE)/renderer.o paramgl.o
ifeq ($(ENGINE),smoke)
	OBJ += tr.o
	OBJ += smoke/GLSLProgram.o
	OBJ += smoke/SmokeShaders.cpp
	OBJ += smoke/framebufferObject.o
#	OBJ += smoke/depthSort.cu_o
	OBJ += smoke/depthSort.o
endif

PROG= renderer

all: $(PROG)

$(PROG): $(OBJ)
	$(LD) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE_FLAGS) -c $< -o $@

%.cu_o: %.cu
	$(NVCC) $(NVCCFLAGS) $(INCLUDE_FLAGS) -c $< -o $@

clean:
	/bin/rm -rf $(OBJ) $(PROG)

clean_bak:
	find . -name '*~' -exec /bin/rm -rf '{}' \;



$(OBJ): anyoption.h  paramgl.h  param.h  read_tipsy.h  RendererData.h  $(ENGINE)/renderloop.h  $(ENGINE)/renderer.h  tipsydefs.h  vector_math.h



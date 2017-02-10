CC=gcc
#For GNU compiler
ifeq ($(CC),gcc)
	CC=gcc 
	CXX=g++ 
	CFLAGS=  -O3 -g -fopenmp -std=gnu99 -Wall -ffast-math  -I/glb/home/uscwc6/include
	CXXFLAGS= -O3 -g -fopenmp -Wall -ffast-math  -I/glb/home/uscwc6/include
	LDFLAGS= -O3 -g -fopenmp -ffast-math  -L/glb/home/uscwc6/lib
	LDLIBS=-lfftw3 -lm -lfftw3_threads -lpthread -lstdc++ #-lfftw3f 
endif

#For Intel compiler
ifeq ($(CC),icc)
	CC=icc
	CXX=icpc
	CFLAGS=  -std=gnu99 -O3 -openmp -fp-model fast=1 -finline -ipo -I/glb/home/uscwc6/include#-no-prec-div
	CXXFLAGS= -O3  -fp-model fast=1 -finline -ipo  #-no-prec-div
	LDFLAGS=  -O3 -openmp  -fp-model fast=1 -finline -ipo -L/glb/home/uscwc6/lib#-no-prec-div -finline
	LDLIBS= -lm -lfftw3 -lfftw3_threads
endif

#For MA
ifeq ($(OPT),ma)
	CFLAGS= -std=gnu99 -O3 -openmp  -Wall -I/usr/include -xAVX -ipo -no-prec-div \
		   	-I/glb/data/comp_mod/Tools/fftw/include -fp-model fast=1 -finline #-fno-stack-limit 
	LDLIBS= -O3 -openmp -L/glb/data/comp_mod/Tools/fftw/lib -lfftw3 -lfftw3_threads -lm 
endif

#For MA
ifeq ($(OPT),ma)
	CFLAGS= -std=gnu99 -O3 -openmp -g -Wall -I/usr/include -xAVX -ipo -no-prec-div \
		   	-I/glb/data/comp_mod/Tools/fftw/include -fp-model fast=1 -finline #-fno-stack-limit 
	LDLIBS= -O3 -openmp -L/glb/data/comp_mod/Tools/fftw/lib -lfftw3 -lfftw3_threads -lm 
	CXXFLAGS= -O3  -fp-model fast=1 -finline -openmp -g -Wall -I/usr/include -xAVX -ipo -no-prec-div \
		   	-I/glb/data/comp_mod/Tools/fftw/include -fp-model fast=1 -finline #-fno-stack-limit 
	LDFLAGS= -O3 -openmp  -fp-model fast=1 -finline #-ipo -no-prec-div -finline
endif

#For PGI compiler
#ifeq ($(CC),pgcc)
#	PGI_CUDA_PATH =/opt/pgi/linux86-64/2013/cuda/4.2
#	FFTW_PATH =/glb/home/uscwc6
#	CC=pgcc
#	CXX=pgCC
#	CFLAGS= -I$(FFTW_PATH)/include -I$(PGI_CUDA_PATH)/include -c99 -O3 -mp -fast -Mipa=fast,inline -Minline  -acc -Minfo=accel -ta=nvidia:fastmath,keepgpu,time,cc35 -Mcuda
#	CXXFLAGS= -c99 -O3 -fast -Mipa=fast,inline -Minline
#	LDLIBS= -L$(FFTW_PATH)/lib -L$(PGI_CUDA_PATH)/lib64 -lfftw3 -lm -lrt -lfftw3_threads -lpthread -lstdc++ -mp -acc -lcudart -lcufft 
#endif

#For PGI compiler
ifeq ($(CC),pgcc)
	FFTW_PATH =/glb/home/uscwc6
	CC=pgcc
	CXX=pgCC
	CFLAGS= -I$(FFTW_PATH)/include -O3 -mp -fast -Mipa=fast,inline -Minline  #-acc -Minfo=accel -ta=nvidia:fastmath,keepgpu,time,cc35 -Mcuda
	CXXFLAGS= -c99 -O3 -fast -Mipa=fast,inline -Minline
	LDLIBS= -L$(FFTW_PATH)/lib  -lfftw3 -lm -lrt -lfftw3_threads -lpthread -lstdc++ -mp #-acc -lcudart -lcufft
endif


#For NVCC Compiler
CUDA_DIR=/opt/cuda/5.5
NVCC= $(CUDA_DIR)/bin/nvcc
NVCCFLAGS = -O3 -I$(CUDA_DIR)/include -arch=sm_35
NVCCLIBS = -L$(CUDA_DIR)/lib64 -lcufft -lcudart

ifeq ($(CC),nvcc)
	CC=$(NVCC)
	CFLAGS = $(NVCCFLAGS) -Xcompiler -fopenmp --use_fast_math --compiler-options "-std=gnu99"
	CXX = $(NVCC)
	CXXFLAGS = -O3 --use_fast_math 
	LDLIBS = $(NVCCLIBS) -lgomp -lfftw3_threads -lfftw3 
endif

C_SRCS       = $(wildcard *.c)
C_OBJS       = $(C_SRCS:.c=.c.o)

CXX_SRCS     = $(wildcard *.cc)
CXX_OBJS     = $(CXX_SRCS:.cc=.cc.o)

CU_SRCS		 = $(wildcard *.cu)
CU_OBJS		 = $(CU_SRCS:.cu=.cu.o)

SRCS         = $(C_SRCS) $(CXX_SRCS) $(CU_SRCS)
OBJS         = $(C_OBJS) $(CXX_OBJS) $(CU_OBJS)

.PHONY: all	clean


.SUFFIXES: .c .cc .cu .o

%.c.o: %.c
	 $(CC) $(CFLAGS) -c $< -o $@

%.cu.o: %.cu
	 $(NVCC) $(NVCCFLAGS) -c $< -o $@

%.cc.o: %.cc
	 $(CXX) $(CXXFLAGS) -c $< -o $@

experiment: $(OBJS)
	 $(CC) -o $@ $(OBJS) $(NVCCLIBS) $(LDFLAGS) $(LDLIBS) 

all:	experiment


clean:
	$(RM) $(OBJS) *.gpu ./experiment

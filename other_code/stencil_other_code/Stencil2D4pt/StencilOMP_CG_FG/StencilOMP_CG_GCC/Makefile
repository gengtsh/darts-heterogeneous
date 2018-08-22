CC=gcc
CXX=g++
#SHELL:=/bin/bash

OMPFLAGS=-fopenmp 
MATHFLAGS=-fexcess-precision=standard -freciprocal-math 
#CFLAGS=-std=c11 -D_XOPEN_SOURCE -D_POSIX_C_SOURCE=200809L -Wall -Wextra -pedantic -Werror -O3 $(OMPFLAGS) 
OPTIM=-O3
DEBUG=#-ggdb
CFLAGS=-std=c11 -march=native -D_XOPEN_SOURCE -D_POSIX_C_SOURCE=200809L -Wall -Wextra -Wno-unknown-pragmas -pedantic $(DEBUG) $(OPTIM) 
CXXFLAGS=-std=c++11 -D_XOPEN_SOURCE -D_POSIX_C_SOURCE=200809L -Wall -Werror $(MATHFLAGS)
CFLAGSOMP:=$(CFLAGS) $(OMPFLAGS)
CXXFLAGSOMP:=$(CXXFLAGS) $(OMPFLAGS)
LDFLAGS=-lm -lrt -lpthread
MODULE_LOAD=. ${MODULESHOME}/init/bash && module load

.PHONY: all 

#all: OMPLaplaceSolver 
#all:  SeqLaplaceSolver
all: OMPLaplaceSolver SeqLaplaceSolver
#all: LoadIntelLib OMPLaplaceSolver SeqLaplaceSolver
.SECONDEXPANSION:
LoadIntelLib:
	${MODULE_LOAD} intel_libomp_oss; \
	echo "LIBIOMP5_PATH= " $${LIBIOMP5_PATH}


INCS_INTEL=-I${LIBIOMP5_PATH}lin_32e/include -I${LIBIOMP5_PATH}lin_32e/include_compat -I${LIBIOMP5_PATH}common/include -I${LIBIOMP5_PATH}common/include_compat

LIB_INTEL=-L${LIBIOMP5_PATH}lin_32e/lib
LD_INTEL=-liomp5
INTEL=
OMPLINKINTEL=$(LIB_INTEL) $(INCS_INTEL) $(LD_INTEL)

Echo:
	echo ${INCS_INTEL}
	echo ${OMPLINKINTEL}

ifdef INTEL
OMPLINK=$(OMPLINKINTEL)
else
OMPLINK= $(OMPFLAGS)
endif

OMPLaplaceSolver: stencil_omp.o main_omp.o timer_omp.o
	$(CC) $(OMPLINK) -o $@ $^ $(LDFLAGS)

SeqLaplaceSolver: stencil_seq.o main_seq.o timer_seq.o
	$(CC) -o $@ $^ $(LDFLAGS)


timer_omp.o: timer.c
	$(CC) $(CFLAGSOMP) -c -o $@ $<

timer_seq.o: timer.c
	$(CC) 			  $(CFLAGS) -c -o $@ $<

main_omp.o: main.c
	$(CC)  $(CFLAGSOMP) -c -o $@ $<

main_seq.o: main.c
	$(CC)  $(CFLAGS) -c -o $@ $<

stencil_seq.o:stencil_omp.c 
	$(CC) $(CFLAGS) -c -o $@ $<

stencil_omp.o:stencil_omp.c
	$(CC) $(CFLAGSOMP) -c -o $@ $<


%_omp.o: %.c %.h 
		$(CC) $(CFLAGSOMP) -c -o $@ $<

clean:
	rm -f *.o OMPLaplaceSolver SeqLaplaceSolver

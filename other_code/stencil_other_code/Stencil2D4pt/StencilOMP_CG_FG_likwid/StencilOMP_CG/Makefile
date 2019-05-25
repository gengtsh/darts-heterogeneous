#CC=gcc
#CXX=g++

CC=/software/apps/gcc/8.1/bin/gcc
CXX=/software/apps/gcc/8.1/bin/g++
LIB=-L/software/apps/gcc/8.1/lib64 -L/software/apps/likwid/4.3.2/lib
INC=-I/software/apps/gcc/8.1/include -I/software/apps/likwid/4.3.2/include

DLIKWIDFLAG=-DLIKWID_PERFMON

LIKWIDFLAGS= $(LIB) $(INC) $(DLIKWIDFLAG)

OMPFLAGS=-fopenmp
MATHFLAGS=-fexcess-precision=standard -freciprocal-math 
OPTIM=-O3
DEBUG=#-ggdb
CFLAGS=-std=c11 -march=native -D_XOPEN_SOURCE -D_POSIX_C_SOURCE=200809L -Wall -Wextra -Wno-unknown-pragmas -pedantic $(DEBUG) $(OPTIM) $(LIKWIDFLAGS) 
CXXFLAGS=-std=c++11 -D_XOPEN_SOURCE -D_POSIX_C_SOURCE=200809L -Wall -Werror $(MATHFLAGS) $(LIKWIDFLAGS)

LDFLAGS=-lm -lrt -lpthread $(LIB) $(DLIKWIDFLAG)  -llikwid


all: OMPLaplaceSolver SeqLaplaceSolver

OMPLaplaceSolver: stencil_omp.o main_omp.o timer_omp.o
	$(CC) $(OMPFLAGS) -o $@ $^ $(LDFLAGS) 

SeqLaplaceSolver: stencil_seq.o main_seq.o timer_seq.o
	$(CC) -o $@ $^ $(LDFLAGS)


timer_omp.o: timer.c
	$(CC) $(OMPFLAGS) $(CFLAGS) -c -o $@ $<

timer_seq.o: timer.c
	$(CC) 			  $(CFLAGS) -c -o $@ $<

main_omp.o: main.c
	$(CC) $(OMPFLAGS) $(CFLAGS) -c -o $@ $<

main_seq.o: main.c
	$(CC) 			  $(CFLAGS) -c -o $@ $<

stencil_seq.o:stencil_omp.c 
	$(CC) 			  $(CFLAGS) -c -o $@ $<

stencil_omp.o:stencil_omp.c
	$(CC) $(OMPFLAGS) $(CFLAGS) -c -o $@ $<


%_omp.o: %.c %.h 
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o OMPLaplaceSolver SeqLaplaceSolver

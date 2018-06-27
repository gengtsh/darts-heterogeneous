CC=gcc
CXX=g++
OMPFLAGS=-fopenmp
MATHFLAGS=-fexcess-precision=standard -freciprocal-math 
#CFLAGS=-std=c11 -D_XOPEN_SOURCE -D_POSIX_C_SOURCE=200809L -Wall -Wextra -pedantic -Werror -O3 $(OMPFLAGS) 
OPTIM=-O3
DEBUG=#-ggdb
CFLAGS=-std=c11 -march=native -D_XOPEN_SOURCE -D_POSIX_C_SOURCE=200809L -Wall -Wextra -Wno-unknown-pragmas -pedantic $(DEBUG) $(OPTIM) 
CXXFLAGS=-std=c++11 -D_XOPEN_SOURCE -D_POSIX_C_SOURCE=200809L -Wall -Werror $(MATHFLAGS)
LDFLAGS=-lm -lrt -lpthread

#all: OMPLaplaceSolver SeqLaplaceSolver
all: OMP_FGV2 

OMP_FGV2: stencil_omp.o main_omp.o timer_omp.o
	$(CC) $(OMPFLAGS) -o $@ $^ $(LDFLAGS)

#SeqLaplaceSolver: stencil_seq.o main_seq.o timer_seq.o
#	$(CC) -o $@ $^ $(LDFLAGS)


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
	$(CC) $(OMPFLAGS) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o OMP_FGV2 SeqLaplaceSolver

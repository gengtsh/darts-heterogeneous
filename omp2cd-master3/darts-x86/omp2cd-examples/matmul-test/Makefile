include ../Makefile.def

SHELL=/bin/bash

TARGET=matmul.exe
SOURCES=matmul.c

.PHONY: all darts omp clean-darts clean

all: omp 

darts:
	$(OMP2CD) $(OMP2CDFLAGS) $(SOURCES) -- $(CLANGFLAGS)
	@clang-format -i -style=Webkit ./*.output.darts.*

omp: $(SOURCES:.c=.o)
	$(CC) $(SOURCES:.c=.o) $(LDFLAGS_OPENMP) -o $(TARGET)

%.o: %.c
	$(CC) -c $< $(CFLAGS_OPENMP) -o $@

clean-darts:
	rm -f ./*.output.darts.*
	
clean:
	rm -f *~ *.o *.obj $(TARGET)
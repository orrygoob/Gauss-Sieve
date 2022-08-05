CC=g++
OPTFLAGS= -march=native -O3 -pipe
STRICT= -std=c++17 -pedantic -Wall -fno-builtin
CFLAGS   = $(STRICT) $(OPTFLAGS)

UNAME_S = $(shell uname -s)
ifeq ($(UNAME_S),Darwin) # if OSX use brew paths
	LIBS = -I/usr/local/opt/OpenBLAS/include/ -L/usr/local/opt/OpenBLAS/lib -lopenblas
else # Otherwise assume on path i.e. windows
	LIBS = -lopenblas
endif

MAIN  = main.cc
OUTPUT = main

all: $(MAIN)
	$(CC) $(MAIN) -o $(OUTPUT) $(LIBS) $(CFLAGS)

clean:
	rm main
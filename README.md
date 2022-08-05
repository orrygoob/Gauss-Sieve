## CONTENTS OF THIS FILE
---------------------

 * Introduction
 * Prerequisites
 * Installation
 * Usage
 * Performance Notes
 * To Do

## INTRODUCTION
------------

Author: Orry Gooberman

This project is a C++ implementation of the Gaussian Sieve originally proposed by Daniele Micciancio and Panagiotis Voulgaris in their paper "Faster Exponential Time Algorithms for the Shortest Vector Problem". It finds an approximation for the shortest vector in a given integer lattice. This implementation is a heavily modified version of Voulgaris' 2011 algorithm and has been updated to use current linear algebra libraries (CBLAS) instead of the depreciated NTL 5.0 (As a result program no longer supports arbitaty length integers and using double precision floats for most calculations). The aim of this simple single-threaded algorithm is to provide a simple and readable example that can be built on top of.


## PREREQUISITES
------------

OpenBLAS
- On Mac OS this can be installed using `brew install openblas`
- For windows this can be installed using binaries from https://github.com/xianyi/OpenBLAS. Ensure the installation location is on the systems path so the compiler can find it.

g++ Compiler (Bundled with the GCC compiler)


## INSTALLATION
------------

1. Download and unzip the repository
2. Open the folder using your choice of command line (i.e. '`cd /PATH_TO_FOLDER`')
3. Run '`make`' to generate the executable file '`./main`'

If the make command fails for some reason you can use the following command to manually build the project:

```shell
  g++ main.cc -o main -I/PATH_TO_BLAS/include/ -L/PATH_TO_BLAS/lib -lopenblas -std=c++17 -pedantic -Wall -fno-builtin -march=native -O3 -pipe
```

Note that the path to the openblas library will need to be replaced with your systems path.


## USAGE
----------------

To run the program use the following command from inside the project folder where basis.txt is an LLL reduced lattice basis in the standard fplll format:

```shell
 ./main -f challenges/basis.txt
```
There are further optional flags that can be set to the programs resources and parameters.

```shell
-v # Enter verbose mode
```
```shell
-l INTEGER # Specify the limit on the number of vectors that can be stored at any one time (DEFAULT 100000)
```
```shell
-s INTEGER # Specify the limit on the number of vectors in the stack waiting to be reduced (DEFAULT 1000)
```
```shell
-h # Display help on the available parameters
```

For example the following command will run the program on "challenges/svp50LLL" in verbose mode with a database made of up to 10000 vectors and a maximum of 100 vectors in the stack at a time.
```shell
 ./main -f challenges/svp50LLL.txt -l 10000 -s 100 -v
```

## PERFORMANCE NOTES
-------------

Analysis has been done on how this algorithm performs with different dimensions, available memory, and collisions before execution is stopped.


## TO DO
-------------

- [ ] Add in plots of performance and vector quality along with descriptions
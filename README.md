## CONTENTS OF THIS FILE
---------------------

 * Introduction
 * Prerequisites
 * Installation
 * Usage
 * Branches
 * Performance Notes
 * To Do

## INTRODUCTION
------------

Author: Orry Gooberman

This project is a C++ implementation of the Gaussian Sieve originally proposed by Daniele Micciancio and Panagiotis Voulgaris in their paper "Faster Exponential Time Algorithms for the Shortest Vector Problem". It finds an approximation for the shortest vector in a given integer lattice. This implementation is a heavily modified version of Voulgaris' 2011 algorithm and has been updated to use current linear algebra libraries (CBLAS) instead of the depreciated NTL 5.0 (As a result program no longer supports arbitary length integers and uses double precision floats for most calculations). The aim of this simple single-threaded algorithm is to provide a simple and readable example that can be built on top of.


## PREREQUISITES
------------

OpenBLAS
- On Mac OS this can be installed using `brew install openblas`
- For windows this can be installed using binaries from https://github.com/xianyi/OpenBLAS. Ensure the installation location is on the systems path so the compiler can find it.

g++ Compiler (Bundled with the GCC compiler)


## INSTALLATION
------------

1. Download and unzip the repository
2. Open the specific sieve folder using your choice of command line (i.e. '`cd /PATH_TO_FOLDER`')
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
 ./main -f ../challenges/basis.txt
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

For example the following command will run the program on "../challenges/svp50LLL" in verbose mode with a database made of up to 10000 vectors and a maximum of 100 vectors in the stack at a time.
```shell
 ./main -f ../challenges/svp50LLL.txt -l 10000 -s 100 -v
```

## BRANCHES
-------------

In this repository there are different folders for slightly modified Gauss Sieves. They are listed below with descriptions.

- simple gauss sieve - The stable standard branch. Implements the Gauss Sieve exactly as Micciancio and Voulgaris' paper describes.
- combined vector db reduction - Combines the separate steps of reducing the vector and reducing the database into a single for loop. First the algorithm tries to reduce $v$ using the database vector $w$, if and only if $v$ has not been reduced it is then used to attempt to reduce $w$. Leads to reductions in the database occurring as soon as possible and results in less inner products being calculated.

## PERFORMANCE NOTES
-------------

Analysis has been done on how this algorithm performs with different dimensions, available memory, and collisions before execution is stopped. The data for this analysis was all taken from a the simple Gauss Sieve implementation.

### Short Vector Quality
![Lattice Challenge Comparison](graphs/norm%20of%20only%20seed%200.png)

This compares the results of this SVP solver with the historical hall of fame records from latticechallenge.org. Proves that this solver is at least as good as all previous records for seed 0. (Records for other seeds were discounted as users have searched for better lattices with shorter vectors by using different seeds).

-------------
### Run Time
![Logarithmic Run Time](graphs/runTimeLogarithmic.png)

This is a plot of the logarithmic run time vs lattice dimension and clearly shows that the run time of the program is exponentially tied to the dimension of the lattice.

-------------
### Collision Limit
![Collision Limit vs Short Vector Quality](graphs/Linear%20Normalised%20Norm%20Collision%20Limit%20Heatmap.png)
![Collision Limit vs Run Time](graphs/Linear%20Run%20Time%20Collision%20Count%20Heatmap.png)

As you can see, changing the collision limit does not affect the run time of the program, only the quality of the output vector.

-------------
### Maximum Number of Vectors
![Vector Limit vs Short Vector Quality](graphs/Linear%20Norm%20Vector%20Limit%20Heatmap.png)
![Vector Limit vs Run Time](graphs/Log%20Run%20Time%20Vector%20Limit%20Heatmap.png)

As you can see, changing the collision limit only affects the run time of the program (especially when very small) and not the quality of the output vector. Additionally you can also see that there is no apparent downside to allowing many more vectors to be stored than strictly necessary.

## TO DO
-------------

- [ ] Add implementation of sieve with 64 bit simhash as implemented by Ducas
- [ ] Add implementation of hash sieve using angular locality sensitive hashing as described by Laarhoven
- [ ] Add implementation of progressive sieving as described by Laarhoven
- [ ] Add implementation of sieving by averages

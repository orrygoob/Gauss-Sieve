#include <getopt.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <chrono>

#include "math.h"
#include "sampler.h"
#include "cblas.h"

// Based off of the original 2011 implementation of gsieve by Panagiotis Voulgaris
// Original can be found at https://cseweb.ucsd.edu//~pvoulgar/impl.html

// Use fplll to LLL reduce basis before running into program. (Or use prereduced ones from svpchallenge.org located in the challenges subfolder)
// Use command 'make' to compile to compile code into './main'. Ensure that gcc and openblas are installed and in your system path.
// Run using: './main -f challenges/svp40LLL.txt'

const char *USAGE =
    "Usage: %s [OPTION]\n"
    "-f: FILE_NAME:\t\tA file containing the input basis.\n"
    "-l: VECTOR LIMIT\tSpecify the maximum size of the vector database.\n"
    "-s: STACK LIMIT\tSpecify the maximum size of the vector stack.\n"
    "-v:\t\t\tVerbose mode, print additional information.\n"
    "-h:\t\t\tPrint this help.\n";

// Program Parameters and their defaults
int dimension = 0;
long arraySize = 100000;
int stackSize = 1000;
char *inputFileName = nullptr;
bool verbose = false;

// db pointers (set index to 1 to include 0 vector)
double *dbPreimage = nullptr;
double *db = nullptr;
double *dbSqrNorm = nullptr;
int indexL = 1;

// stack pointers
double *stackPreimage = nullptr;
double *stack = nullptr;
int indexS = 0;

// Count number of lines in file and use this to get dimension of lattice
int getDimension()
{
  std::ifstream inputFile(inputFileName);
  if (!inputFile.good())
    return 0;

  // Get dimensions based on number of lines in file
  int numberOfLines = 0;
  std::string line;
  while (std::getline(inputFile, line))
    ++numberOfLines;

  inputFile.close();
  return numberOfLines - 1;
}

// Read in file, return dimension of basis
bool readInFile(double *B)
{
  std::ifstream inputFile(inputFileName);
  if (!inputFile.good())
    return false;

  // Get file contents
  int i = 0;
  std::string basisString;
  while (std::getline(inputFile, basisString))
  {
    std::stringstream ss(basisString);

    double number;
    char c;

    // Skip open brackets on line to the first number
    while (ss >> c && c == '[')
      ;
    ss.putback(c);

    // Iterate through each number on line
    while (ss >> number)
    {
      B[i] = number;
      i++;
    }
  }
  inputFile.close();
  return true;
}

bool stackFullErrorFlag = false; // Keep track of stack full error so only reported if new occurrence
void addToStack(double *srcPreimage, double *srcV)
{
  // Add to stack
  if (indexS < stackSize)
  {
    cblas_dcopy(dimension, srcV, 1, &stack[indexS * dimension], 1);
    cblas_dcopy(dimension, srcPreimage, 1, &stackPreimage[indexS * dimension], 1);
    indexS++;
  }
  else
  {
    if (!stackFullErrorFlag)
      std::cerr << "Stack Full" << std::endl;
    stackFullErrorFlag = true;
  }
}

void popFromStack(double *destPreimage, double *destV)
{
  if (indexS > 0)
  {
    indexS--;
    cblas_dcopy(dimension, &stack[indexS * dimension], 1, destV, 1);
    cblas_dcopy(dimension, &stackPreimage[indexS * dimension], 1, destPreimage, 1);
    stackFullErrorFlag = false;
  }
  else
  {
    std::cerr << "Stack Empty" << std::endl;
  }
}

void addToDB(double *srcPreimage, double *srcV, double sqrNorm)
{
  // Add to database
  int dest = -1;
  if (indexL < arraySize)
  {
    dest = indexL;
    indexL++;
  }
  else
  {
    // If array has run out of space purge largest vector
    // std::cout << "Vector database full, removing largest vector." << std::endl;
    int max = -1;
    for (int i = 0; i < indexL; i++)
    {
      if (dbSqrNorm[i] > max)
      {
        max = dbSqrNorm[i];
        dest = i;
      }
    }
  }

  // Copy into designated index
  cblas_dcopy(dimension, srcPreimage, 1, &dbPreimage[dest * dimension], 1);
  cblas_dcopy(dimension, srcV, 1, &db[dest * dimension], 1);
  dbSqrNorm[dest] = sqrNorm;
}

void removeFromDB(int i)
{
  // Remove element in index i from db (Replace with last in array if not last element)
  indexL--;
  cblas_dcopy(dimension, &db[indexL * dimension], 1, &db[i * dimension], 1);
  cblas_dcopy(dimension, &dbPreimage[indexL * dimension], 1, &dbPreimage[i * dimension], 1);
  dbSqrNorm[i] = dbSqrNorm[indexL];
}

bool reduce(double *vPreimage, double *v, double *vNorm, double *wPreimage, double *w, double *wNorm, double vw)
{
  // Attempt to reduce v using w
  long vwAbsDouble = 2 * (vw > 0 ? vw : -vw);

  // Reduce v with w if possible
  if (vwAbsDouble > *wNorm)
  {
    // v -= q * db[i] where q = number of times projected vector fits into db[i]
    long q = round((double)vw / *wNorm);
    cblas_daxpy(dimension, -q, w, 1, v, 1);
    cblas_daxpy(dimension, -q, wPreimage, 1, vPreimage, 1);
    *vNorm += q * q * (*wNorm) - 2 * q * vw;
    return true;
  }
  return false;
}

int main(int argc, char **argv)
{
  // Get command line arguments
  if (verbose)
    std::cout << "Getting command line arguments" << std::endl;
  int option;
  while ((option = getopt(argc, argv, "vf:l:s:h?")) != -1)
  {
    switch (option)
    {
    case 'v':
      verbose = true;
      break;
    case 'f':
      inputFileName = optarg;
      break;
    case 'l':
      arraySize = atoi(optarg);
      break;
    case 's':
      stackSize = atoi(optarg);
      break;
    case 'h':
    case '?':
      fprintf(stderr, USAGE, argv[0]);
      return -1;
    }
  }

  // --------- READ IN BASIS FILE --------
  if (inputFileName == nullptr)
  {
    std::cerr << "Must specify lattice basis file name by using -f option" << std::endl;
    return 1;
  }

  if (verbose)
    std::cout << "Reading in lattice file" << std::endl;

  // Get dimension and allocate memory for basis
  dimension = getDimension();
  if (dimension < 2)
  {
    std::cerr << "Error getting dimension of file" << std::endl;
    return 1;
  }

  // Allocate memory for basis
  double *B = (double *)malloc(sizeof(double) * dimension * dimension);
  // Read in basis into matrix
  if (!readInFile(B))
  {
    std::cerr << "Error reading in file" << std::endl;
    return 1;
  }

  // ---------------------------------------- INITIALISE PROGRAM --------------------------------------

  if (verbose)
  {
    std::cout << "Initialising program" << std::endl;
    std::cout << "Configuring sampler" << std::endl;
  }
  // Set seed for random number generator to current time and init sampler
  std::srand(time(0));
  KleinSampler sampler;
  sampler.Init(B, dimension, dimension);

  if (verbose)
    std::cout << "Allocating database memory" << std::endl;

  // Vector database array - Index of next free space is 1 as arrays already contain 0 vector
  // dbPreimage = preimage, db = lattice point, dbSqrNorm = Squared Norm of lattice point
  dbPreimage = (double *)calloc(dimension * arraySize, sizeof(double));
  db = (double *)calloc(dimension * arraySize, sizeof(double));
  dbSqrNorm = (double *)calloc(arraySize, sizeof(double));

  if (verbose)
    std::cout << "Allocating stack memory" << std::endl;

  // Stack of vectors that have changed so are waiting to be reduced again. Start with all basis in stack
  stackPreimage = (double *)calloc(dimension * stackSize, sizeof(double));
  stack = (double *)calloc(dimension * stackSize, sizeof(double));

  indexS = dimension;
  memcpy(stack, B, dimension * dimension * sizeof(double));
  for (int i = 0; i < dimension; i++)
    for (int j = 0; j < dimension; j++)
      stackPreimage[i * dimension + j] = (j == i) ? 1 : 0;

  if (verbose)
    std::cout << "Configuring other parameters" << std::endl;

  long long collisionLimit = 0.1 * indexL + 200;
  long long collisionCount = 0;

  double *vPreimage = (double *)calloc(dimension, sizeof(double));
  double *v = (double *)calloc(dimension, sizeof(double));
  double vNorm = 0;
  // ------------------------- START SIEVING HERE ---------------------------
  if (verbose)
    std::cout << "START SIEVING" << std::endl;

  auto start = std::chrono::high_resolution_clock::now();

  while (collisionCount < collisionLimit)
  {
    if (indexS > 0)
    {
      // If a vector is waiting to be reduced pop it from stack
      popFromStack(vPreimage, v);
    }
    else
    {
      // Otherwise randomly sample a vector (returns preimage of lattice point)
      sampler.Sample(vPreimage);

      // Convert from preimage to real coordinates (v = B * vPreimage)
      cblas_dgemv(CblasColMajor, CblasNoTrans, dimension, dimension, 1, B, dimension, vPreimage, 1, 0, v, 1);
    }

    // Get squared norm of vector
    vNorm = cblas_ddot(dimension, v, 1, v, 1);

    // Reduce v
    bool reduced = false;
    for (int i = 0; i < indexL; i++)
    {
      double *w = &db[i * dimension];
      double *wPreimage = &dbPreimage[i * dimension];
      double *wNorm = &dbSqrNorm[i];
      double vw = cblas_ddot(dimension, v, 1, w, 1);

      // Attempt to reduce v using w, recording if a reduction took place
      reduced |= reduce(vPreimage, v, &vNorm, wPreimage, w, wNorm, vw);
    }

    // Check for a collision where v was not independent of vectors in database
    if (vNorm == 0)
    {
      collisionCount++;
      continue;
    }

    // If vector was reduced then add to stack of vectors to be reduced on next round
    if (reduced)
    {
      addToStack(vPreimage, v);
      continue;
    }
      
    // v is now fully reduced so loop through db and reduce all vectors using v
    for (int i = 0; i < indexL; i++)
    {
      double *w = &db[i * dimension];
      double *wPreimage = &dbPreimage[i * dimension];
      double *wNorm = &dbSqrNorm[i];
      double vw = cblas_ddot(dimension, v, 1, w, 1);

      if (reduce(wPreimage, w, wNorm, vPreimage, v, &vNorm, vw))
      {
        // Reduce w with v and add to stack if successful
        // Because this reduction only takes place if w wasn't used to reduce v we know that vectors are independent
        addToStack(wPreimage, w);
        removeFromDB(i);
        i--;
      }
    }

    addToDB(vPreimage, v, vNorm);
    
    // Update collision limit based on how many vectors are in database
    collisionLimit = 0.1 * indexL + 1000;
  }

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

  if (verbose)
    std::cout << "FINISHED" << std::endl;

  // Print out smallest vector in database
  int minIndex = -1;
  int min = INT_MAX;
  for (int i = 0; i < indexL; i++)
  {
    if (dbSqrNorm[i] < min && dbSqrNorm[i] > 0)
    {
      min = dbSqrNorm[i];
      minIndex = i;
    }
  }

  if (verbose)
  {
    std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;
    std::cout << "Number Of Vectors: " << indexL << std::endl;
    std::cout << "Collisions: " << collisionCount << std::endl;
    std::cout << std::fixed << "Squared Norm: " << dbSqrNorm[minIndex] << std::endl;
    std::cout << std::fixed << "Norm: " << sqrt(dbSqrNorm[minIndex]) << std::endl
              << std::endl;
    std::cout << std::defaultfloat << "Lattice Point: ";
    for (int x = 0; x < dimension; x++)
      std::cout << db[minIndex * dimension + x] << "  ";
    std::cout << std::endl
              << std::endl;
    std::cout << "Preimage: ";
    for (int x = 0; x < dimension; x++)
      std::cout << dbPreimage[minIndex * dimension + x] << "  ";
    std::cout << std::endl;
  }
  else
  {
    // Dimension, Klein Sample Range Divisor, Run Time (Microseconds), Smallest Norm Squared
    std::cout << "'Dimension: " << dimension << "',";
    std::cout << "'Sieve Time (Microseconds): " << duration.count() << "',";
    std::cout << "'Number Of Vectors: " << indexL << "',";
    std::cout << "'Smallest Vector Norm: " << sqrt(dbSqrNorm[minIndex]) << "'" << std::endl;
  }
}

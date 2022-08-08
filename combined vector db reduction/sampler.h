// Based off of the original 2011 implementation of gsieve by Panagiotis Voulgaris
// Can be found at https://cseweb.ucsd.edu//~pvoulgar/impl.html

// Sampling follows the work in:
// "Trapdoors for Hard Lattices and New Cryptographic Constructions"
// Craig Gentry, Chris Peikert, Vinod Vaikuntanathan
// See there for explanation of the variables s, sPrime, t etc...

#ifndef __SAMPLER__
#define __SAMPLER__

#include <stdio.h>
#include <math.h>
#include <cblas.h>

class KleinSampler
{
public:
  void Init(double *BIn, long rows, long columns);
  void Exit();
  void Sample(double *v);

private:
  long rows_;
  long columns_;
  double t_;
  double *mu_;
  double *coef_;
  double *sPrimeSquare_;
  long SampleZ(double c, double sSquare);
  void ComputeGramSchmidt(const double *B, double *mu, double *c);
};

void KleinSampler::Init(double *BIn, long rows, long columns)
{
  rows_ = rows;
  columns_ = columns;

  mu_ = (double *)calloc(sizeof(double), rows_ * columns_);
  double *bStarSquare = (double *)calloc(sizeof(double), rows_);

  ComputeGramSchmidt(BIn, mu_, bStarSquare);

  coef_ = (double *)malloc(sizeof(double) * rows_);
  sPrimeSquare_ = (double *)malloc(sizeof(double) * rows_);

  long maxStarSqrNorm = 0;
  for (int i = 0; i < rows_; ++i)
  {
    if ((long)bStarSquare[i] > maxStarSqrNorm)
      maxStarSqrNorm = (long)bStarSquare[i];
  }

  t_ = log(rows_);
  double sSquare = maxStarSqrNorm * log(rows_);

  for (int i = 0; i < rows_; ++i)
    sPrimeSquare_[i] = sSquare / (double)bStarSquare[i];

  free(bStarSquare);
}

void KleinSampler::Exit()
{
  free(mu_);
  free(coef_);
  free(sPrimeSquare_);
}

// Samples a vector in the preimage of the basis (Premultiply by basis to get lattice point)
void KleinSampler::Sample(double *v)
{
  memset(coef_, 0, rows_ * sizeof(double));

  for (int i = rows_ - 1; i >= 0; --i)
  {
    coef_[i] = SampleZ(coef_[i], sPrimeSquare_[i]);

    for (int j = 0; j < i; j++)
      coef_[j] -= (coef_[i] * mu_[i * rows_ + j]);
  }

  for (size_t x = 0; x < rows_; x++)
    v[x] = coef_[x];
}

long KleinSampler::SampleZ(double c, double sSquare)
{
  double s = sqrt(sSquare);
  long minC = floor(c - s * t_);
  long maxC = ceil(c + s * t_);
  long x;
  double rho;
  while (true)
  {
    x = minC + round((maxC - minC) * double(rand()) / RAND_MAX);
    rho = exp(-M_PI * (x - c) * (x - c) / sSquare);
    if ((double(rand()) / RAND_MAX) <= rho)
      return x;
  }
}

// Mat, mat vec
void KleinSampler::ComputeGramSchmidt(const double *B, double *mu, double *c)
{
  double *B1 = (double *)malloc(sizeof(double) * rows_ * columns_);
  double *b = (double *)malloc(sizeof(double) * columns_);
  double *buf = (double *)calloc(sizeof(double), columns_);

  // Copy B into B1
  for (int i = 0; i < columns_; i++)
    for (int j = 0; j < rows_; j++)
      B1[i * rows_ + j] = B[i * rows_ + j];

  for (int i = 0; i < columns_; i++)
    b[i] = cblas_ddot(rows_, &B1[i * rows_], 1, &B1[i * rows_], 1);

  const double precision = 150;
  const double bound = pow(2, 2 * long(0.15 * precision));
  const double bound2 = pow(2, 2 * precision);

  for (int k = 0; k < columns_; k++)
  {
    double s, t;

    if (k > 0)
      buf[0] = mu[k * rows_] * c[0];

    for (int j = 0; j < k; j++)
    {
      s = cblas_ddot(rows_, &B1[k * rows_], 1, &B1[j * rows_], 1);

      double t1 = s * s * bound;
      t = b[k] * b[j];

      if (t >= bound2 && t >= t1)
        s = cblas_ddot(rows_, &B[k * rows_], 1, &B[j * rows_], 1);

      t1 = 0;
      for (int i = 0; i < j; i++)
      {
        t = mu[j * rows_ + i] * buf[i];
        t1 += t;
      }

      t = s - t1;
      buf[j] = t;
      mu[k * rows_ + j] = t / c[j];
    }

    s = 0;
    for (int j = 0; j < k; j++)
      s += mu[k * rows_ + j] * buf[j];

    c[k] = b[k] - s;
  }

  free(B1);
  free(b);
  free(buf);
}

#endif
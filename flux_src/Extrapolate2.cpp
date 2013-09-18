// Extrapolate2.cpp

#include "Flux.h"

// ##############################################
// 1D extrapolation function with nonuniform grid
// ##############################################
double Flux::Extrapolate2 (int I, double *X, gsl_matrix *Y, double x, int k)
{
  int    i0 = I-2;
  double sm = (x-X[i0  ]) * (x-X[i0+1]) /(X[i0-1]-X[i0  ]) /(X[i0-1]-X[i0+1]);
  double s0 = (x-X[i0-1]) * (x-X[i0+1]) /(X[i0  ]-X[i0-1]) /(X[i0  ]-X[i0+1]);
  double s1 = (x-X[i0-1]) * (x-X[i0  ]) /(X[i0+1]-X[i0-1]) /(X[i0+1]-X[i0  ]);

  return 
    + sm * gsl_matrix_get (Y, i0-1, k) 
    + s0 * gsl_matrix_get (Y, i0,   k) 
    + s1 * gsl_matrix_get (Y, i0+1, k); 
}

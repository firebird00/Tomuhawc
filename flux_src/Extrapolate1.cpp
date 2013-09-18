// Extrapolate1.cpp

#include "Flux.h"

// ##############################################
// 1D extrapolation function with nonuniform grid
// ##############################################
double Flux::Extrapolate1 (int I, double *X, double *Y, double x)
{
  int    i0 = I-2;
  double sm = (x-X[i0  ]) * (x-X[i0+1]) /(X[i0-1]-X[i0  ]) /(X[i0-1]-X[i0+1]);
  double s0 = (x-X[i0-1]) * (x-X[i0+1]) /(X[i0  ]-X[i0-1]) /(X[i0  ]-X[i0+1]);
  double s1 = (x-X[i0-1]) * (x-X[i0  ]) /(X[i0+1]-X[i0-1]) /(X[i0+1]-X[i0  ]);
      
  return sm * Y[i0-1] + s0 * Y[i0] + s1 * Y[i0+1];
}

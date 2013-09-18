// FTCos.cpp

#include "Flux.h"

// ###########################################
// Function to return cosine Fourier transform
// ###########################################
void Flux::FTCos (int j, double in[], gsl_matrix *out)
{
  double hh = 1. /double (K-1);

  for (int m = 0; m < K/2; m++)
    {
      double sum = 0.;
      
      for (int k = 0; k < K; k++)
	{
	  double fac = (k == 0 || k == K-1) ? 0.5: 1.;
	  double   t = hh * double(k);

	  sum += fac * in[k] * cos(double(m) * th[k]) 
	    * M_PI * (0.5*M_PI/tc)*cos(0.5*M_PI*t/tc)/sin(0.5*M_PI/tc); 
	}
      sum *= hh * (2./M_PI);
      gsl_matrix_set (out, j, m, sum);
    }
  gsl_matrix_set (out, j, 0, gsl_matrix_get (out, j, 0)/2.);
}


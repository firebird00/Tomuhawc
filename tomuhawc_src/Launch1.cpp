// Launch1.cpp

#include "Tomuhawc.h"

// ####################################################################################
// Function to initialize single solution vector at flux surface label r
//
// r       ... flux surface label
// i       ... dominant poloidal harmonic is mpol[i]
// Y[dim1] ... single solution vector
//
// side                .. number of sideband harmonics
// vac                 .. total number of rational surfaces (including vacuum surfaces)
// num                 .. number of rational surfaces in plasma
// dim  = vac + 2*side .. dimension of coupling matrices
// dim1 = 2*dim        .. dimension of single solution vector 
//
// ####################################################################################
void Thawc::Launch1 (double r, int i, double Y[])
{
  double R   = r;
  double Q   = gsl_spline_eval (Sq, R, Aq);
  double M   = double(mpol[i]);
  double N   = double(ntor);
  double MNQ = M - N *Q;

  for (int k = 0; k < dim1; k++)
    Y[k] = 0.;

  if (mpol[i] == 0)
    Y[dim+i] = 1.;
  else if (mpol[i] > 0)
    {
      Y[i]     = 1.;
      Y[dim+i] = + Y[i] *MNQ /M;
    }
  else if (mpol[i] < 0)
    {
      Y[i]     = 1.;
      Y[dim+i] = - Y[i] *MNQ /M;
    }
}

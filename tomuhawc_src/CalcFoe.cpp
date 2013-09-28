// CalcFoe.cpp

#include "Tomuhawc.h"

// ##############################################################################################
// Function to calculate Foe-matrix
//
// dPsi(vac, dim+vac) .. coefficient of small solution
//  Psi(i, j)            ... coefficient at ith rational surface due to jth solution
// x(dim, vac)        .. boundary condition vector
// Foe(vac, vac)      .. Foe matrix
//
// side                  .. number of sideband harmonics
// vac                   .. total number of rational surfaces (including vacuum surfaces)
// num                   .. number of rational surfaces in plasma
// dim  = vac + 2*side   .. dimension of coupling matrices
// dim1 = 2*dim          .. dimension of single solution vector 
//
// ##############################################################################################
void Thawc::CalcFoe (gsl_matrix *dPsi, gsl_matrix *x, gsl_matrix *Foe)
{
  for (int j = 0; j < vac; j++)
    {
      for (int i = 0; i < vac; i++)
	{
	  double psi = gsl_matrix_get (dPsi, i, dim+j);
	  for (int k = 0; k < dim; k++)
	    psi += gsl_matrix_get (dPsi, i, k) * gsl_matrix_get (x, k, j);
	  gsl_matrix_set (Foe, i, j, psi);
	}
    }
}

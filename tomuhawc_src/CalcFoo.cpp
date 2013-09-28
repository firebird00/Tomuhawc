// CalcFoo.cpp

#include "Tomuhawc.h"

// ##############################################################################################
// Function to calculate Foo-matrix
//
// Psi(vac, dim+vac) .. coefficient of large solution
//  Psi(i, j)           ... coefficient at ith rational surface due to jth solution
// x(dim, vac)       .. boundary condition vector
// Foo(vac, vac)     .. Foo matrix
//
// side                  .. number of sideband harmonics
// vac                   .. total number of rational surfaces (including vacuum surfaces)
// num                   .. number of rational surfaces in plasma
// dim  = vac + 2*side   .. dimension of coupling matrices
// dim1 = 2*dim          .. dimension of single solution vector 
//
// ##############################################################################################
void Thawc::CalcFoo (gsl_matrix *Psi, gsl_matrix *x, gsl_matrix *Foo)
{
  for (int j = 0; j < vac; j++)
    {
      for (int i = 0; i < vac; i++)
	{
	  double psi = gsl_matrix_get (Psi, i, dim+j);
	  for (int k = 0; k < dim; k++)
	    psi += gsl_matrix_get (Psi, i, k) * gsl_matrix_get (x, k, j);
	  gsl_matrix_set (Foo, i, j, psi);
	}
    }
}

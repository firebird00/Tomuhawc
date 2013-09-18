// Launch.cpp

#include "Tomuhawc.h"

// ####################################################################################
// Function to initialize multiple solution vectors at flux surface label r
//
// r                               .. flux surface label
// YY  (dim1, dim+vac)             .. solution vector at radius r
//  YY(i=0,  dim -1;k=0,dim-1)       ... psi: i     .. poloidal harmonic index
//  YY(i=dim,dim1-1;k=0,dim-1)       ... Z  : i-dim .. poloidal harmonic index
//                                            k     .. index of independent solutions
//                                                      launched from axis
//  YY(i=0,  dim -1;j=dim,dim+vac-1) ... psi: i     .. poloidal harmonic index
//  YY(i=dim,dim1-1;j=dim,dim+vac-1) ... Z  : i-dim .. poloidal harmonic index
//                                            j     .. index of small solutions
//                                                      launched from rational surfaces 
//
// side                .. number of sideband harmonics
// vac                 .. total number of rational surfaces (including vacuum surfaces)
// num                 .. number of rational surfaces in plasma
// dim  = vac + 2*side .. dimension of coupling matrices
// dim1 = 2*dim        .. dimension of single solution vector 
//
// ####################################################################################
void Thawc::Launch (double r, gsl_matrix *YY)
{
  double *Y = new double[dim1];

  // Launch independent solutions from magnetic axis
  for (int i = 0; i < dim; i++)
    {
      Launch1 (r, i, Y);

      for (int k = 0; k < dim1; k++)
	gsl_matrix_set (YY, k, i, Y[k]);
    }

  // Initialize small solutions launched from rational surfaces to zero
  for (int i = 0; i < vac; i++)
    {
      double zero = 0.;
      
      for (int k = 0; k < dim1; k++)
	gsl_matrix_set (YY, k, dim+i, zero);
    }

  delete[] Y;
}

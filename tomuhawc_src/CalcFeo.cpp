// CalcFeo.cpp

#include "Tomuhawc.h"

// ##############################################################################################
// Function to calculate Feo-matrix
//
// dPsi(vac, dim+vac) .. coefficient of small solution
//  Psi(i, j)                         ... coefficient at ith rational surface due to jth solution
// x(dim, vac)        .. boundary condition vector
// Feo(vac, vac)      .. Feo matrix
// YY(dim1, dim+vac)  .. solution vector
//  YY(i=0,  dim -1;k=0,dim-1)        ... psi: i     .. poloidal harmonic index
//  YY(i=dim,dim1-1;k=0,dim-1)        ... Z  : i-dim .. poloidal harmomic index 
//                                             k     .. index of independent solutions
//                                                       launched from axis
//  YY(i=0,  dim -1;j=dim,dim+vac-1)  ... psi: i     .. poloidal harmonic index
//  YY(i=dim,dim1-1;j=dim,dim+vac-1)  ... Z  : i-dim .. poloidal harmonic index
//                                             j     .. index of small solutions
//                                                       launched from rational surfaces
// Boundary condition:
// Sum_k = 0,dim-1 YY(i, k) * x (k, j) = - YY(i, dim+j) for i = 0,dim-1 and j = 0,vac-1
//
// side                  .. number of sideband harmonics
// vac                   .. total number of rational surfaces (including vacuum surfaces)
// num                   .. number of rational surfaces in plasma
// dim  = vac + 2*side   .. dimension of coupling matrices
// dim1 = 2*dim          .. dimension of single solution vector 
//
// ##############################################################################################
void Thawc::CalcFeo (gsl_matrix *dPsi, gsl_matrix *x, gsl_matrix *Feo)
{
  for (int j = 0; j < vac; j++)
    {
      for (int i = 0; i < vac; i++)
	{
	  double psi = gsl_matrix_get (dPsi, i, dim+j);
	  for (int k = 0; k < dim; k++)
	    psi += gsl_matrix_get (dPsi, i, k) * gsl_matrix_get (x, k, j);
	  gsl_matrix_set (Feo, i, j, psi);
	}
    }
}

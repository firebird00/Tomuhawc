// Jump.cpp

#include "Tomuhawc.h"

// ####################################################################################
// Function to evolve multiple solution vectors across jth rational surface
// while calculating reconnected fluxes
//
// r                               .. flux-surface label
// j                               .. rational surface index (j = 0, vac-1)
// YY  (dim1, dim+vac)             .. solution vector at radius r
//  YY(i=0,  dim -1;k=0,dim-1)       ... psi: i     .. poloidal harmonic index
//  YY(i=dim,dim1-1;k=0,dim-1)       ... Z  : i-dim .. poloidal harmonic index
//                                            k     .. index of independent solutions
//                                                      launched from axis
//  YY(i=0,  dim -1;j=dim,dim+vac-1) ... psi: i     .. poloidal harmonic index
//  YY(i=dim,dim1-1;j=dim,dim+vac-1) ... Z  : i-dim .. poloidal harmonic index
//                                            j     .. index of small solutions
//                                                      launched from rational surfaces
// Psi (vac,  dim+vac)             .. coefficients of large solution
//  Psi (j=0, vac-1;i=0+dim+vac-1)   ... Psi at jth rational surface due to ith
//                                        independent solution or small solution
//                                        launched from (i-dim)th rational surface
// dPsi(vac,  dim+vac)             .. coefficients of small solution
//  dPsi(j=0, vac-1;i=0+dim+vac-1)   ... dPsi at jth rational surface due to ith
//                                        independent solution or small solution
//                                        launched from (i-dim)th rational surface 
// _tear                           .. 1/0 for tearing/twisting parity
// interactive                     .. set if in interactive mode
//
// side                .. number of sideband harmonics
// vac                 .. total number of rational surfaces (including vacuum surfaces)
// num                 .. number of rational surfaces in plasma
// dim  = vac + 2*side .. dimension of coupling matrices
// dim1 = 2*dim        .. dimension of single solution vector 
//
// ####################################################################################
void Thawc::Jump (double& r, int j, gsl_matrix *YY, gsl_matrix *Psi, gsl_matrix *dPsi, int _tear, int interactive)
{
  double *Y = new double[dim1];

  double delta = Rres[j] - r;
  for (int i = 0; i < dim+vac; i++)
    {
      for (int k = 0; k < dim1; k++)
	Y[k] = gsl_matrix_get (YY, k, i);
      
      double psi, dpsi;
      Jump1 (delta, i, j, Y, psi, dpsi, _tear, interactive);
      gsl_matrix_set (Psi,  j, i, psi);
      gsl_matrix_set (dPsi, j, i, dpsi);
      
      for (int k = 0; k < dim1; k++)
	gsl_matrix_set (YY, k, i, Y[k]);
    }

  r += 2.*delta;
  
  delete[] Y;
}


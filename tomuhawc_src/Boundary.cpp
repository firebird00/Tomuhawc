// Boundary.cpp

#include "Tomuhawc.h"

// #####################################################################################
// Function to apply boundary conditions at wall
//
// Ya(dim1, dim+vac)                 ... solution vector at r = a
//  Ya(i=0,  dim -1;k=0,dim-1)        ... psi: i     .. poloidal harmonic index
//  Ya(i=dim,dim1-1;k=0,dim-1)        ... Z  : i-dim .. poloidal harmomic index 
//                                             k     .. index of independent solutions
//                                                       launched from axis
//  Ya(i=0,  dim -1;j=dim,dim+vac-1)  ... psi: i     .. poloidal harmonic index
//  Ya(i=dim,dim1-1;j=dim,dim+vac-1)  ... Z  : i-dim .. poloidal harmonic index
//                                             j     .. index of small solutions
//                                                       launched from rational surfaces
// QP(dim1, dim+vac)                 ... q, p coefficients of solution vector at r = a
//  QP(i=0,  dim -1;k=0,dim-1)        ... q  : i     .. poloidal harmonic index
//  QP(i=dim,dim1-1;k=0,dim-1)        ... p  : i-dim .. poloidal harmonic index
//                                             k     .. index of independent solutions
//                                                      launched from axis
//  QP(i=0,  dim -1;j=dim,dim+vac-1)  ... q  : i     .. poloidal harmonic index
//  QP(i=dim,dim1-1;j=dim,dim+vac-1)  ... p  : i-dim .. poloidal harmonic index
//                                             j     .. index of small solutions
//                                                      launched from rational surfaces 
// x(dim, vac)                       ... boundary condition vector
// interactive                       ... set if in interactive mode
//
// Fixed boundary condition:
// Sum_k = 0,dim-1 Ya(i, k) * x(k, j) = - Ya(i, dim+j) for i = 0,dim-1 and j = 0,vac-1
//
// Free boundary condition:
// Sum_k = 0,dim-1 QP(i, k) * x(k, j) = - QP(i, dim+j) for i = 0,dim-1 and j = 0,vac-1
//
// side                .. number of sideband harmonics
// vac                 .. total number of rational surfaces (including vacuum surfaces)
// num                 .. number of rational surfaces in plasma
// dim  = vac + 2*side .. dimension of coupling matrices
// dim1 = 2*dim        .. dimension of single solution vector 
//
// ######################################################################################
void Thawc::Boundary (gsl_matrix *Ya, gsl_matrix *QP, gsl_matrix *x, int interactive)
{
  gsl_matrix      *Mat = gsl_matrix_alloc (dim, dim);
  gsl_matrix      *Y1  = gsl_matrix_alloc (dim, vac);
  gsl_vector      *Rhs = gsl_vector_alloc (dim);
  gsl_vector      *xx  = gsl_vector_alloc (dim);
  gsl_permutation *p   = gsl_permutation_alloc (dim);

  for (int i = 0; i < dim; i++)
    for (int k = 0; k < dim; k++)
      if (free)
	gsl_matrix_set (Mat, i, k, gsl_matrix_get (QP, i, k));
      else
	gsl_matrix_set (Mat, i, k, gsl_matrix_get (Ya, i, k));

  int s;
  gsl_linalg_LU_decomp (Mat, p, &s);

  for (int j = 0; j < vac; j++)
    {
      for (int i = 0; i < dim; i++)
	if (free)
	  gsl_vector_set (Rhs, i, - gsl_matrix_get (QP, i, dim+j));
	else
	  gsl_vector_set (Rhs, i, - gsl_matrix_get (Ya, i, dim+j));
     
      gsl_linalg_LU_solve (Mat, p, Rhs, xx);

      for (int k = 0; k < dim; k++)
	gsl_matrix_set (x, k, j, gsl_vector_get (xx, k));

      for (int i = 0; i < dim; i++)
	{
	  double y;
	  if (free)
	    y = gsl_matrix_get (QP, i, dim+j);
	  else
	    y = gsl_matrix_get (Ya, i, dim+j);
	    
	  for (int k = 0; k < dim; k++)
	    if (free)
	      y += gsl_matrix_get (QP, i, k) * gsl_matrix_get (x, k, j);
	    else
	      y += gsl_matrix_get (Ya, i, k) * gsl_matrix_get (x, k, j);
	  y = (fabs(y) < Eta) ? 0. : y;
	  gsl_matrix_set (Y1, i, j, y);
	}
    }
 
  if (interactive) LogY1 (Y1);

  gsl_matrix_free (Mat);
  gsl_matrix_free (Y1);
  gsl_vector_free (Rhs);
  gsl_vector_free (xx);
  gsl_permutation_free (p);
}

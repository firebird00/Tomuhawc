// Rhs.cpp

#include "Tomuhawc.h"

// ######################################################################################
// Function to evaluate right-hand sides of mode equations
//
// r          .. flux surface label
// y[dimN]    .. solution vector
// dydr[dimN] .. derivative of solution vector
//
//   i = 0,dim-1         .. independent solutions launched from axis
//   i = dim,dim+vac-1   .. small solutions launched from rational surfaces
//
//      y[dim1*i     +0] .. y[dim1*i     +dim-1]  - psi
//      y[dim1*i+dim +0] .. y[dim1*i+dim +dim-1]  - Z
//
// side                  .. number of sideband harmonics
// vac                   .. total number of rational surfaces (including vacuum surfaces)
// num                   .. number of rational surfaces in plasma
// dim  = vac + 2*side   .. dimension of coupling matrices
// dim1 = 2*dim          .. dimension of single solution vector 
// dimN = dim1*(dim+vac) .. dimension of composite solution vector
//
// ######################################################################################
int Thawc::Rhs (double r, const double y[], double dydr[], void *)
{
  Couple (r);
  double q, qp, nq, nqp;
  q   = gsl_spline_eval (Sq,  r, Aq);
  qp  = gsl_spline_eval (Sqp, r, Aqp);
  if (r > r1 && r < ra-eps)
    qp = (gsl_spline_eval (Sq, r+eps, Aq) - gsl_spline_eval (Sq, r-eps, Aq)) /(2.*eps);
  nq  = double(ntor) * q;
  nqp = double(ntor) * qp;
  
  for (int i = 0; i < dim+vac; i++)
    for (int k = 0; k < dim; k++)
      {
	double mkmnq = double(mpol[k]) - nq;
	
	// Set dpsi/dr
	dydr[dim1*i+k] = 0.;
	for (int j = 0; j < dim; j++)
	  {
	    double mjmnq = double(mpol[j]) - nq;
	    dydr[dim1*i+k] += (gsl_matrix_get (Mmat, k, j) * y[dim1*i+j] 
				+ gsl_matrix_get (Lmat, k, j) * y[dim1*i+dim+j]) /mjmnq/r;
	  }
	
	// Set dZ/dr
	dydr[dim1*i+dim+k] = - nqp * y[dim1*i+dim+k] /mkmnq;
	for (int j = 0; j < dim; j++)
	  {
	    double mjmnq = double(mpol[j]) - nq;
	    dydr[dim1*i+dim+k] += (gsl_matrix_get (Pmat, k, j) * y[dim1*i+j] 
				    + gsl_matrix_get (Nmat, k, j) * y[dim1*i+dim+j]) /mjmnq/r;
	  }
      }

  return GSL_SUCCESS;
}

// Segment.cpp

#include "Tomuhawc.h"

// ######################################################################################
// Function to integrate multiple solution vectors to r = rx
// via adaptive step method
//
// r                     .. flux-surface label
// rx                    .. integrate until r = rx
// YY(dim1, dim+vac)     .. solution vector at radius r
//  YY(i=0,  dim -1;k=0,dim-1)       ... psi: i     .. poloidal harmonic index
//  YY(i=dim,dim1-1;k=0,dim-1)       ... Z  : i-dim .. poloidal harmonic index
//                                            k     .. index of independent solutions
//                                                      launched from axis
//  YY(i=0,  dim -1;j=dim,dim+vac-1) ... psi: i     .. poloidal harmonic index
//  YY(i=dim,dim1-1;j=dim,dim+vac-1) ... Z  : i-dim .. poloidal harmonic index
//                                            j     .. index of small solutions
//                                                      launched from rational surfaces 
// flag                  .. truncation error flag
// edge                  .. set if rx is wall radius
// interactive           .. set if in iteractive mode
// _adap                 .. name of adaptive step monitor file
// _soln                 .. name of solution file
// _edge                 .. name of edge data file
//
// side                  .. number of sideband harmonics
// vac                   .. total number of rational surfaces (including vacuum surfaces)
// num                   .. number of rational surfaces in plasma
// dim  = vac + 2*side   .. dimension of coupling matrices
// dim1 = 2*dim          .. dimension of single solution vector 
// dimN = dim1*(dim+vac) .. dimension of composite solution vector
//
// ######################################################################################
void Thawc::Segment (double &r, double rx, gsl_matrix *YY, 
		     int flag, int edge, int interactive,
		     char *_adap, char *_soln, char *_edge)
{
  double *Y = new double[dimN];
  
  for (int i = 0; i < dim+vac; i++)
    for (int k = 0; k < dim1; k++)
      Y[dim1*i+k] = gsl_matrix_get (YY, k, i);
  
  Segment1 (r, rx, Y, flag, edge, interactive, _adap, _soln, _edge);
  
  for (int i = 0; i < dim+vac; i++)
    for (int k = 0; k < dim1; k++)
      gsl_matrix_set (YY, k, i, Y[dim1*i+k]);
  
  delete[] Y;
}


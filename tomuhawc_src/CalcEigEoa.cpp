// CalcEigEoa.cpp

#include "Tomuhawc.h"

// ########################################################
// Function to calculate eigenfunctions of Eo-matrix at r=a
// ########################################################
void Thawc::CalcEigEoa (gsl_matrix *Eo)
{
  FILE *zfo  = OpenFiler ("Stage3/ZFoa.out");
  FILE *xio  = OpenFile  ("Stage3/xio.out");

  double r;
  gsl_matrix *ZFO = gsl_matrix_alloc (dim, vac);
  fscanf  (zfo, "%lf", &r);
  
  // Read in eigenfunctions of F-matrix
  for (int j = 0; j < vac; j++)
    for (int k = 0; k < dim; k++)
      {
	double val;
	fscanf (zfo, "%lf", &val);
	gsl_matrix_set (ZFO, k, j, val); 
      }
  
  // Calculate eigenfunctions of E-matrix
  double sumz;
  for (int j = 0; j < num; j++)
    for (int k = 0; k < dim; k++)
      {
	double valp = 0., valz = 0.;
	for (int n = 0; n < num; n++)
	  valz += gsl_matrix_get (ZFO, k, n) * gsl_matrix_get (Eo, n, j);
	valz *= (a/qa);
	fprintf (xio, " %e", valz);
      }	
  fprintf (xio, "\n");
  
  gsl_matrix_free (ZFO); 

  fclose (xio); 
}

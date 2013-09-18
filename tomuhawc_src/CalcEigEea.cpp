// CalcEigEea.cpp

#include "Tomuhawc.h"

// #########################################################
// Function to calculate eigenfunctions of Ee-matrix at wall
// #########################################################
void Thawc::CalcEigEea (gsl_matrix *Ee)
{
  FILE *zfe  = OpenFiler ("Stage3/ZFea.out");
  FILE *xie  = OpenFile  ("Stage3/xie.out");

  double r;
  gsl_matrix *ZFE = gsl_matrix_alloc (dim, vac);
  fscanf (zfe, "%lf", &r);
  
  // Read in eigenfunctions of F-matrix
  for (int j = 0; j < vac; j++)
    for (int k = 0; k < dim; k++)
      {
	double val;
	fscanf (zfe, "%lf", &val);
	gsl_matrix_set (ZFE, k, j, val); 
      }
  
  // Calculate eigenfunctions of E-matrix
  double sumz;
  for (int j = 0; j < num; j++)
    for (int k = 0; k < dim; k++)
      {
	double valz = 0.;
	for (int n = 0; n < num; n++)
	  valz += gsl_matrix_get (ZFE, k, n) * gsl_matrix_get (Ee, n, j);
	valz *= (a/qa);
	fprintf (xie, " %e", valz);
      }	
  fprintf (xie, "\n");
  
  gsl_matrix_free (ZFE); 

  fclose (xie); 
}

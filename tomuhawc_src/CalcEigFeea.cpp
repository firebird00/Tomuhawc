// CalcEigFeea.cpp

#include "Tomuhawc.h"

// #########################################################
// Function to calculate eigenfunctions of Fee-matrix at r=a
// #########################################################
void Thawc::CalcEigFeea (gsl_matrix *x, gsl_matrix *Ya)
{
  double  r;
  double *Y = new double[dimN];
  
  // Form linear combination of solution vectors 
  // that satisfies boundary condition at r = a
  FILE *psi = OpenFile ("Stage3/psiFea.out");
  FILE *Z   = OpenFile ("Stage3/ZFea.out");

  gsl_matrix *YYY = gsl_matrix_alloc (dim1, dim+vac);
  gsl_matrix *PFU = gsl_matrix_alloc (dim,  vac);
  gsl_matrix *ZFU = gsl_matrix_alloc (dim,  vac);
  
  fprintf (psi, "%e", r);
  fprintf (Z,   "%e", r);
  
  double q  = gsl_spline_eval (Sq, r, Aq);
  double nq = double(ntor) * q;
  
  // Calculate eigenfunctions
  for (int j = 0; j < vac; j++)
    for (int k = 0; k < dim; k++)
      {
	double valp = gsl_matrix_get (Ya, k,     dim+j);
	double valz = gsl_matrix_get (Ya, dim+k, dim+j);
	for (int i = 0; i < dim; i++)
	  {
	    valp += gsl_matrix_get (Ya, k,     i) * gsl_matrix_get (x, i, j);
	    valz += gsl_matrix_get (Ya, dim+k, i) * gsl_matrix_get (x, i, j);
	  }
	gsl_matrix_set (PFU, k, j, valp);
	gsl_matrix_set (ZFU, k, j, valz);
	fprintf (psi, " %e", gsl_matrix_get (PFU, k, j));
	fprintf (Z,   " %e", gsl_matrix_get (ZFU, k, j));
      }
  fprintf (psi, "\n"); fprintf (Z, "\n");
  
  gsl_matrix_free (PFU);
  gsl_matrix_free (ZFU);
  fclose (psi); fclose (Z); 
  
  delete[] Y;
}

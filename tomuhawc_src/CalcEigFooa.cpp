// CalcEigFooa.cpp

#include "Tomuhawc.h"

// ##########################################################
// Function to calculate eigenfunctions of Foo-matrix at wall
// ##########################################################
void Thawc::CalcEigFooa (gsl_matrix *x)
{
  double  r;
  double *Y = new double[dimN];
  
  // Read in solution vectors
  FILE *soln = OpenFiler ("Stage3/Yoa.out");
  fscanf (soln, "%lf", &r);
  for (int k = 0; k < dimN; k++)
    fscanf (soln, "%lf", &Y[k]);
  fclose (soln);

  // Form linear combination of solution vectors 
  // that satisfies boundary condition at r = a
  FILE *psi = OpenFile ("Stage3/psiFoa.out");
  FILE *Z   = OpenFile ("Stage3/ZFoa.out");

  gsl_matrix *YYY = gsl_matrix_alloc (dim1, dim+vac);
  gsl_matrix *PFU = gsl_matrix_alloc (dim,  vac);
  gsl_matrix *ZFU = gsl_matrix_alloc (dim,  vac);
  
  fprintf (psi, "%e", r);
  fprintf (Z,   "%e", r);
  
  double q  = gsl_spline_eval (Sq, r, Aq);
  double nq = double(ntor) * q;
  
  for (int i = 0; i < dim+vac; i++)
    for (int k = 0; k < dim1; k++)
      gsl_matrix_set (YYY, k, i, Y[dim1*i+k]);
  
  // Calculate eigenfunctions
  for (int j = 0; j < vac; j++)
    for (int k = 0; k < dim; k++)
      {
	double valp = gsl_matrix_get (YYY, k,     dim+j);
	double valz = gsl_matrix_get (YYY, dim+k, dim+j);
	for (int i = 0; i < dim; i++)
	  {
	    valp += gsl_matrix_get (YYY, k,     i) * gsl_matrix_get (x, i, j);
	    valz += gsl_matrix_get (YYY, dim+k, i) * gsl_matrix_get (x, i, j);
	  }
	gsl_matrix_set (PFU, k, j, valp);
	gsl_matrix_set (ZFU, k, j, valz);
	fprintf (psi, " %e", gsl_matrix_get (PFU, k, j));
	fprintf (Z,   " %e", gsl_matrix_get (ZFU, k, j));
      }
  fprintf (psi, "\n"); fprintf (Z, "\n");
  
  gsl_matrix_free (YYY);
  gsl_matrix_free (PFU);
  gsl_matrix_free (ZFU);
  fclose (psi); fclose (Z); 
  
  delete[] Y;
}

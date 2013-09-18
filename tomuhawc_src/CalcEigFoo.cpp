// CalcEigFoo.cpp

#include "Tomuhawc.h"

// ##################################################
// Function to calculate eigenfunctions of Foo-matrix
// ##################################################
void Thawc::CalcEigFoo (gsl_matrix *x)
{
  // Find number of radial grid points
  FILE *soln = OpenFiler ("Stage3/solno.out");
  int c, npts = 0;
  do 
    {
      c = fgetc (soln);
      if (c == '\n') 
	npts++;
    }   
  while (c != EOF);
  fclose (soln);

  double     *rr = new double[npts];
  gsl_matrix *Y  = gsl_matrix_alloc (npts, dimN);
  
  // Read in solution vectors
  soln = OpenFiler ("Stage3/solno.out");
  for (int l = 0; l < npts; l++)
    {
      fscanf (soln, "%lf", &rr[l]);
      for (int k = 0; k < dimN; k++)
	  {
	    double val;
	    fscanf (soln, "%lf", &val);
	    gsl_matrix_set (Y, l, k, val); 
	  }
    }
  fclose (soln);

  // At each grid point, form linear combination of solution vectors 
  // that satisfies boundary condition at r = a
  FILE *psi = OpenFile ("Stage3/psiFo.out");
  FILE *Z   = OpenFile ("Stage3/ZFo.out");
  for (int l = 0; l < npts; l++)
    {
      gsl_matrix *YYY = gsl_matrix_alloc (dim1, dim+vac);
      gsl_matrix *PFO = gsl_matrix_alloc (dim,  vac);
      gsl_matrix *ZFO = gsl_matrix_alloc (dim,  vac);

      fprintf (psi, "%e", rr[l]);
      fprintf (Z,   "%e", rr[l]);
  
      double q  = gsl_spline_eval (Sq,  rr[l], Aq);
      double nq = double(ntor) * q;

      for (int i = 0; i < dim+vac; i++)
	for (int k = 0; k < dim1; k++)
	  gsl_matrix_set (YYY, k, i, gsl_matrix_get (Y, l, dim1*i+k));

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
	    gsl_matrix_set (PFO, k, j, valp);
	    gsl_matrix_set (ZFO, k, j, valz);
	    fprintf (psi, " %e", gsl_matrix_get (PFO, k, j));
	    fprintf (Z,   " %e", gsl_matrix_get (ZFO, k, j));
	  }
      fprintf (psi, "\n"); fprintf (Z, "\n");

      gsl_matrix_free (YYY);
      gsl_matrix_free (PFO);
      gsl_matrix_free (ZFO);
    }
  fclose (psi); fclose (Z);

  delete[] rr;
  gsl_matrix_free (Y);
}

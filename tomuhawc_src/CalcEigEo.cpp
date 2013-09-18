// CalcEigEo.cpp

#include "Tomuhawc.h"

// #################################################
// Function to calculate eigenfunctions of Eo-matrix
// #################################################
void Thawc::CalcEigEo (gsl_matrix *Eo)
{
  // Find number of radial grid points
  FILE *pfu = OpenFiler ("Stage3/PsiFo.out");
  int c, npts = 0;
  do 
    {
      c = fgetc (pfu);
      if (c == '\n') 
	npts++;
    }   
  while (c != EOF);
  fclose (pfu);

  FILE   *pfo = OpenFiler ("Stage3/PsiFo.out");
  FILE   *zfo = OpenFiler ("Stage3/ZFo.out");
  FILE   *puo = OpenFile  ("Stage3/PsiUo.out");
  FILE   *zuo = OpenFile  ("Stage3/ZUo.out");
  double *rr  = new double[npts];
  for (int l = 0; l < npts; l++)
    {
      gsl_matrix *PFO = gsl_matrix_alloc (dim, vac);
      gsl_matrix *ZFO = gsl_matrix_alloc (dim, vac);

      fscanf  (pfo, "%lf", &rr[l]);
      fscanf  (zfo, "%lf", &rr[l]);
      fprintf (puo , "%e",  rr[l]);
      fprintf (zuo , "%e",  rr[l]);

      // Read in eigenfunctions of F-matrix
      for (int j = 0; j < vac; j++)
	for (int k = 0; k < dim; k++)
	  {
	    double val;
	    fscanf (pfo, "%lf", &val);
	    gsl_matrix_set (PFO, k, j, val); 
	    fscanf (zfo, "%lf", &val);
	    gsl_matrix_set (ZFO, k, j, val); 
	  }

      // Calculate eigenfunctions of E-matrix
      for (int j = 0; j < num; j++)
	{
	  for (int k = 0; k < dim; k++)
	    {
	      double valp = 0., valz = 0.;
	      for (int n = 0; n < num; n++)
		{
		  valp += gsl_matrix_get (PFO, k, n) * gsl_matrix_get (Eo, n, j);
		  valz += gsl_matrix_get (ZFO, k, n) * gsl_matrix_get (Eo, n, j);
		}
	      fprintf (puo, " %e", valp);
	      fprintf (zuo, " %e", valz);
	    }
	}
      fprintf (puo, "\n"); fprintf (zuo, "\n");
      
      gsl_matrix_free (PFO); 
      gsl_matrix_free (ZFO); 
    }
  fclose (puo); fclose (zuo);
  
  delete[] (rr);
}


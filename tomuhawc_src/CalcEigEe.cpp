// CalcEigEe.cpp

#include "Tomuhawc.h"

// #################################################
// Function to calculate eigenfunctions of Ee-matrix
// #################################################
void Thawc::CalcEigEe (gsl_matrix *Ee)
{
  // Find number of radial grid points
  FILE *pfu = OpenFiler ("Stage3/PsiFe.out");
  int c, npts = 0;
  do 
    {
      c = fgetc (pfu);
      if (c == '\n') 
	npts++;
    }   
  while (c != EOF);
  fclose (pfu);

  FILE   *pfe = OpenFiler ("Stage3/PsiFe.out");
  FILE   *zfe = OpenFiler ("Stage3/ZFe.out");
  FILE   *pue = OpenFile ("Stage3/PsiUe.out");
  FILE   *zue = OpenFile ("Stage3/ZUe.out");
  double *rr  = new double[npts];
  for (int l = 0; l < npts; l++)
    {
      gsl_matrix *PFE = gsl_matrix_alloc (dim, vac);
      gsl_matrix *ZFE = gsl_matrix_alloc (dim, vac);

      fscanf  (pfe, "%lf", &rr[l]);
      fscanf  (zfe, "%lf", &rr[l]);
      fprintf (pue , "%e",  rr[l]);
      fprintf (zue , "%e",  rr[l]);

      // Read in eigenfunctions of F-matrix
      for (int j = 0; j < vac; j++)
	for (int k = 0; k < dim; k++)
	  {
	    double val;
	    fscanf (pfe, "%lf", &val);
	    gsl_matrix_set (PFE, k, j, val); 
	    fscanf (zfe, "%lf", &val);
	    gsl_matrix_set (ZFE, k, j, val); 
	  }

      // Calculate eigenfunctions of E-matrix
      for (int j = 0; j < num; j++)
	{
	  for (int k = 0; k < dim; k++)
	    {
	      double valp = 0., valz = 0.;
	      for (int n = 0; n < num; n++)
		{
		  valp += gsl_matrix_get (PFE, k, n) * gsl_matrix_get (Ee, n, j);
		  valz += gsl_matrix_get (ZFE, k, n) * gsl_matrix_get (Ee, n, j);
		}
	      fprintf (pue, " %e", valp);
	      fprintf (zue, " %e", valz);
	    }
	}	
      fprintf (pue, "\n"); fprintf (zue, "\n");
      
      gsl_matrix_free (PFE); 
      gsl_matrix_free (ZFE); 
    }
  fclose (pue); fclose (zue);

  delete[] (rr);
}

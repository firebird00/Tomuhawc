// Fixup.cpp

#include "Tomuhawc.h"

// ##############################################################################################
// Function to perform fixup of multiple solution vectors
//
// YY(dim1, dim+vac)               .. solution vector 
//  YY(i=0,  dim -1;k=0,dim-1)        ... psi: i     .. poloidal harmonic index
//  YY(i=dim,dim1-1;k=0,dim-1)        ... Z  : i-dim .. poloidal harmonic index
//                                             k     .. index of independent solutions
//                                                       launched from axis
// Psi(vac, dim+vac)               .. coefficient of large solution
//  Psi(i, j)                         ... coefficient at ith rational surface due to jth solution
// dPsi(vac, dim+vac)              .. coefficient of small solution
//  Psi(i, j)                         ... coefficient at ith rational surface due to jth solution
// interactive                     .. set if in interactive mode
// _soln                           .. name of solution file
//
// side                            .. number of sideband harmonics
// vac                             .. total number of rational surfaces (including vacuum surfaces)
// num                             .. number of rational surfaces in plasma
// dim  = vac + 2*side             .. dimension of coupling matrices
// dim1 = 2*dim                    .. dimension of single solution vector 
// dimN = dim1*(dim+vac)           .. dimension of composite solution vector
//
// ##############################################################################################
void Thawc::Fixup (gsl_matrix *YY, gsl_matrix *Psi, gsl_matrix *dPsi, int interactive, 
		   char *_soln)
{
  FILE       *soln;
  int         c, npts;
  double      r, *rr, *Ye;
  gsl_matrix *Y;

  if (interactive)
    {
      // Find number of radial grid points in solution file
      soln = OpenFiler (_soln);
      npts = 0;
      do 
	{
	  c = fgetc (soln);
	  if (c == '\n') 
	    npts++;
	}   
      while (c != EOF);
      fclose (soln);
      
      rr = new double[npts];
      Y  = gsl_matrix_alloc (npts, dimN);
      
      // Read in solution vectors from solution file
      soln = OpenFiler (_soln);
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
    }

  // Find index of m=0 poloidal harmonic
  int i0 = -1;
  for (int i = 0; i < dim; i++)
    if (mpol[i] == 0) i0 = i;
  
  // Range of poloidal mode numbers does not include m=0
  if (i0 == -1)
    {
      for (int i = dim-1; i >= 0; i--)
	{
	  double yii = gsl_matrix_get (YY, i, i);
	  for (int j = i-1; j >= 0; j--)
	    {
	      double yij = gsl_matrix_get (YY, i, j);

	      for (int k = 0; k < dim1; k++)
		{
		  double yki = gsl_matrix_get (YY, k, i);
		  double ykj = gsl_matrix_get (YY, k, j);

		  // Fixup solution vector
		  double val = ykj - (yij /yii) *yki;
		  if (fabs(val) < Eta) val = 0.;
		  gsl_matrix_set (YY, k, j, val);
		  
		  // Fixup solution vectors from solution file
		  if (interactive)
		    for (int l = 0; l < npts; l++)
		      {
			double Yki = gsl_matrix_get (Y, l, dim1*i+k);
			double Ykj = gsl_matrix_get (Y, l, dim1*j+k);
			val = Ykj - (yij /yii) *Yki;
			if (fabs(val) < Eta) val = 0.;
			gsl_matrix_set (Y, l, dim1*j+k, val);
		      }
		}

	      // Fixup Psi and dPsi matrices
	      for (int l = 0; l < vac; l++)
		{
		  double val;
		  val = gsl_matrix_get (Psi,  l, j) - (yij /yii) *gsl_matrix_get (Psi,  l, i);
		  gsl_matrix_set (Psi, l, j, val);
		  val = gsl_matrix_get (dPsi, l, j) - (yij /yii) *gsl_matrix_get (dPsi, l, i);
		  gsl_matrix_set (dPsi, l, j, val);
		}
	    }
	}
    }
  // Range of poloidal mode numbers does include m=0
  else
    {
      for (int i = dim-1; i > i0; i--)
	{
	  double yii = gsl_matrix_get (YY, i, i);
	  for (int j = i-1; j >= 0; j--)
	    {
	      double yij = gsl_matrix_get (YY, i, j);

	      for (int k = 0; k < dim1; k++)
		{
		  double yki = gsl_matrix_get (YY, k, i);
		  double ykj = gsl_matrix_get (YY, k, j);

		  // Fixup solution vector
		  double val = ykj - (yij /yii) *yki;
		  if (fabs(val) < Eta) val = 0.;
		  gsl_matrix_set (YY, k, j, val);

		  // Fixup solution vectors from solution file
		  if (interactive) 
		    for (int l = 0; l < npts; l++)
		      {
			double Yki = gsl_matrix_get (Y, l, dim1*i+k);
			double Ykj = gsl_matrix_get (Y, l, dim1*j+k);
			val = Ykj - (yij /yii) *Yki;
			if (fabs(val) < Eta) val = 0.;
			gsl_matrix_set (Y, l, dim1*j+k, val);
		      }	
		}
	      
	      // Fixup Psi and dPsi matrices
	      for (int l = 0; l < vac; l++)
		{
		  double val;
		  val = gsl_matrix_get (Psi,  l, j) - (yij /yii) *gsl_matrix_get (Psi,  l, i);
		  gsl_matrix_set (Psi, l, j, val);
		  val = gsl_matrix_get (dPsi, l, j) - (yij /yii) *gsl_matrix_get (dPsi, l, i);
		  gsl_matrix_set (dPsi, l, j, val);
		}
	    }
	}
      for (int i = 0; i < i0; i++)
	{
	  double yii = gsl_matrix_get (YY, i, i);
	  for (int j = i+1; j < dim; j++)
	    {
	      double yij = gsl_matrix_get (YY, i, j);

	      for (int k = 0; k < dim1; k++)
		{
		  double yki = gsl_matrix_get (YY, k, i);
		  double ykj = gsl_matrix_get (YY, k, j);

		  // Fixup solution vector
		  double val = ykj - (yij /yii) *yki;
		  if (fabs(val) < Eta) val = 0.;
		  gsl_matrix_set (YY, k, j, val);
		  		  
		  // Fixup solution vectors from solution file
		  if (interactive)
		    for (int l = 0; l < npts; l++)
		      {
			double Yki = gsl_matrix_get (Y, l, dim1*i+k);
			double Ykj = gsl_matrix_get (Y, l, dim1*j+k);
			val = Ykj - (yij /yii) *Yki;
			if (fabs(val) < Eta) val = 0.;
			gsl_matrix_set (Y, l, dim1*j+k, val);
		      }
		}

	      // Fixup Psi and dPsi matrices
	      for (int l = 0; l < vac; l++)
		{
		  double val;
		  val = gsl_matrix_get (Psi,  l, j) - (yij /yii) *gsl_matrix_get (Psi,  l, i);
		  gsl_matrix_set (Psi,  l, j, val);
		  val = gsl_matrix_get (dPsi, l, j) - (yij /yii) *gsl_matrix_get (dPsi, l, i);
		  gsl_matrix_set (dPsi, l, j, val);
		}
	    }
	}
    }
  
  // Renormalize independent solution vectors launched from axis
  for (int i = 0; i < dim; i++)
    {
      double yii = (mpol[i] == 0) ? gsl_matrix_get (YY, dim+i, i)
	: gsl_matrix_get (YY, i, i);
      for (int k = 0; k < dim1; k++)
	{
	  double yki = gsl_matrix_get (YY, k, i);

	  // Renormalize solution vector
	  gsl_matrix_set (YY, k, i, yki /yii);
	
	  // Renormalize solution vectors from solution file
	  if (interactive)
	    for (int l = 0; l < npts; l++)
	      {
		double Yki = gsl_matrix_get (Y, l, dim1*i+k);
		gsl_matrix_set (Y, l, dim1*i+k, Yki /yii);
	      }
	}

      // Renormalize Psi and dPsi matrices
      for (int l = 0; l < vac; l++)
	{
	  double val;
	  val = gsl_matrix_get (Psi, l, i)  /yii;
	  gsl_matrix_set (Psi, l, i, val);
	  val = gsl_matrix_get (dPsi, l, i) /yii;
	  gsl_matrix_set (dPsi, l, i, val);
	}
    }
  
  // Rewrite solution file
  if (interactive)
    {
      soln = OpenFile (_soln);
      for (int l = 0; l < npts; l++)
	{
	  fprintf (soln, "%e", rr[l]);
	  for (int i = 0; i < dimN; i++)
	    fprintf (soln, " %20.15e", gsl_matrix_get (Y, l, i));
	  fprintf (soln, "\n");
	}
      fclose (soln);
    }

  if (interactive)
    {
      delete[] rr;
      gsl_matrix_free (Y);
    }
}

// Vacuum.cpp

#include "Tomuhawc.h"

// ############################################
// Function to calculate vacuum matching matrix
// ############################################
void Thawc::Vacuum (int interactive)
{
  // Calculate vacuum matrix
  Couple (1.);

  gsl_matrix      *Mat = gsl_matrix_alloc      (dim, dim);
  gsl_matrix      *Inv = gsl_matrix_alloc      (dim, dim);
  gsl_permutation *p   = gsl_permutation_alloc (dim);
  int s;

  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      gsl_matrix_set (Mat, i, j, gsl_matrix_get (Pmat, i, j));

  gsl_linalg_LU_decomp (Mat, p, &s);
  gsl_linalg_LU_invert (Mat, p, Inv);

  gsl_matrix_free      (Mat);
  gsl_permutation_free (p);

  double *rhs = new double [dim];
  double *xi  = new double [dim];

  for (int l = 0; l < dim; l++)
    {
      double fac = cos (double (ntor+1) * M_PI) 
	*gsl_sf_gamma (double (- mpol[l] - ntor) + 0.5)
	/gsl_sf_gamma (double (- mpol[l] + ntor) + 0.5);
      int    ml = abs (mpol[l]);
      double fl = (mpol[l] < 0) ? 1. : -1.;

      for (int j = 0; j < dim; j++)
	{
	  int    mj  = abs (mpol[j]);
	  double fcj = (mj == 0)      ? 1.  :  0.5;
	  double fsj = (mpol[j] > 0)  ? 0.5 : -0.5;

	  double sum = 0.;
	  for (int jj = 0; jj < dim; jj++)
	    {
	      int    mjj  = abs (mpol[jj]);
	      double fcjj = (mjj == 0)     ? 1.  :  0.5;
	      double fsjj = (mpol[jj] > 0) ? 0.5 : -0.5;

	      sum += gsl_matrix_get (Nmat, j, jj) 
		* (fcjj * gsl_matrix_get (Qc, mjj, ml) + fl*fsjj * gsl_matrix_get (Qs, mjj, ml));
	    }
	 
	  rhs[j] = (double (mpol[j]) - double (ntor) * qa) 
	    * (fcj * gsl_matrix_get (dQcdr, mj, ml) + fl*fsj * gsl_matrix_get (dQsdr, mj, ml)) - sum;
	}

      for (int j = 0; j < dim; j++)
	{
	  xi[j] = 0.;
	  
	  for (int jj = 0; jj < dim; jj++)
	    xi[j] += gsl_matrix_get (Inv, j, jj) * rhs[jj];
	}	

      for (int j = 0; j < dim; j++)
	{
	  int    mj  = abs (mpol[j]);
	  double fcj = (mj == 0)     ? 1.  :  0.5;
	  double fsj = (mpol[j] > 0) ? 0.5 : -0.5;

	  gsl_matrix_set (Vmat, j,     l, - fac * (fcj * gsl_matrix_get (Qc, mj, ml) + fl*fsj * gsl_matrix_get (Qs, mj, ml)));
	  gsl_matrix_set (Vmat, dim+j, l,   fac * xi[j]);
	}
    }	

  for (int l = 0; l < dim; l++)
    {
      double fac = - cos (double (ntor+1) * M_PI) 
	* gsl_sf_gamma (double (- mpol[l] - ntor) + 0.5)
	/ gsl_sf_gamma (double (- mpol[l] + ntor) + 0.5);
      int    ml = abs (mpol[l]);
      double fl = (mpol[l] < 0) ? 1. : -1.;

       for (int j = 0; j < dim; j++)
	{
	  int    mj  = abs (mpol[j]);
	  double fcj = (mj == 0)     ? 1.  :  0.5;
	  double fsj = (mpol[j] > 0) ? 0.5 : -0.5;

	  double sum = 0.;
	  for (int jj = 0; jj < dim; jj++)
	    {
	      int    mjj  = abs (mpol[jj]);
	      double fcjj = (mjj == 0)     ? 1.  :  0.5;
	      double fsjj = (mpol[jj] > 0) ? 0.5 : -0.5;

	      sum += gsl_matrix_get (Nmat, j, jj) 
		* (fcjj * gsl_matrix_get (Pc, mjj, ml) + fl*fsjj * gsl_matrix_get (Ps, mjj, ml));
	    }
	 
	  rhs[j] = (double (mpol[j]) - double (ntor) * qa) 
	    * (fcj * gsl_matrix_get (dPcdr, mj, ml) + fl*fsj * gsl_matrix_get (dPsdr, mj, ml)) - sum;
	}

      for (int j = 0; j < dim; j++)
	{
	  xi[j] = 0.;
	  
	  for (int jj = 0; jj < dim; jj++)
	    xi[j] += gsl_matrix_get (Inv, j, jj) * rhs[jj];
	}	

      for (int j = 0; j < dim; j++)
	{
	  int    mj  = abs (mpol[j]);
	  double fcj = (mj == 0)      ? 1.  :  0.5;
	  double fsj = (mpol[j]  > 0) ? 0.5 : -0.5;

	  gsl_matrix_set (Vmat, j,     dim+l, - fac * (fcj * gsl_matrix_get (Pc, mj, ml) + fl*fsj * gsl_matrix_get (Ps, mj, ml)));
	  gsl_matrix_set (Vmat, dim+j, dim+l,   fac * xi[j]);
	}
    }	

  gsl_matrix_free (Inv);
  delete[] rhs; delete[] xi;

  // Diagonalize vacuum matrix
  printf ("Vacuum matrix residual: %11.4e\n", VacuumResidual ());
  double res;
  do 
    {
      for (int i = 0; i < dim; i++)
	{  
	  double sumii = 0.;
	  for (int k = 0; k < dim; k++)
	    sumii += (+ gsl_matrix_get (Vmat, k,     i) * gsl_matrix_get (Vmat, dim+k, dim+i)
		      - gsl_matrix_get (Vmat, dim+k, i) * gsl_matrix_get (Vmat, k,     dim+i)) 
	      * (double (mpol[k]) - double (ntor) * qa);
	  
	  for (int j = 0; j < dim; j++)
	    if (j != i)
	      {
		double sumij = 0.;
		for (int k = 0; k < dim; k++)
		  sumij += (+ gsl_matrix_get (Vmat, k,     i) * gsl_matrix_get (Vmat, dim+k, dim+j)
			    - gsl_matrix_get (Vmat, dim+k, i) * gsl_matrix_get (Vmat, k,     dim+j)) 
		    * (double (mpol[k]) - double (ntor) * qa);
		
		for (int k = 0; k < dim1; k++)
		  {
		    double val = gsl_matrix_get (Vmat, k, dim+j) - (sumij/sumii) * gsl_matrix_get (Vmat, k, dim+i);
		    gsl_matrix_set (Vmat, k, dim+j, val);
		  }
	      }
	}
      for (int i = 0; i < dim; i++)
	{  
	  double sumii = 0.;
	  for (int k = 0; k < dim; k++)
	    sumii += (+ gsl_matrix_get (Vmat, k,     i) * gsl_matrix_get (Vmat, dim+k, dim+i)
		      - gsl_matrix_get (Vmat, dim+k, i) * gsl_matrix_get (Vmat, k,     dim+i)) 
	      * (double (mpol[k]) - double (ntor) * qa);
	  
	  for (int j = 0; j < dim; j++)
	    if (j != i)
	      {
		double sumij = 0.;
		for (int k = 0; k < dim; k++)
		  sumij += (+ gsl_matrix_get (Vmat, k,     i) * gsl_matrix_get (Vmat, dim+k, j)
			    - gsl_matrix_get (Vmat, dim+k, i) * gsl_matrix_get (Vmat, k,     j)) 
		    * (double (mpol[k]) - double (ntor) * qa);
		
		for (int k = 0; k < dim1; k++)
		  {
		    double val = gsl_matrix_get (Vmat, k, j) - (sumij/sumii) * gsl_matrix_get (Vmat, k, dim+i);
		    gsl_matrix_set (Vmat, k, j, val);
		  }
	      }
	}
      for (int i = 0; i < dim; i++)
	{  
	  double sumii = 0.;
	  for (int k = 0; k < dim; k++)
	    sumii += (+ gsl_matrix_get (Vmat, k,     dim+i) * gsl_matrix_get (Vmat, dim+k, i)
		      - gsl_matrix_get (Vmat, dim+k, dim+i) * gsl_matrix_get (Vmat, k,     i)) 
	      * (double (mpol[k]) - double (ntor) * qa);
	  
	  for (int j = 0; j < dim; j++)
	    if (j != i)
	      {
		double sumij = 0.;
		for (int k = 0; k < dim; k++)
		  sumij += (+ gsl_matrix_get (Vmat, k,     dim+i) * gsl_matrix_get (Vmat, dim+k, dim+j)
			    - gsl_matrix_get (Vmat, dim+k, dim+i) * gsl_matrix_get (Vmat, k,     dim+j)) 
		    * (double (mpol[k]) - double (ntor) * qa);
		
		for (int k = 0; k < dim1; k++)
		  {
		    double val = gsl_matrix_get (Vmat, k, dim+j) - (sumij/sumii) * gsl_matrix_get (Vmat, k, i);
		    gsl_matrix_set (Vmat, k, dim+j, val);
		  }
	      }
	}
      res =  VacuumResidual ();
      printf ("Vacuum matrix residual: %11.4e\n", res);
    
    } while (res > 1.e-6);
  
  // Log vacuum matrix
  if (interactive) 
    LogVmat ();
}

// ########################################################
// Function to calculate residual of vacuum matching matrix
// ########################################################
double Thawc::VacuumResidual ()
{  
  double residual = 0.;

  for (int i = 0; i < dim1; i++)
    for (int j = 0; j < dim1; j++)
      {
	double sum = 0.;

	for (int k = 0; k < dim; k++)
	  sum += (+ gsl_matrix_get (Vmat, k,     i) * gsl_matrix_get (Vmat, dim+k, j)
		  - gsl_matrix_get (Vmat, dim+k, i) * gsl_matrix_get (Vmat, k,     j)) 
	  * (double (mpol[k]) - double (ntor) * qa);
	
	if (i - j == dim || j - i == dim)
	  sum = 0.;

	if (fabs (sum) > residual)
	  residual = fabs (sum);
      }

  return residual;
}

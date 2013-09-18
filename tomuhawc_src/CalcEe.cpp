// CalcEe.cpp

#include "Tomuhawc.h"

// #######################################
// Function to calculate Eo and G matrices
// #######################################
void Thawc::CalcEe (gsl_matrix *Fee, gsl_matrix *Foe, gsl_matrix *Ee, gsl_matrix *Gp)
{
  // Calculate Ee matrix
  gsl_matrix      *Fmat = gsl_matrix_alloc      (num, num);
  gsl_permutation *p    = gsl_permutation_alloc (num);

  for (int i = 0; i < num; i++)
    for (int k = 0; k < num; k++)
      gsl_matrix_set (Fmat, i, k, gsl_matrix_get (Fee, i, k));

  int s;
  gsl_linalg_LU_decomp (Fmat, p, &s);
  gsl_linalg_LU_invert (Fmat, p, Ee);

  // Calculate Gp matrix
  for (int i = 0; i < num; i++)
    for (int j = 0; j < num; j++)
      {
	double sum = 0.;
	for (int k = 0; k < num; k++)
	  sum += gsl_matrix_get (Foe, i, k) * gsl_matrix_get (Ee, k, j);
	
	gsl_matrix_set (Gp, i, j, sum);
      }

  gsl_matrix_free      (Fmat);
  gsl_permutation_free (p); 
}

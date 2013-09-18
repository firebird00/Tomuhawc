// CalcEo.cpp

#include "Tomuhawc.h"

// #######################################
// Function to calculate Eo and G matrices
// #######################################
void Thawc::CalcEo (gsl_matrix *Foo, gsl_matrix *Feo, gsl_matrix *Eo, gsl_matrix *G)
{
  // Calculate Eo matrix
  gsl_matrix      *Fmat = gsl_matrix_alloc      (num, num);
  gsl_permutation *p    = gsl_permutation_alloc (num);

  for (int i = 0; i < num; i++)
    for (int k = 0; k < num; k++)
      gsl_matrix_set (Fmat, i, k, gsl_matrix_get (Foo, i, k));

  int s;
  gsl_linalg_LU_decomp (Fmat, p, &s);
  gsl_linalg_LU_invert (Fmat, p, Eo);

  // Calculate G matrix
  for (int i = 0; i < num; i++)
    for (int j = 0; j < num; j++)
      {
	double sum = 0.;
	for (int k = 0; k < num; k++)
	  sum += gsl_matrix_get (Feo, i, k) * gsl_matrix_get (Eo, k, j);
	
	gsl_matrix_set (G, i, j, sum);
      }

  gsl_matrix_free      (Fmat);
  gsl_permutation_free (p); 
}


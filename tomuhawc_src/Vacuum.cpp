// Vacuum.cpp

#include "Tomuhawc.h"

// ############################################
// Function to calculate vacuum matching matrix
// ############################################
void Thawc::Vacuum ()
{
  gsl_matrix      *Mat = gsl_matrix_alloc      (dim1, dim1);
  gsl_permutation *p   = gsl_permutation_alloc (dim1);
  int s;

  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      {
	int mi = abs(mpol[i]);
	int mj = abs(mpol[j]);						
	//gsl_matrix_set (Mat, i,     j,     gsl_matrix_get (Ptor0, mi, mj));
	//gsl_matrix_set (Mat, i+dim, j,     gsl_matrix_get (Ptor1, mi, mj));
	//gsl_matrix_set (Mat, i,     j+dim, gsl_matrix_get (Qtor0, mi, mj));
	//gsl_matrix_set (Mat, i+dim, j+dim, gsl_matrix_get (Qtor1, mi, mj));
      }

  //  gsl_linalg_LU_decomp (Mat, p, &s);
  //gsl_linalg_LU_invert (Mat, p, Vmat);

  gsl_matrix_free      (Mat);
  gsl_permutation_free (p);
}

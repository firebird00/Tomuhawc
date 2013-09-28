// Solv.cpp

#include "Tomuhawc.h"

// ####################################################################################
// Function to calculate edge data
//
// YY(dim1, dim+vac)     .. solution vector at r = a
//  YY(i=0,  dim -1;k=0,dim-1)       ... psi: i     .. poloidal harmonic index
//  YY(i=dim,dim1-1;k=0,dim-1)       ... Z  : i-dim .. poloidal harmonic index
//                                            k     .. index of independent solutions
//                                                      launched from axis
//  YY(i=0,  dim -1;j=dim,dim+vac-1) ... psi: i     .. poloidal harmonic index
//  YY(i=dim,dim1-1;j=dim,dim+vac-1) ... Z  : i-dim .. poloidal harmonic index
//                                            j     .. index of small solutions
//                                                      launched from rational surfaces 
// Ya(dim1, dim+vac)     .. solution vector at r = a
//  Ya(i=0,  dim -1;k=0,dim-1)       ... psi: i     .. poloidal harmonic index
//  Ya(i=dim,dim1-1;k=0,dim-1)       ... Z  : i-dim .. poloidal harmonic index
//                                            k     .. index of independent solutions
//                                                      launched from axis
//  Ya(i=0,  dim -1;j=dim,dim+vac-1) ... psi: i     .. poloidal harmonic index
//  Ya(i=dim,dim1-1;j=dim,dim+vac-1) ... Z  : i-dim .. poloidal harmonic index
//                                            j     .. index of small solutions
//                                                      launched from rational surfaces 
// QP(dim1, dim+vac)     .. q, p coefficients of solution vector at r = a
//  QP(i=0,  dim -1;k=0,dim-1)       ... q  : i     .. poloidal harmonic index
//  QP(i=dim,dim1-1;k=0,dim-1)       ... p  : i-dim .. poloidal harmonic index
//                                            k     .. index of independent solutions
//                                                      launched from axis
//  QP(i=0,  dim -1;j=dim,dim+vac-1) ... q  : i     .. poloidal harmonic index
//  QP(i=dim,dim1-1;j=dim,dim+vac-1) ... p  : i-dim .. poloidal harmonic index
//                                            j     .. index of small solutions
//                                                      launched from rational surfaces 
// ####################################################################################
void Thawc::Edge (gsl_matrix *YY, gsl_matrix *Ya, gsl_matrix *QP, int interactive)
{
  // Calculate Ya
  for (int i = 0; i < dim1; i++)
    for (int j = 0; j < dim+vac; j++)
      gsl_matrix_set (Ya, i, j, gsl_matrix_get (YY, i, j));

  // Calculate QP
  for (int i = 0; i < dim+vac; i++)
    {
      for (int l = 0; l < dim; l++)
	{
	  double p = 0, q = 0.;
	  
	  for (int j = 0; j < dim; j++)
	    {
	      p +=
		+ gsl_matrix_get (Vmat, j,     l)     * gsl_matrix_get (Ya, j,     i) 
		+ gsl_matrix_get (Vmat, dim+j, l)     * gsl_matrix_get (Ya, dim+j, i);
	      q +=
		+ gsl_matrix_get (Vmat, j,     dim+l) * gsl_matrix_get (Ya, j,     i) 
		+ gsl_matrix_get (Vmat, dim+j, dim+l) * gsl_matrix_get (Ya, dim+j, i);
	    }

	  gsl_matrix_set (QP, l,     i, q);
	  gsl_matrix_set (QP, dim+l, i, p);
	}
    }

  if (interactive)
    {
      FILE *file = OpenFile ("Stage3/Edge.out");
      for (int i = 0; i < dim; i++)
	for (int j = 0; j < dim; j++)
	  fprintf (file, "%3d %3d %11.4e %11.4e %11.4e %11.4e\n",
		  mpol[i], mpol[j], gsl_matrix_get (Ya, j, i), gsl_matrix_get (Ya, dim+j, i),
		   gsl_matrix_get (QP, j, i), gsl_matrix_get (QP, dim+j, i));
      for (int i = 0; i < vac; i++)
	for (int j = 0; j < dim; j++)
	  fprintf (file, "%3d %3d %11.4e %11.4e %11.4e %11.4e\n",
		   Mres[i], mpol[j], gsl_matrix_get (Ya, j, dim+i), gsl_matrix_get (Ya, dim+j, dim+i),
		   gsl_matrix_get (QP, j, dim+i), gsl_matrix_get (QP, dim+j, dim+i));
      fclose (file);
    }
}

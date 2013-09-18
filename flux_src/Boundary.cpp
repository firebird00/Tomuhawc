// Boundary.cpp

#include "Flux.h"

// ##################################
// Routine to calculate boundary data
// ##################################
void Flux::Boundary ()
{
  gsl_matrix *PPa = gsl_matrix_alloc (K/2, K/2);
  gsl_matrix *QQa = gsl_matrix_alloc (K/2, K/2);
  gsl_matrix *PPb = gsl_matrix_alloc (K/2, K/2);
  gsl_matrix *QQb = gsl_matrix_alloc (K/2, K/2);

  double *Pa = new double[K];
  double *Qa = new double[K];
  double *Pb = new double[K];
  double *Qb = new double[K];

  double *Pm = new double[K/2];
  double *Qm = new double[K/2];

  for (int m = 0; m < K/2; m++)
    {
      for (int k = 0; k < K; k++)
	{
	  double rra = 1.;
	  double RRa = Extrapolate2 (J, rP, Rst, rra, k);
	  double ZZa = Extrapolate2 (J, rP, Zst, rra, k);
	  double Ra  = Raxis + rra * RRa;
	  double Za  =         rra * ZZa;
	  double d1a = sqrt((Ra + 1.) * (Ra + 1.) + Za*Za);
	  double d2a = sqrt((Ra - 1.) * (Ra - 1.) + Za*Za);
	  double cma = 0.5 * (d1a/d2a + d2a/d1a);
	  double eta; 
	  if (k == 0)
	    eta = M_PI;
	  else if (k == K-1)
	    eta = 0.;
	  else
	    eta = acos ((d1a*d1a + d2a*d2a - 4.) /2./d1a/d2a);
	  
	  Pa[k] = sqrt(cma - cos(eta)) * ToroidalP (ntor, m, cma) * cos (double (m) * eta);
	  Qa[k] = sqrt(cma - cos(eta)) * ToroidalQ (ntor, m, cma) * cos (double (m) * eta);
	    
	  double rrb = rP[J-1];
	  double RRb = gsl_matrix_get (Rst, J-1, k);
	  double ZZb = gsl_matrix_get (Zst, J-1, k);
	  double Rb  = Raxis + rrb * RRb;
	  double Zb  =         rrb * ZZb;
	  double d1b = sqrt((Rb + 1.) * (Rb + 1.) + Zb*Zb);
	  double d2b = sqrt((Rb - 1.) * (Rb - 1.) + Zb*Zb);
	  double cmb = 0.5 * (d1b/d2b + d2b/d1b);
	  double etb;
	  if (k == 0)
	    etb = M_PI;
	  else if (k == K-1)
	    etb = 0.;
	  else
	    etb = acos ((d1b*d1b + d2b*d2b - 4.) /2./d1b/d2b);

	  Pb[k] = sqrt(cmb - cos(etb)) * ToroidalP (ntor, m, cmb) * cos (double (m) * etb);
	  Qb[k] = sqrt(cmb - cos(etb)) * ToroidalQ (ntor, m, cmb) * cos (double (m) * etb);
	}
      
      Pm[m] = -1.; Qm[m] = -1.;
      for (int k = 0; k < K; k++)
	{
	  if (fabs (Pb[k]) > Pm[m]) Pm[m] = fabs (Pb[k]); 
	  if (fabs (Qa[k]) > Qm[m]) Qm[m] = fabs (Qa[k]); 
	}

      for (int k = 0; k < K; k++)
	{
	  Pa[k] /= Pm[m];
	  Qa[k] /= Qm[m];
	  Pb[k] /= Pm[m];
	  Qb[k] /= Qm[m];
	}	
      
      FTCos (m, Pa, PPa);
      FTCos (m, Qa, QQa);
      FTCos (m, Pb, PPb);
      FTCos (m, Qb, QQb);
    }

  FILE *file1 = OpenFile ("Stage2/P.out");
  FILE *file2 = OpenFile ("Stage2/Q.out");
  FILE *file3 = OpenFile ("Stage2/PP.out");
  FILE *file4 = OpenFile ("Stage2/QQ.out");
  FILE *file5 = OpenFile ("Stage2/PQm.out");
  for (int m = 0; m < K/2; m++)
    {
      for (int mm = 0; mm < K/2; mm++)
	{
	  double P  = gsl_matrix_get (PPa, mm, m);
	  double Q  = gsl_matrix_get (QQa, mm, m);
	  double PP = (P - gsl_matrix_get (PPb, mm, m)) /(1. - rP[J-1]);
	  double QQ = (Q - gsl_matrix_get (QQb, mm, m)) /(1. - rP[J-1]);

	  fprintf (file1, "%e ", P);
	  fprintf (file2, "%e ", Q);
	  fprintf (file3, "%e ", PP);
	  fprintf (file4, "%e ", QQ);
	}
      fprintf (file1, "\n"); fprintf (file2, "\n"); 
      fprintf (file3, "\n"); fprintf (file4, "\n"); 

      fprintf (file5, "%e %e\n", Pm[m], Qm[m]);
    }
  fclose (file1); fclose (file2); fclose (file3); fclose (file4); fclose (file5);

  delete[] Pa; delete[] Qa; delete[] Pb; delete[] Qb;
  delete[] Pm; delete[] Qm;
  gsl_matrix_free (PPa); gsl_matrix_free (QQa); gsl_matrix_free (PPb); gsl_matrix_free (QQb);

}

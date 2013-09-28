// Boundary.cpp

#include "Flux.h"

// ##################################
// Routine to calculate boundary data
// ##################################
void Flux::Boundary ()
{
  gsl_matrix *PPca = gsl_matrix_alloc (K/2, K/2);
  gsl_matrix *QQca = gsl_matrix_alloc (K/2, K/2);
  gsl_matrix *PPsa = gsl_matrix_alloc (K/2, K/2);
  gsl_matrix *QQsa = gsl_matrix_alloc (K/2, K/2);
  gsl_matrix *PPcb = gsl_matrix_alloc (K/2, K/2);
  gsl_matrix *QQcb = gsl_matrix_alloc (K/2, K/2);
  gsl_matrix *PPsb = gsl_matrix_alloc (K/2, K/2);
  gsl_matrix *QQsb = gsl_matrix_alloc (K/2, K/2);

  double *Pca = new double[K];
  double *Qca = new double[K];
  double *Psa = new double[K];
  double *Qsa = new double[K];
  double *Pcb = new double[K];
  double *Qcb = new double[K];
  double *Psb = new double[K];
  double *Qsb = new double[K];

  double dr = 1.e-6;
  for (int m = 0; m < K/2; m++)
    {
      for (int k = 0; k < K; k++)
	{
	  double rra = 1.;
	  double RRa = Extrapolate2 (J, rP, Rst, rra, k);
	  double ZZa = Extrapolate2 (J, rP, Zst, rra, k);
	  double Ra  = Raxis + rra * RRa;
	  double Za  =         rra * ZZa;
	  double d1a = sqrt((Ra + 1.) * (Ra + 1.) + Za * Za);
	  double d2a = sqrt((Ra - 1.) * (Ra - 1.) + Za * Za);
	  double cma = 0.5 * (d1a/d2a + d2a/d1a);
	  double sma = sqrt (cma*cma - 1.);
	  double eta; 
	  if (k == 0)
	    eta = M_PI;
	  else if (k == K-1)
	    eta = 0.;
	  else
	    eta = acos ((d1a*d1a + d2a*d2a - 4.) /2./d1a/d2a);
	  
	  Pca[k] = sqrt(cma - cos(eta)) * ToroidalP (ntor, m, cma) * cos (double (m) * eta);
	  Qca[k] = sqrt(cma - cos(eta)) * ToroidalQ (ntor, m, cma) * cos (double (m) * eta);
	  Psa[k] = sqrt(cma - cos(eta)) * ToroidalP (ntor, m, cma) * sin (double (m) * eta);
	  Qsa[k] = sqrt(cma - cos(eta)) * ToroidalQ (ntor, m, cma) * sin (double (m) * eta);
	    
	  double rrb = 1. - dr;
	  double RRb = Extrapolate2 (J, rP, Rst, rrb, k);
	  double ZZb = Extrapolate2 (J, rP, Zst, rrb, k);
	  double Rb  = Raxis + rrb * RRb;
	  double Zb  =         rrb * ZZb;
	  double d1b = sqrt((Rb + 1.) * (Rb + 1.) + Zb * Zb);
	  double d2b = sqrt((Rb - 1.) * (Rb - 1.) + Zb * Zb);
	  double cmb = 0.5 * (d1b/d2b + d2b/d1b);
	  double etb;
	  if (k == 0)
	    etb = M_PI;
	  else if (k == K-1)
	    etb = 0.;
	  else
	    etb = acos ((d1b*d1b + d2b*d2b - 4.) /2./d1b/d2b);

	  double dcmdr = (cma - cmb) / dr;
	  double detdr = (eta - etb) / dr;

	  Pcb[k] = 
	    + (0.5/sqrt(cma - cos(eta))) * (dcmdr + sin(eta)*detdr) 
	    * ToroidalP (ntor, m, cma) * cos (double (m) * eta) 
	    + sqrt(cma - cos(eta)) * ToroidalPp (ntor, m, cma) * dcmdr * cos (double(m) * eta)
	    - sqrt(cma - cos(eta)) * ToroidalP (ntor, m, cma) * double (m) * sin (double (m) * eta) * detdr;
	  Qcb[k] =
	    + (0.5/sqrt(cma - cos(eta))) * (dcmdr + sin(eta)*detdr) 
	    * ToroidalQ (ntor, m, cma) * cos (double (m) * eta) 
	    + sqrt(cma - cos(eta)) * ToroidalQp (ntor, m, cma) * dcmdr * cos (double (m) * eta)
	    - sqrt(cma - cos(eta)) * ToroidalQ (ntor, m, cma) * double (m) * sin (double (m) * eta) * detdr;
	  Psb[k] = 
	    + (0.5/sqrt(cma - cos(eta))) * (dcmdr + sin(eta)*detdr) 
	    * ToroidalP (ntor, m, cma) * sin (double (m) * eta) 
	    + sqrt(cma - cos(eta)) * ToroidalPp (ntor, m, cma) * dcmdr * sin (double(m) * eta)
	    + sqrt(cma - cos(eta)) * ToroidalP (ntor, m, cma) * double (m) * cos (double (m) * eta) * detdr;
	  Qsb[k] =
	    + (0.5/sqrt(cma - cos(eta))) * (dcmdr + sin(eta)*detdr) 
	    * ToroidalQ (ntor, m, cma) * sin (double (m) * eta) 
	    + sqrt(cma - cos(eta)) * ToroidalQp (ntor, m, cma) * dcmdr * sin (double (m) * eta)
	    + sqrt(cma - cos(eta)) * ToroidalQ (ntor, m, cma) * double (m) * cos (double (m) * eta) * detdr;
	}
      
      double Pm = 0., Qm = 0.;
      for (int k = 0; k < K; k++)
	{
	  if (fabs (Pca[k]) > Pm) Pm = fabs (Pca[k]); 
	  if (fabs (Qca[k]) > Qm) Qm = fabs (Qca[k]); 
	}

      for (int k = 0; k < K; k++)
	{
	  Pca[k] /= sqrt(Pm/Qm);
	  Qca[k] *= sqrt(Pm/Qm);
	  Psa[k] /= sqrt(Pm/Qm);
	  Qsa[k] *= sqrt(Pm/Qm);
	  Pcb[k] /= sqrt(Pm/Qm);
	  Qcb[k] *= sqrt(Pm/Qm);
	  Psb[k] /= sqrt(Pm/Qm);
	  Qsb[k] *= sqrt(Pm/Qm);
	}	
      
      FTCos (m, Pca, PPca);
      FTCos (m, Qca, QQca);
      FTCos (m, Pcb, PPcb);
      FTCos (m, Qcb, QQcb);
      FTSin (m, Psa, PPsa);
      FTSin (m, Qsa, QQsa);
      FTSin (m, Psb, PPsb);
      FTSin (m, Qsb, QQsb);
    }

  FILE *file1 = OpenFile ("Stage2/Pc.out");
  FILE *file2 = OpenFile ("Stage2/Qc.out");
  FILE *file3 = OpenFile ("Stage2/dPcdr.out");
  FILE *file4 = OpenFile ("Stage2/dQcdr.out");
  FILE *file5 = OpenFile ("Stage2/Ps.out");
  FILE *file6 = OpenFile ("Stage2/Qs.out");
  FILE *file7 = OpenFile ("Stage2/dPsdr.out");
  FILE *file8 = OpenFile ("Stage2/dQsdr.out");
  for (int m = 0; m < K/2; m++)
    {
      for (int mm = 0; mm < K/2; mm++)
	{
	  double Pc    = gsl_matrix_get (PPca, mm, m);
	  double Qc    = gsl_matrix_get (QQca, mm, m);
	  double dPcdr = gsl_matrix_get (PPcb, mm, m);
	  double dQcdr = gsl_matrix_get (QQcb, mm, m);
	  double Ps    = gsl_matrix_get (PPsa, mm, m);
	  double Qs    = gsl_matrix_get (QQsa, mm, m);
	  double dPsdr = gsl_matrix_get (PPsb, mm, m);
	  double dQsdr = gsl_matrix_get (QQsb, mm, m);

	  fprintf (file1, "%17.10e ", Pc);
	  fprintf (file2, "%17.10e ", Qc);
	  fprintf (file3, "%17.10e ", dPcdr);
	  fprintf (file4, "%17.10e ", dQcdr);
	  fprintf (file5, "%17.10e ", Ps);
	  fprintf (file6, "%17.10e ", Qs);
	  fprintf (file7, "%17.10e ", dPsdr);
	  fprintf (file8, "%17.10e ", dQsdr);
	}
      fprintf (file1, "\n"); fprintf (file2, "\n"); 
      fprintf (file3, "\n"); fprintf (file4, "\n"); 
      fprintf (file5, "\n"); fprintf (file6, "\n"); 
      fprintf (file7, "\n"); fprintf (file8, "\n"); 
    }
  fclose (file1); fclose (file2); fclose (file3); fclose (file4); 
  fclose (file5); fclose (file6); fclose (file7); fclose (file8); 

  delete[] Pca; delete[] Qca; delete[] Pcb; delete[] Qcb;
  delete[] Psa; delete[] Qsa; delete[] Psb; delete[] Qsb;
  gsl_matrix_free (PPca); gsl_matrix_free (QQca); gsl_matrix_free (PPcb); gsl_matrix_free (QQcb);
  gsl_matrix_free (PPsa); gsl_matrix_free (QQsa); gsl_matrix_free (PPsb); gsl_matrix_free (QQsb);
}

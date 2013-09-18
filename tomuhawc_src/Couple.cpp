// Couple.cpp

#include "Tomuhawc.h"

// ###################################################
// Function to calculate elements of coupling matrices
//
// R ... flux surface label
//
// ###################################################
void Thawc::Couple (double R)
{
  // Get profile data
  double Q   = gsl_spline_eval (Sq,      R, Aq);
  double QP  = gsl_spline_eval (Sqp,     R, Aqp);
  double P1  = gsl_spline_eval (Sprof1,  R, Aprof1);
  double P2  = gsl_spline_eval (Sprof2,  R, Aprof2);
  double P3  = gsl_spline_eval (Sprof3,  R, Aprof3);
  double P4  = gsl_spline_eval (Sprof4,  R, Aprof4);
  double P1P = gsl_spline_eval (Sprof1p, R, Aprof1p);
  double P2P = gsl_spline_eval (Sprof2p, R, Aprof2p);
  double P3P = gsl_spline_eval (Sprof3p, R, Aprof3p);

  // Get metric data
  double *MET1  = new double [diag+1];
  double *MET2  = new double [diag+1];
  double *MET3  = new double [diag+1];
  double *MET4  = new double [diag+1];
  double *MET5  = new double [diag+1];
  double *MET6  = new double [diag+1];
  double *MET7  = new double [diag+1];
  double *MET1P = new double [diag+1];
  double MET3P, MET4P;
  for (int j = 0; j <= diag; j++)
    {
      MET1 [j] = gsl_spline_eval (SM1 [j], R, AM1 [j]);
      MET2 [j] = gsl_spline_eval (SM2 [j], R, AM2 [j]);
      MET3 [j] = gsl_spline_eval (SM3 [j], R, AM3 [j]);
      MET4 [j] = gsl_spline_eval (SM4 [j], R, AM4 [j]);
      MET5 [j] = gsl_spline_eval (SM5 [j], R, AM5 [j]);
      MET6 [j] = gsl_spline_eval (SM6 [j], R, AM6 [j]);
      MET7 [j] = gsl_spline_eval (SM7 [j], R, AM7 [j]);
      MET1P[j] = gsl_spline_eval (SM1P[j], R, AM1P[j]);
    }	
  MET3P = gsl_spline_eval (SM3P, R, AM3P);
  MET4P = gsl_spline_eval (SM4P, R, AM4P);

  // Calculate common quantities
  double N  = double(ntor);
  double NQ = N*Q;
  double S  = R*QP/Q;

  // Evaluate element of B, C, D, and E matrices
  double *lam  = new double[dim]; 
  double *lamp = new double[dim];
  double Bij, Cij, Dij, Eij;
  for (int i = 0; i < dim; i++)
    {
      int m1 = mpol[i];

      double M1   = double(m1);
      double M1NQ = M1 - NQ;

      for (int j = 0; j < dim; j++)
	{
	  int m2 = mpol[j];
	  
	  double M2   = double(m2);
	  double M2NQ = M2 - NQ;
	  
	  int k = m1 - m2;
	  if (k == 0)
	    {
	      Bij = M1*M2*MET3[k] + N*N*P1;
	      Cij =   M1*M2NQ*P2*MET3[k]/N - M1*P3*MET4[k];
	      Dij = - M1NQ*M2*P2*MET3[k]/N + M2*P3*MET4[k];
	      Eij = M1NQ*M2NQ*(MET2[k] - P2*P2*MET3[k]/N/N)
		- M1NQ*R*P2P/N
		+ P3*(P2*(M1NQ+M2NQ)*MET4[k]/N + R*MET1P[k] - P4*MET1[k] - P3*MET5[k]); 
	      lam[i]  = - Cij/Bij;
	      lamp[i] = 
		- (- M1*QP*P2*MET3[k] + M1*M1NQ*(P2P*MET3[k] + P2*MET3P)/N - M1*(P3P*MET4[k] + P3*MET4P)) 
		/ (M1*M2*MET3[k] + N*N*P1)
		- lam[i] * (M1*M2*MET3P + N*N*P1P) 
		/ (M1*M2*MET3[k] + N*N*P1);
	    }
	  else if (k > 0 && k <= diag)
	    {
	      Bij = M1*M2*0.5*MET3[k];
	      Cij =   M1*M2NQ*(- 0.5*MET6[k]  + P2*0.5*MET3[k]/N)  - M1*P3*0.5*MET4[k];
	      Dij = - M1NQ*M2*(+ 0.5*MET6[k]  + P2*0.5*MET3[k]/N)  + M2*P3*0.5*MET4[k];
	      Eij = M1NQ*M2NQ*(0.5*MET2[k] - P2*P2*0.5*MET3[k]/N/N)
		+ P3*(+ P2*(M1NQ+M2NQ)*0.5*MET4[k]/N - P3*0.5*MET5[k] + (M1-M2)*0.5*MET7[k] 
		      + R*0.5*MET1P[k] - P4*0.5*MET1[k]);
	    }	
	  else if (k < 0 && k >= -diag)
	    {
	      Bij = M1*M2*0.5*MET3[-k];
	      Cij =   M1*M2NQ*(+ 0.5*MET6[-k] + P2*0.5*MET3[-k]/N) - M1*P3*0.5*MET4[-k];
	      Dij = - M1NQ*M2*(- 0.5*MET6[-k] + P2*0.5*MET3[-k]/N) + M2*P3*0.5*MET4[-k];
	      Eij = M1NQ*M2NQ*(0.5*MET2[-k] - P2*P2*0.5*MET3[-k]/N/N)
		+ P3*(+ P2*(M1NQ+M2NQ)*0.5*MET4[-k]/N - P3*0.5*MET5[-k] - (M1-M2)*0.5*MET7[-k] 
		      + R*0.5*MET1P[-k] - P4*0.5*MET1[-k]);
	    }
	  else
	    {
	      Bij = 0.; Cij = 0.; Dij = 0.; Eij = 0.;
	    }

	  gsl_matrix_set (bmat, i, j, Bij);
	  gsl_matrix_set (cmat, i, j, Cij);
	  gsl_matrix_set (dmat, i, j, Dij);
	  gsl_matrix_set (emat, i, j, Eij);
	}
    }

  // Evaluate element of L, M, N, and P matrices
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      {
	double Bij = gsl_matrix_get (bmat, i, j);
	double Cij = gsl_matrix_get (cmat, i, j);
	double Dij = gsl_matrix_get (dmat, i, j);
	double Eij = gsl_matrix_get (emat, i, j);

	double Lij = Bij;
	double Mij = Cij + Lij * lam[j];
	double Nij = Dij - Lij * lam[i];
	double Pij = Eij - lam[i]*Mij + lam[j]*Nij + lam[i]*lam[j]*Lij;

	if (i == j)
	  {
	    double M1   = double(mpol[i]);
	    double M1NQ = M1 - NQ;
	    Pij += - lam[i]*NQ*S - M1NQ*R*lamp[i];
	  }	

	gsl_matrix_set (Lmat, i, j, Lij);
	gsl_matrix_set (Mmat, i, j, Mij);
	gsl_matrix_set (Nmat, i, j, Nij);
	gsl_matrix_set (Pmat, i, j, Pij);
      }
  
  // Clean up
  delete[] lam;  delete[] lamp;
  delete[] MET1; delete[] MET2; delete[] MET3; delete[] MET4;
  delete[] MET5; delete[] MET6; delete[] MET7; delete[] MET1P;
}



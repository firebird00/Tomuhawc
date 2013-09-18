// L1kP1k.cpp

#include "Tomuhawc.h"

// #####################################################################
// Function to calculate L1k, M1k, N1k, P1k coefficients at given radius
//
//  R  - flux surface label 
//  M  - resonant poloidal mode number
//  km - index of resonant harmonic
//
// #####################################################################
void Thawc::L1kP1k (double R, double M, int km, double *L1k, double *M1k, double *N1k, double *P1k)
{
  double Q0, Q1, Q2, QP, QPP, S;
  Q0 = gsl_spline_eval (Sq, R,     Aq);
  Q1 = gsl_spline_eval (Sq, R+eps, Aq);
  Q2 = gsl_spline_eval (Sq, R-eps, Aq);
  
  QP  = (Q1 - Q2) /(2.*eps);
  QPP = (Q1 - 2.*Q0 + Q2) /(eps*eps);
  S   = R*QP/Q0;

  double LM, MM, NM, PM, LMP, MMP, NMP, PMP;
  for (int k = 0; k < dim; k++)
    {
      Couple (R);
      LM  = gsl_matrix_get (Lmat, k, km);
      MM  = gsl_matrix_get (Mmat, k, km);
      NM  = gsl_matrix_get (Nmat, k, km);
      PM  = gsl_matrix_get (Pmat, k, km);

      Couple (R+eps);
      LMP = gsl_matrix_get (Lmat, k, km);
      MMP = gsl_matrix_get (Mmat, k, km);
      NMP = gsl_matrix_get (Nmat, k, km);
      PMP = gsl_matrix_get (Pmat, k, km);

      Couple (R-eps);
      LMP -= gsl_matrix_get (Lmat, k, km);
      MMP -= gsl_matrix_get (Mmat, k, km);
      NMP -= gsl_matrix_get (Nmat, k, km);
      PMP -= gsl_matrix_get (Pmat, k, km);
    
      LMP /= (2.*eps);
      MMP /= (2.*eps);
      NMP /= (2.*eps);
      PMP /= (2.*eps);
      
      L1k[k] = (R/M/S) * (0.5*(QPP/QP) * LM - LMP);
      M1k[k] = (R/M/S) * (0.5*(QPP/QP) * MM - MMP);
      N1k[k] = (R/M/S) * (0.5*(QPP/QP) * NM - NMP);
      P1k[k] = (R/M/S) * (0.5*(QPP/QP) * PM - PMP);
    }
}

// L1P1T1.cpp

#include "Tomuhawc.h"

// #################################################################
// Function to calculate coefficients L1, P1, and T1 at given radius
//
//  R  - flux surface label
//  M  - resonant poloidal mode number
//  km - index of resonant harmonic
//
// #################################################################
void Thawc::L1P1T1 (double R, double M, int km, double& L1, double& P1, double& T1, double& S)
{
  double Q0, Q1, Q2, QP, QPP, LM, LMP, PM, PMP;
  Q0  = gsl_spline_eval (Sq, R,     Aq);
  Q1  = gsl_spline_eval (Sq, R+eps, Aq);
  Q2  = gsl_spline_eval (Sq, R-eps, Aq);
  QP  = (Q1 - Q2) /(2.*eps);
  QPP = (Q1 - 2.*Q0 + Q2) /(eps*eps);
  S   = R*QP/Q0; 
  T1  = 1. + 0.5*R*QPP/QP;

  Couple (R);
  LM  = gsl_matrix_get (Lmat, km, km);
  PM  = gsl_matrix_get (Pmat, km, km);

  Couple (R+eps);
  LMP = gsl_matrix_get (Lmat, km, km);
  PMP = gsl_matrix_get (Pmat, km, km);

  Couple (R-eps);
  LMP -= gsl_matrix_get (Lmat, km, km);
  PMP -= gsl_matrix_get (Pmat, km, km);
  
  LMP /= (2.*eps);
  PMP /= (2.*eps);
  
  L1 = (R/M/S) * (0.5*(QPP/QP) * LM - LMP);
  P1 = (R/M/S) * (0.5*(QPP/QP) * PM - PMP);
}

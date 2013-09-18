// Misc.cpp

#include "Flux.h"

// ##############################
// Function to evaluate Psi (R, Z)
// ##############################
double Flux::GetPsi (double r, double z)
{
  return InterpolatePsi (r, z, 0);
}

// ###################################
// Function to evaluate dPsi/dR (R, Z)
// ###################################
double Flux::GetPsiR (double r, double z)
{
  return InterpolatePsi (r, z, 1);
}

// ###################################
// Function to evaluate dPsi/dZ (R, Z)
// ###################################
double Flux::GetPsiZ (double r, double z)
{
  return InterpolatePsi (r, z, 2);
}

// #######################################
// Function to evaluate d^2Psi/dR^2 (R, Z)
// #######################################
double Flux::GetPsiRR (double r, double z)
{
  return InterpolatePsi (r, z, 3);
}

// #######################################
// Function to evaluate d^2Psi/dZ^2 (R, Z)
// #######################################
double Flux::GetPsiZZ (double r, double z)
{
  return InterpolatePsi (r, z, 4);
}

// #####################################################
// Function to evaluate right-hand sides of q/g equation
// #####################################################
int Flux::Rhs1 (double r, const double y[], double dydr[], void *)
{
  // y[0] - R     
  // y[1] - Z
  // y[2] - q/g

  double PsiR = GetPsiR (y[0], y[1]);
  double PsiZ = GetPsiZ (y[0], y[1]);
  double Grad = sqrt (PsiR*PsiR + PsiZ*PsiZ);

  dydr[0] = - PsiZ /Grad;
  dydr[1] = + PsiR /Grad;
  dydr[2] = 1./M_PI /fabs(psic) /y[0] /Grad;

  return GSL_SUCCESS;
}

// ###################################################
// Function to evaluate right-hand sides of r equation
// ###################################################
int Flux::Rhs2 (double r, const double y[], double dydr[], void *)
{
  // y[0] - (epsa rP)^2    

  double qg = Interpolate (J, S, QGP, sqrt(1.-r), 0);

  dydr[0] = 2.*fabs(psic) * qg;

  return GSL_SUCCESS;
}

// #######################################################
// Function to evaluate right-hand sides of theta equation
// #######################################################
int Flux::Rhs3 (double r, const double y[], double dydr[], void *)
{
  // y[0] - R     
  // y[1] - Z

  double PsiR = GetPsiR (y[0], y[1]);
  double PsiZ = GetPsiZ (y[0], y[1]);

  dydr[0] = - qgp * y[0] * PsiZ;
  dydr[1] = + qgp * y[0] * PsiR;

  return GSL_SUCCESS;
}

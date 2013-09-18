// Flux.h

// Program to input CHEASE EQDSK equilibrium data and output Tomuhawc equilibrium data. 
//
// Radial grid is 
//
// Psi_j = 1. - s * tanh(s/sc)/tanh(1./sc)  for j = 0,J-1
//
// where s = j /(J-1),
//
// Poloidal grid is
//
// theta_k = M_PI * sin(0.5*M_PI*t/tc) /sin(0.5*M_PI*1./tc)  for k = 0,K-1
//
// where t = k /(K-1).
//
// Calculation control parameters in flux.in.
//
// Input data in EQDSK_COCOS_02.OUT  
// Intermediate data in folder Stage1/
// Final data in folder Stage2/ (see README)
//
// If called with no argument program prints usage message.
// If called with argument 0 program reads equilibrium data and exits.
// If called with non-zero argument program reads equilibrium data, constructs flux coordinate system, and exits.

#ifndef FLUX
#define FLUX

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// Maximum EQDSK resolution
#define N 1000

// Pointers to right-hand side function for adaptive integration
extern "C" int pRhs1  (double, const double [], double [], void *);
extern "C" int pRhs2  (double, const double [], double [], void *);
extern "C" int pRhs3  (double, const double [], double [], void *);

// ############
// Class header
// ############
class Flux
{
private:

  // Control parameters read from flux.in
  int    J;            // Number of points in r grid
  int    K;            // Number of points in theta grid
  double sc;           // Controls central concentration of r grid points
  double tc;           // Controls inboard concentration of th grid points
  double Gamma;        // Plasma ratio of specific heats
  double h0;           // Initial integration step-length
  double acc;          // Integration accuracy

  // Control parameter read from thawc.in
  int    ntor;         // Toroidal mode number

  // Equilibrium parameters
  double QA;           // CHEASE safety factor at wall
  double IASPCT;       // CHEASE aspect-ratio at wall
  double BETA;         // CHEASE beta
  double epsa;         // Tomuhawc inverse aspect-ratio at wall
  double epsb;         // Tomuhawc inverse aspect-ratio at plasma boundary
  double qb;           // Safety factor at plasma boundary
  double psib;         // Psi at plasma boundary
  double Raxis;        // R coordinate of magnetic axis
  double Zaxis;        // Z coordinate of magnetic axis
  double psic;         // Psi on magnetic axis
  int    ic;           // R grid index of magnetic axis
  int    jc;           // Z grid index of magnetic axis
  double qgp;

  // EQDSK parameters
  int         II;      // Number of R points
  double     *R;       // R array
  int         JJ;      // Number of Z points
  double     *Z;       // Z array
  gsl_matrix *Psi;     // Psi(R, Z) array
  int         KK;      // Number of boundary points
  double     *Rb;      // R on boundary
  double     *Zb;      // Z on boundary
  
  // Stage 1 profile parameters
  double *PSI;         // PSI array (from EQDSK)
  double *g;           // g(PSI)
  double *p;           // p(PSI)
  double *ggp;         // g dg/dPSI
  double *pp;          // dp/dPSI
  double *q;           // q(PSI)

  int     L;           // Number of points in Psi(R,0) array
  double *s;           // Array of s = sqrt[1 - Psi(R,0)] values
  double *Rs;          // Array of R(s) values

  // Stage 2 profile parameters
  double *P ;          // Psi array
  double *RP;          // R(Psi)
  double *rP;          // r(Psi)
  double *GP;          // g(Psi)
  double *QGP;         // q(psi)/g(psi) 
  double *QP;          // q(Psi)
  double *PP;          // P(Psi)
  double *GPP;         // dg/dPsi
  double *PPP;         // dP/dPsi
  double *S;           // sqrt(1 - Psi)
  double *QX;          // q(Psi) from EQDSK

  // Tomuhawc metric data
  double     *th;      // theta array
  gsl_matrix *Rst;     // [R(r,theta) - Raxis] /r
  gsl_matrix *Zst;     // [Z(r,theta) - Zaxis] /r
  gsl_matrix *Rt;      // (dR/dtheta) /r
  gsl_matrix *Zt;      // (dR/dtheta) /r
  gsl_matrix *Rr;      // (dR/dr) 
  gsl_matrix *Zr;      // (dZ/dr) 
  gsl_matrix *igradr2; // |grad r|^{-2}
  gsl_matrix *gradrt;  // r grad r . grad th / |grad r|^{-2}
  gsl_matrix *R2;      // (R^2 - R_axis^2) /r

public:

  // Constructor
  Flux ();
  // Set global parameters
  void SetParameters ();
  // Input EQDSK data and output intermediate data
  void Stage1 (int flag);
  // Input stage1 data and output Tomuhawc data
  void Stage2 ();
  // Evaluate right-hand sides of q/g equation
  int Rhs1 (double r, const double y[], double dydr[], void *);
  // Evaluate right-hand sides of r equation
  int Rhs2 (double r, const double y[], double dydr[], void *);
  // Evaluate right-hand sides of theta equation
  int Rhs3 (double r, const double y[], double dydr[], void *);

private:

  // Read EQDSK file.
  int EQDSK_Read (char *FILENAME,
		  int &INRBOX, int &INZBOX, 
		  double &RBOXLEN, double &ZBOXLEN, double &RBOXLFT, 
		  double &R0EXP, double &B0EXP,
		  double &RAXIS, double &ZAXIS, double &PSIAXIS, 
		  double xg[], double xP[], double xggp[], double xPp[], double xq[], 
		  double psi[][N],
		  int &NBPTS, double RBPTS[], double ZBPTS[],
		  double &IP, double &QC, double &QA, double &BETAP, double &LI);

  // Calculate q(P)/g(P) profile
  void CalcQGP ();
  // Calculate r(P) profile
  void CalcrP ();
  // Calculate straight angle
  void Calctst ();
  // Evaluate profile functions
  void Profile ();
  // Evaluate metric functions
  void Metric ();

  // Calculate boundary data
  void Boundary ();
  // Calculate toroidal ring function
  double ToroidalP (int m, int n, double z);
  // Calculate toroidal ring function
  double ToroidalQ (int m, int n, double z);

  // Interpolate Psi on uniform 2D grid
  double InterpolatePsi (double RR, double ZZ, int order);
  // 1D interpolation function with nonuniform grid
  double Interpolate (int I, double *X, double *Y, double x, int order);
  // 1D interpolation function with nonuniform grid
  double Interpolate1 (int I, double *X, gsl_matrix *Y, double x, int m, int order);
  // 1D extrapolation function with nonuniform grid
  double Extrapolate1 (int I, double *X, double *Y, double x);
  // 1D extrapolation function with nonuniform grid
  double Extrapolate2 (int I,  double *X, gsl_matrix *Y, double x, int k);
  // 1D extrapolation function with nonuniform grid
  double Extrapolate3 (int I,  double *X, gsl_matrix *Y, double x, int k);

  // Return cosine Fourier transform
  void FTCos (int j, double in[], gsl_matrix *out);
  // Return sine Fourier transform
  void FTSin (int j, double in[], gsl_matrix *out);

  // Open new file for writing
  FILE* OpenFile (char *filename);
  // Open existing file for reading
  FILE* OpenFiler (char *filename);

  // ...........
  // In Misc.cpp
  // ...........

  // Evaluate Psi (R, Z)
  double GetPsi (double r, double z);
  // Evaluate dPsi/dR (R, Z)
  double GetPsiR (double r, double z);
  // Evaluate dPsi/dZ (R, Z)
  double GetPsiZ (double r, double z);
  // Evaluate d^2Psi/dR^2 (R, Z)
  double GetPsiRR (double r, double z);
  // Evaluate d^2Psi/dZ^2 (R, Z)
  double GetPsiZZ (double r, double z);
};

#endif //FLUX

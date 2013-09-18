// Tomuhawc.h

// Tomuhawc (TOkamak-MUlti-Harmonic-tearing-stability-Analyis-With-Coupling).

// Class to calculate tearing/twisting stability matrix in tokamak plasma given
// equilibrium metric data.

// Reference: R. Fitzpatrick, R.J. Hastie, T.J. Martin, and C.M. Roach, 
// Nucl. Fusion 33, 1533 (1993).

// Normalization:         All lengths normalized to major radius of magnetic axis, R0.
//                        All magnetic field-strengths normalized to vacuum toroidal 
//                        field-strength on magnetic axis, B0. 
//                        All plasma pressures normalized to B0^2/\mu_0.

// Coordinates:           (r, theta, phi)
//                        r is flux-surface label (r=0 at magnetic axis, r=a at wall).
//                        theta is straight poloidal angle, phi is toroidal angle.
//                        Jacobian = r R^2, where R is major radius.

// Toroidal field:        Toroidal magnetic field-strength: 
//                          B_\phi = g(r)/R

// Pressure:              P(r) 

// Boundary condition:    Plasma surrounded by perfectly conducting wall at r=a.

// --------------------------------------------------------------------------------------------

// Integration control parameters:

//     flg = 0 ... truncation error is absolute
//     flg = 1 ... truncation error is relative
//     flg = 2 ... truncation error is mixed (default)

//    meth = 0 ... 4th order (classical) Runge-Kutta method
//    meth = 1 ... Embedded Runge-Kutta-Fehlberg (4,5) method
//    meth = 2 ... Embedded Runge-Kutta Cash-Karp (4,5) method
//    meth = 3 ... Embedded Runge-Kutta Prince-Dormand (8,9) method (default)
//    meth = 4 ... Implicit 4th order Runge-Kutta at Gaussian points
//    meth = 5 ... M=2 implicit Gear method

// Program uses GNU SCIENTIFIC LIBRARY (GSL) (http://www.gnu.org/software/gsl/) for 
// a) Variable dimension matrices and vectors.
// b) Interpolation.
// c) Linear algebra.
// d) Adaptive integration of systems of coupled ordinary differential equations.

// All interpolation uses cubic splines.

// --------------------------------------------------------------------------------------------

// Default values of class parameters specified in thawc.in

// --------------------------------------------------------------------------------------------

// Metric data (in folder Stage2). 

// q(r)     ... safety factor
// qp(r)    ... a dq/dr
// p1(r)    ... r^2
// p2(r)    ... (q/r/g) (dg/dr)
// p3(r)    ... (q^2/r/g^2) (dP/dr) 
// p4(r)    ... (q/g)(r d/dr)(g/q)
// p5(r)    ... r^2 g^2 /q^2
// p6(r)    ... Gamma P(r)
// p1p(r)   ... a d(p1)/dr
// p2p(r)   ... a d(p2)/dr
// p3p(r)   ... a d(p3)/dr

// M1    ... (R/R_0)^2
// M2    ... a^2 |Grad r|^(-2) R^(-2)
// M3    ... a^2 |Grad r|^(-2) 
// M4    ... a^2 |Grad r|^(-2) R^2
// M5    ... a^2 |Grad r|^(-2) R^4
// M6    ... |Grad r|^(-2) (r Grad r.Grad theta) 
// M7    ... |Grad r|^(-2) (r Grad r.Grad theta) R^2
// M8    ... a^(-2) |Grad r|^2
// M9    ... R^4
// M1P   ... a d(M1)/dr
// M3P   ... a d(M3)/dr
// M4P   ... a d(M4)/dr

// info.out
// Gamma, abar, a, rb

// grid.out
// nint, diag

// radial grid:          r/a[i] (for i = 0, nint-1)
// poloidal mode number: m    (for m = 0, diag)

// profile.out
// r/a[i], q[i], qp[i], p1[i], p2[i], p3[i], p4[i], p5[i], p6[i], p1p[i], p2p[i], p3p[i]\n 
// (for i = 0, nint-1)

// M1.out
// M1[i][m] (for m = 0, diag)\n (for i = 0, nint-1)

// M2.out
// M2[i][m] (for m = 0, diag)\n (for i = 0, nint-1)

// M3.out
// M3[i][m] (for m = 0, diag)\n (for i = 0, nint-1)

// M4.out
// M4[i][m] (for m = 0, diag)\n (for i = 0, nint-1)

// M5.out
// M5[i][m] (for m = 0, diag)\n (for i = 0, nint-1)

// M6.out
// M6[i][m] (for m = 0, diag)\n (for i = 0, nint-1)

// M7.out
// M7[i][m] (for m = 0, diag)\n (for i = 0, nint-1)

// M8.out
// M8[i][0]\n (for i = 0, nint-1)

// M9.out
// M9[i][0]\n (for i = 0, nint-1)

// M1P.out
// M1P[i][m] (for m = 0, diag)\n (for i = 0, nint-1)

// M3P.out
// M3P[i][0]\n (for i = 0, nint-1)

// M4P.out
// M4P[i][0]\n (for i = 0, nint-1)

// P.out
// Ptor0[i][j] (for j = 0, diag)\n (for i = 0, diag)

// Q.out
// Qtor0[i][j] (for j = 0, diag)\n (for i = 0, diag)

// PP.out
// Ptor1[i][j] (for j = 0, diag)\n (for i = 0, diag)

// QQ.out
// Qtor1[i][j] (for j = 0, diag)\n (for i = 0, diag)

// PQm.out
// Pm[j], Qm[j]\n (for j = 0, diag)

// --------------------------------------------------------

// Output files (in folder Stage3)

// logfile          ... logs details of calculation
// rational.out     ... num then m, r/a, s, E+F, H, K*E-G, K*F+G, K*H, f_A, f_S 
//                      (for each rational surface)
// ratsur.out       ... rs/a, qs
// Delta.out        ... num then diagonal elements of Ee matrix for each rational surface
// Fee.out          ... Fee matrix
// Foo.out          ... Foo matrix
// Foe.out          ... Foe matrix
// Feo.out          ... Feo matrix
// Ee.out           ... Ee  matrix
// Eo.out           ... Eo  matrix
// G.out            ... G   matrix
// Gp.out           ... Gp  matrix
// Bmat.out         ... coupling matrix B at r/a=rb
// Cmat.out         ... coupling matrix C at r/a=rb
// Dmat.out         ... coupling matrix D at r/a=rb
// Emat.out         ... coupling matrix E at r/a=rb
// Lmat.out         ... coupling matrix L at r/a=rb
// Mmat.out         ... coupling matrix M at r/a=rb
// Nmat.out         ... coupling matrix N at r/a=rb
// Pmat.out         ... coupling matrix P at r/a=rb
// info.out         ... dim, num, vac, off, b/a
// adape.out        ... r/a, log10(h), log10(Psi_max) for tearing parity solution
// adapo.out        ... r/a, log10(h), log10(Psi_max) for twisting parity solution
// GGJ.out          ... r/a, E+F, H, K*E-G, K*F+G, K*H, F_A, F_R 
// Yea.out          ... r/a, Y(r=a) for tearing parity solution
// Yoa.out          ... r/a, Y(r=a) for twisting parity solution
// psiFe.out        ... r/a, Psi(tearing  parity, fully reconnected)
// ZFe.out          ... r/a, Z  (tearing  parity, fully reconnected)
// psiFea.out       ... r/a, Psi(tearing  parity, fully reconnected) at r=a
// ZFea.out         ... r/a, Z  (tearing  parity, fully reconnected) at r=a
// psiUe.out        ... r/a, Psi(tearing  parity, unreconnected)
// ZUe.out          ... r/a, Z  (tearing  parity, unreconnected)
// psiFo.out        ... r/a, Psi(twisting parity, fully reconnected)
// ZFo.out          ... r/a, Z  (twisting parity, fully reconnected)
// psiFoa.out       ... r/a, Psi(twisting parity, fully reconnected) at r=a
// ZFoa.out         ... r/a, Z  (twisting parity, fully reconnected) at r=a
// psiUo.out        ... r/a, Psi(twisting parity, unreconnected)
// ZUo.out          ... r/a, Z  (twisting parity, unreconnected)
// xie.out          ... Error field response function xie (tearing  parity)
// xio.out          ... Error field response function xio (twisting parity)
// modes.out        ... poloidal mode numbers

// -----------------------------------------------------------

#ifndef THAWK
#define THAWK

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>

// Pointers to right-hand side function for adaptive integration
extern "C" int pRhs  (double, const double [], double [], void *);

// Global variables
extern int step, func;

// ############
// Class header
// ############
class Thawc
{
private:
  // ---------------------------------------
  // Metric iterpolation parameters and data
  // ---------------------------------------
  int                nint;    // Number of radial grid points 
  int                diag;    // Number of off-diagonal elements in coupling matrices
  double             abar;    // True minor radius
  double             a;       // Tomuhawc minor radius
  double             r0;      // Solutions launched from r/a = r0
  double             ra;      // Wall at r/a = ra = 1.
  double             qa;      // Safety-factor at wall
  double             rb;      // Plasma/vacuum boundary at r/a = b/a = rb
  double             Gamma;   // Plasma ratio of specific heats

  gsl_spline        *Sq;      // Interpolated safety-factor profile
  gsl_spline        *Sqp;     // Interpolated first derivative of safety-factor profile
  gsl_spline        *Sprof1;  // Interpolated first profile function
  gsl_spline        *Sprof2;  // Interpolated second profile function
  gsl_spline        *Sprof3;  // Interpolated third profile function
  gsl_spline        *Sprof4;  // Interpolated fourth profile function
  gsl_spline        *Sprof5;  // Interpolated fifth profile function
  gsl_spline        *Sprof6;  // Interpolated sixth profile function
  gsl_spline        *Sprof1p; // Interpolated derivative of first profile function
  gsl_spline        *Sprof2p; // Interpolated derivative of second profile function
  gsl_spline        *Sprof3p; // Interpolated derivative of third profile function
  gsl_spline       **SM1;     // Interpolated first metric element
  gsl_spline       **SM2;     // Interpolated second metric element
  gsl_spline       **SM3;     // Interpolated third metric element
  gsl_spline       **SM4;     // Interpolated fourth metric element
  gsl_spline       **SM5;     // Interpolated fifth metric element
  gsl_spline       **SM6;     // Interpolated sixth metric element
  gsl_spline       **SM7;     // Interpolated seventh metric element
  gsl_spline        *SM8;     // Interpolated eighth metric element
  gsl_spline        *SM9;     // Interpolated ninth metric element
  gsl_spline       **SM1P;    // Interpolated derivative of first metric element
  gsl_spline        *SM3P;    // Interpolated derivative of third metric element
  gsl_spline        *SM4P;    // Interpolated derivative of fourth metric element

  gsl_interp_accel  *Aq;      // Accelerator for interpolated q-profile
  gsl_interp_accel  *Aqp;     // Accelerator for interpolated qp-profile 
  gsl_interp_accel  *Aprof1;  // Accelerator for interpolated first profile function
  gsl_interp_accel  *Aprof2;  // Accelerator for interpolated second profile function
  gsl_interp_accel  *Aprof3;  // Accelerator for interpolated third profile function
  gsl_interp_accel  *Aprof4;  // Accelerator for interpolated fourth profile function
  gsl_interp_accel  *Aprof5;  // Accelerator for interpolated fifth profile function
  gsl_interp_accel  *Aprof6;  // Accelerator for interpolated sixth profile function
  gsl_interp_accel  *Aprof1p; // Accelerator for interpolated derivative of first prof. function
  gsl_interp_accel  *Aprof2p; // Accelerator for interpolated derivative of second prof. function
  gsl_interp_accel  *Aprof3p; // Accelerator for interpolated derivative of third prof. function
  gsl_interp_accel **AM1;     // Accelerator for interpolated first metric element
  gsl_interp_accel **AM2;     // Accelerator for interpolated second metric element
  gsl_interp_accel **AM3;     // Accelerator for interpolated third metric element
  gsl_interp_accel **AM4;     // Accelerator for interpolated fourth metric element
  gsl_interp_accel **AM5;     // Accelerator for interpolated fifth metric element
  gsl_interp_accel **AM6;     // Accelerator for interpolated sixth metric element
  gsl_interp_accel **AM7;     // Accelerator for interpolated seventh metric element 
  gsl_interp_accel  *AM8;     // Accelerator for interpolated eighth metric element
  gsl_interp_accel  *AM9;     // Accelerator for interpolated ninth metric element
  gsl_interp_accel **AM1P;    // Accelerator for interpolated derivative of first metric element
  gsl_interp_accel  *AM3P;    // Accelerator for interpolated derivative of second metric element
  gsl_interp_accel  *AM4P;    // Accelerator for interpolated derivative of third metric element

  gsl_matrix        *Ptor0;   // Vacuum solution data
  gsl_matrix        *Qtor0;   // Vacuum solution data
  gsl_matrix        *Ptor1;   // Vacuum solution data
  gsl_matrix        *Qtor1;   // Vacuum solution data
  double            *Pm;      // Vacuum solution data
  double            *Qm;      // Vacuum solution data

  // ---------------
  // Mode parameters
  // ---------------
  int     twist; // Twisting parity calculation enabled when twist = 1
  int     ntor;  // Common toroidal mode number
  int     side;  // Number of sideband harmonics
  int     off;   // Offset in sideband harmonics
  int     num;   // Number of rational surfaces in plasma
  int     vac;   // Total number of rational surfaces (including those in vacuum)

  int    *Mres;  // Resonant poloidal mode numbers
  double *Rres;  // Radii of resonant surfaces
  double *Qres;  // Safety factors at resonant surfaces
  double *Sres;  // Magnetic shears at resonant surfaces
  double *EFres; // GGJ parameters E+F at resonant surfaces
  double *HHres; // GGJ parameters H at resonant surfaces
  double *ETres; // GGJ parameters K*E-G at resonant surfaces
  double *FTres; // GGJ parameters K*F+G at resonant surfaces
  double *HTres; // GGJ parameters K*H at resonant surfaces 
  double *FRres; // GGJ parameters FR at resonant surfaces
  double *FAres; // GGJ parameters FA at resonant surfaces

  // -----------------
  // Coupling matrices
  // -----------------
  int         dim;  // Dimension of coupling matrices
  int         dim1; // Dimension of single solution vector
  int         dimN; // Dimension of composite solution vector

  int        *mpol; // Poloidal mode numbers of coupled harmonics
  gsl_matrix *bmat; // Coupling matrix
  gsl_matrix *cmat; // Coupling matrix
  gsl_matrix *dmat; // Coupling matrix
  gsl_matrix *emat; // Coupling matrix
  gsl_matrix *Lmat; // Coupling matrix
  gsl_matrix *Mmat; // Coupling matrix
  gsl_matrix *Nmat; // Coupling matrix
  gsl_matrix *Pmat; // Coupling matrix
  gsl_matrix *Vmat; // Vacuum matching matrix

  // ----------------
  // Fixup parameters
  // ----------------
  double r1;    // Last core fixup at r/a = r1 
  int    nfix;  // Number of fixups to perform

  // -------------------------------
  // Adaptive integration parameters
  // -------------------------------
  double h0;    // Initial step length 
  double acc;   // Integration accuracy parameter 
  int    flg;   // Controls calculation of truncation error (for mode calculation)
  int    meth;  // Controls choice of stepping method (for mode calculation)
  int    skip;  // Eigenfunctions calculated every skip grid points
  double del;   // Distance of closest approach to rational surfaces
  double nulc;  // Use zero pressure jump conditions when |nu_L| < nulc
  double eps;   // Step length for finite difference determination of derivatives
  int    nflag; // Indicates use of zero pressure jump condition

  // -----------------------
  // Root finding parameters
  // -----------------------
  double Eta;     // Min. magnitude of f at root f(x) = 0
  int    Maxiter; // Maximum number of iterations
  double qval;

  // ----------------------
  // Public class functions
  // ----------------------
public:
  
  // +++++++++++++++
  // in Tomahawk.cpp
  // +++++++++++++++
  // Constructor
  Thawc (); 
  // Destructor
  ~Thawc ();    
  // Status function    
  void Status (int);  

  // Function to solve stability problem 
  int Solv (int);   

  // Evaluate right-hand sides of mode equations
  int Rhs (double, const double [], double [], void *);

  // .................................
  // Functions to set class parameters
  // .................................  
  
  // +++++++++++++++
  // in Tomahawk.cpp
  // +++++++++++++++
  int Settwist (int);
  int Setntor  (int);
  int Setside  (int);
  int Setoff   (int);

  int Setr0    (double);
  int Setr1    (double);
  int Setnfix  (int);

  int SetEta   (double);

  int Seth0    (double);
  int Setacc   (double);
  int Setflg   (int);
  int Setmeth  (int);
  int Setskip  (int);
  int Setdel   (double);
  int Setnulc  (double);
  int Seteps   (double);

  // -----------------------
  // Private class functions
  // -----------------------
private:
 
  // Read in equilibrium data
  void Equilibrium (int, int);
  // Calculate resonant surface data
  void Resonant (int, int);
  // Calculate Glasser, Green, Johnson layer theory parameters
  void GGJCalc (double, double &, double &, double &, double &, double &, double &, double &);
  // Calculate vacuum matching matrix
  void Vacuum ();

  // Initialize single solution vector at magnetic axis
  void Launch1 (double, int, double []);
  // Initialize multiple solution vectors at magnetic axis
  void Launch (double, gsl_matrix *);
  // Integrate multiple solution vectors between specified radii
  // while performing nfix fixups
  void Segment_Fixup (double &, double, int, gsl_matrix *, gsl_matrix *, gsl_matrix *, 
		      int, int, int, int, char *, char *, char *);
  // Perform fixup on mutliple solution vectors
  void Fixup (gsl_matrix *, gsl_matrix *, gsl_matrix *, int, int, char *, char *);

  // Integrate single solution vector between specified radii
  void Segment1 (double &, double, double [], int, int, int, char *, char *, char *);
  // Integrate multiple solution vectors between specified radii
  void Segment (double &, double, gsl_matrix *, int, int, int, char *, char *, char *);
  // Calculate maximum Psi value
  double CalcPsiMax (double []);

  // Calculate elements of coupling matrices at given radius
  void Couple (double);

  // Calculate L1, P1, T1, s coefficients at given radius
  void L1P1T1 (double, double, int, double &, double &, double &, double &);
  // Calculate L1k, M1k, N1k, P1k coefficients at given radius
  void L1kP1k (double, double, int, double *, double *, double *, double *);
  // Evolve single solution vector across given rational surface
  void Jump1 (double, int, int, double [], double &, double &, int, int);
  // Evolve multiple solution vectors across given rational surface
  void Jump (double &, int, gsl_matrix *, gsl_matrix *, gsl_matrix *, int, int);

  // Apply boundary conditions at edge of plasma
  void Boundary (gsl_matrix *, gsl_matrix *, int);

  // Calculate Fee-matrix 
  void CalcFee (gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *);
  // Calculate Foe-matrix
  void CalcFoe (gsl_matrix *, gsl_matrix *, gsl_matrix *);
  // Calculate Foo-matrix 
  void CalcFoo (gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *);
  // Calculate Feo-matrix
  void CalcFeo (gsl_matrix *, gsl_matrix *, gsl_matrix *);

  // Calculate Ee and Gp matrices
  void CalcEe (gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *);
  // Calculate Eo and G matrices
  void CalcEo (gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *);

  // Calculate eigenfunctions of Fee-matrix
  void CalcEigFee (gsl_matrix *);
  // Calculate eigenfunctions of Fee-matrix at wall
  void CalcEigFeea (gsl_matrix *);
  // Calculate eigenfunctions of Foo-matrix
  void CalcEigFoo (gsl_matrix *);
  // Calculate eigenfunctions of Foo-matrix at wall
  void CalcEigFooa (gsl_matrix *);
  // Calculate eigenfunctions of Ee-matrix
  void CalcEigEe (gsl_matrix *);
  // Calculate eigenfunctions of Ee-matrix at wall
  void CalcEigEea (gsl_matrix *);
  // Calculate eigenfunctions of Eo-matrix
  void CalcEigEo (gsl_matrix *);
  // Calculate eigenfunctions of Eo-matrix at wall
  void CalcEigEoa (gsl_matrix *);

  // +++++++++++++++
  // in ZeroFind.cpp
  // +++++++++++++++
  // Target functions for zero finding
  double Feval (double);
  // Zero finding routine
  double RootFind ();
  // Ridder's method for finding root of F(x) = 0
  void Ridder (double, double, double, double, double &);

  // ++++++++++
  // in Log.cpp
  // ++++++++++
  // Log coupling matrices
  void LogMatrices ();
  // Log multiple solution vectors
  void LogVector (double, gsl_matrix *);
  // Log Y1-matrix
  void LogY1 (gsl_matrix *);
  // Log Fee-matrix
  void LogFee (gsl_matrix *);  
  // Log Foe-matrix
  void LogFoe (gsl_matrix *);
  // Log Ee-matrix
  void LogEe (gsl_matrix *, int);
  // Log Gp-matrix
  void LogGp (gsl_matrix *, int); 
  // Log Foo-matrix
  void LogFoo (gsl_matrix *);  
  // Log Feo-matrix
  void LogFeo (gsl_matrix *);
  // Log Eo-matrix
  void LogEo (gsl_matrix *, int);
  // Log G-matrix
  void LogG (gsl_matrix *, gsl_matrix *, int); 
  // Open new file for writing
  FILE* OpenFile (char *);
  // Open existing file for writing
  FILE* OpenFilea (char *);
  // Open existing file for reading
  FILE* OpenFiler (char *);
};

#endif //THAWK

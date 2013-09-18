// GGJ.h

// Class to solve Glasser-Green-Johnson layer equations.
// GGJ equations and parameters are defined in Glasser, Green, Johnson (Phys Fluids 18, 875, 1975).

// Uses finite-difference method of Glasser, Jardin, Tesauro (Phys Fluids 27, 1225, 1984).

// fQ   = fA**(1/3) / fA;
// fD   = fS**nuSL/3
// nuLS = (1. - 4.*(EF + HH))**(1/2);

//     Q_tomuhawc = fQ * Q_GGJ
// Delta_tomuhawc =  2 * Delta_GGJ /fD

// -----------------------------------------------------------------------------------------

#ifndef GGJGGJ
#define GGJGGJ

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

// ############
// Class header
// ############
class GGJ
{
private:
  // ----------------
  // Layer parameters
  // ----------------
  double EF;     // GGJ layer parameter E+F
  double H;      // GGJ layer parameter H
  double KFG;    // GGJ layer parameter K*F+G
  double KEG;    // GGJ layer parameter K*E-G
  double KH;     // GGJ layer parameter K*H
  double Qr;     // Real part of normalized growth-rate
  double Qi;     // Imaginary part of normalized growth-rate
  double fA;     // Normalization factor for Alfven frequency
  double fS;     // Normalization factor for Lundquist number

  // --------------------
  // Numerical parameters
  // --------------------
  int N;         // Number of grid-points 

  // -----------------------
  // Root finding parameters
  // -----------------------
  int    nint;    // Number of search intervals
  double Eta;     // Min. magnitude of f at root f(x) = 0
  int    Maxiter; // Maximum number of iterations

  // ----------------------
  // Public class functions
  // ----------------------
public:

  // Constructors
  GGJ ();
  GGJ (double, double, double, double, double, double, double);
  // Destructor
  ~GGJ ();   

  // Status function
  void Status ();  
  // Solve GGJ layer equations
  void Solv ();
  // Scan Deltae
  void ScanDeltae ();
  // Scan Deltao
  void ScanDeltao (); 
  // Scan Qr
  void ScanQr ();
  // Scan Qi
  void ScanQi ();
  // Find Q when Deltae=E
  void FindQe (double, double&, double&, double, double&, int&);
  // Find Q when Deltao=E
  void FindQo (double, double&, double&, double, double&, int&);
  // Calculate Deltae
  void GetDeltae (double, double, double&, double&);
  // Calculate Deltao
  void GetDeltao (double, double, double&, double&);
  // Find critical Deltae at which Re(Q)=0
  void FindCrite (double, double&, double&, double&);
 
  // ----------------------------------
  // Functions to set class parameters.
  // ----------------------------------
  int SetEF   (double);
  int SetH    (double);
  int SetKFG  (double);
  int SetKEG  (double);
  int SetKH   (double);
  int SetfA   (double);
  int SetfS   (double);
  int SetQr   (double);
  int SetQi   (double); 
  int SetN    (int);

  // -----------------------
  // Private class functions
  // -----------------------
private:
 
  // Find root of Deletae = x
  void FindRoote (double x, complex<double>& Q, complex<double>& Deltae, double& err, int& iter,
		  double Eps, int IterMax, double dQ);  
  // Find root of Deletao = x
  void FindRooto (double x, complex<double>& Q, complex<double>& Deltao, double& err, int& iter,
		  double Eps, int IterMax, double dQ);
  // Solve GGJ layer equations
  void Solv (complex<double>& Deltae, double& epse, double& dele,
	     complex<double>& Deltao, double& epso, double& delo, int flag);
  // Calculate Delta
  void CalcDelta (complex<double> Q, double L, int M, int sign, complex<double>& Delta, int flag);
  
  // Match asymptotic solution to finite-difference solution
  void Match (gsl_matrix_complex *DD, 
	      complex<double> *T3, complex<double> *T3p,
	      complex<double> *T4, complex<double> *T4p,
	      complex<double> *T5, complex<double> *T5p,
	      complex<double> *T6, complex<double> *T6p,
	      complex<double>& Delta, complex<double>& Eplus, complex<double>& Eminu);
  // Calculate asymptotic solutions
  void CalcAsymptotic (complex<double> Q, double x, int M,
		       double p3, double p4,  complex<double>  s5,  complex<double>  s6, 
		       complex<double> *psi3, complex<double> *xi3, complex<double> *v3,
		       complex<double> *psi4, complex<double> *xi4, complex<double> *v4,
		       complex<double> *psi5, complex<double> *xi5, complex<double> *v5,
		       complex<double> *psi6, complex<double> *xi6, complex<double> *v6,
		       complex<double> *T3,   complex<double> *T3p,
		       complex<double> *T4,   complex<double> *T4p,
		       complex<double> *T5,   complex<double> *T5p,
		       complex<double> *T6,   complex<double> *T6p);
  // Calculate asymptotic series
  void CalcSeries (complex<double> Q, int M,
		   double& p3, double& p4, complex<double>& s5,  complex<double>& s6,
		   complex<double> *psi3,  complex<double> *xi3, complex<double> *v3,
		   complex<double> *psi4,  complex<double> *xi4, complex<double> *v4,
		   complex<double> *psi5,  complex<double> *xi5, complex<double> *v5,
		   complex<double> *psi6,  complex<double> *xi6, complex<double> *v6);
  // Calculate A, B, and C matrices
  void CalcABC (complex<double> Q, double x, double h, 
		gsl_matrix_complex *A, gsl_matrix_complex *B, gsl_matrix_complex *C);
  // Calculate U, U', U'', V, etc. matrices
  void CalcUV (complex<double> Q, double x, 
	       gsl_matrix_complex *U,   gsl_matrix_complex *Up,  gsl_matrix_complex *Upp, 
	       gsl_matrix_complex *V,   gsl_matrix_complex *I,   gsl_matrix_complex *UV,  
	       gsl_matrix_complex *VU,  gsl_matrix_complex *VUp, gsl_matrix_complex *V2, 
	       gsl_matrix_complex *V2U, gsl_matrix_complex *U2,  gsl_matrix_complex *V3);

  // Add 3x3 matrices: C = A + B
  void AddMatrix (gsl_matrix_complex *A, gsl_matrix_complex *B, gsl_matrix_complex *C);
  // Add 3x1 vectors: c = a + a
  void AddVector (complex<double> *a, complex<double> *b, complex<double> *c);
  // Multiply 3x3 matrices: C = A*B
  void MultMatrix (gsl_matrix_complex *A, gsl_matrix_complex *B, gsl_matrix_complex *C);
  // Multiply 3x3 matrix and 3x1 vector: c = A*b
  void MultMatVec (gsl_matrix_complex *A, complex<double> *b, complex<double> *c);
  // Invert 3x3 matrix/vector equation: c = g*A^(-1)*b
  void InvMatVec (gsl_matrix_complex *A, complex<double> *b, complex<double> *c, complex<double> g);
  // Invert 3x3 matrix equation: C = g*A^(-1)*B
  void InvMatrix (gsl_matrix_complex *A, gsl_matrix_complex *B, gsl_matrix_complex *C, complex<double> g);
  // Calculate cross product of 3x1 vectors: c = axb
  void CrossProduct (complex<double> *a, complex<double> *b, complex<double> *c);
  // Calculate dot product of 3x1 vectors: c = a.b
  complex<double> DotProduct (complex<double> *a, complex<double> *b);
  // Load 3x1 vector
  void LoadVector (complex<double> *a, complex<double> a0, complex<double> a1, complex<double> a2);
  // Unload 3x1 vector
  void UnloadVec (complex<double> *a, complex<double>& a0, complex<double>& a1, complex<double>& a2);

  // Convert complex<double> into gsl_complex
  gsl_complex GslComplex (complex<double> z);
  // Convert double into gsl_complex
  gsl_complex GslComplex (double x);
  // Convert gsl_complex into complex<double>
  complex<double> ComplexDouble (gsl_complex z);
  // Return modulus of complex<double>
  double hypot (complex<double> z);

  // Target functions for zero finding
  double Feval (double);
  // Zero finding routine
  double RootFind (double, double);
  // Ridder's method for finding root of F(x) = 0
  void Ridder (double, double, double, double, double&);

  // Open file for writing 
  FILE* OpenFile (char *filename);
};

#endif //GGJGGJ

// GGJ.cpp

#include "GGJ.h"

// ############
// Constructors
// ############
GGJ::GGJ ()
{
  // Set default values of class parameters

  // betan = 0.1
  /*
  EF  = -4.228978e-03;
  H   =  1.325602e-03;
  KEG = -4.680152e+02;
  KFG =  4.546381e+02;
  KH  =  4.193132e+00; 
  fA  =  5.730443e-01;
  fS  =  2.289473e+00;
  */

  // betan = 0.25
  EF  = -1.106950e-02;
  H   =  3.766154e-03;
  KEG = -1.872715e+02;
  KFG =  1.817157e+02;
  KH  =  1.890262e+00;
  fA  =  5.749192e-01;
  fS  =  2.289415e+00;

  // betan = 1.5
  /*
  EF  = -9.481735e-02;
  H   =  4.475351e-02;
  KEG = -3.001516e+01;
  KFG =  2.888337e+01;
  KH  =  5.342013e-01;
  fA  =  5.808478e-01;
  fS  =  2.320894e+00;
  */

  Qr      = 0.01;
  Qi      = 0.;
  
  N       = 800; 

  nint    = 10;
  Eta     = 1.e-12;
  Maxiter = 60;
}

GGJ::GGJ (double _EF, double _H, double _KEG, double _KFG, double _KH, double _fA, double _fS)
{
  N       = 800;

  nint    = 10;
  Eta     = 1.e-12;
  Maxiter = 60;

  EF      = _EF;
  H       = _H;
  KEG     = _KEG;
  KFG     = _KFG;
  KH      = _KH;
  fA      = _fA;
  fS      = _fS;
}

// ##########
// Destructor
// ##########
GGJ::~GGJ ()
{}

// ###############
// Status function
// ###############
void GGJ::Status ()
{ 
  printf ("\nEF = %11.4e  H = %11.4e  KFG = %11.4e  KEG = %11.4e  KH = %11.4e  fA = %11.4e  fS = %11.4e\n\n",
	  EF, H, KFG, KEG, KH, fA, fS);
  printf ("Qr = %11.4e  Qi = %11.4e\n\n", 
	  Qr, Qi);
  printf ("N = %4d\n\n",
	  N);
}

// #################################
// Functions to set class parameters
// #################################
int GGJ::SetEF (double _EF)
{
  if (0)
    {
      printf ("GGJ::SetEF: Error - invalid data value\n");
      return 1;
    }
  EF = _EF;
  return 0;
}
int GGJ::SetH (double _H)
{
  if (0)
    {
      printf ("GGJ::SetH: Error - invalid data value\n");
      return 1;
    }
  H = _H;
  return 0;
}
int GGJ::SetKFG (double _KFG)
{
  if (0)
    {
      printf ("GGJ::SetKFG: Error - invalid data value\n");
      return 1;
    }
  KFG = _KFG;
  return 0;
}
int GGJ::SetKEG (double _KEG)
{
  if (0)
    {
      printf ("GGJ::SetKEG: Error - invalid data value\n");
      return 1;
    }
  KEG = _KEG;
  return 0;
}
int GGJ::SetKH (double _KH)
{
  if (0)
    {
      printf ("GGJ::SetKH: Error - invalid data value\n");
      return 1;
    }
  KH = _KH;
  return 0;
}
int GGJ::SetfA (double _fA)
{
  if (_fA < 0.)
    {
      printf ("GGJ::SetfA: Error - invalid data value\n");
      return 1;
    }
  fA = _fA;
  return 0;
}
int GGJ::SetfS (double _fS)
{
  if (_fS < 0.)
    {
      printf ("GGJ::SetfS: Error - invalid data value\n");
      return 1;
    }
  fS = _fS;
  return 0;
}
int GGJ::SetQr (double _Qr)
{
  if (0)
    {
      printf ("GGJ::SetQr: Error - invalid data value\n");
      return 1;
    }
  Qr = _Qr;
  return 0;
}
int GGJ::SetQi (double _Qi)
{
  if (0)
    {
      printf ("GGJ::SetQi: Error - invalid data value\n");
      return 1;
    }
  Qi = _Qi;
  return 0;
}
int GGJ::SetN (int _N)
{
  if (_N < 4)
    {
      printf ("GGJ::SetN: Error - invalid data value\n");
      return 1;
    }
  N = _N;
  return 0;
}

// #######################################################
// Function to find critical TONUHAWC Deltae at which Qr=0
// #######################################################
void GGJ::FindCrite (double x2, double& _Qi, double& Deltar, double& Deltai)
{
  double x1 = -1.e-6;

  _Qi = RootFind (x1, x2);
  GetDeltae (0., _Qi, Deltar, Deltai);
}

// #############################################################
// Function to find TOMUHAWC growth-rate when TOMUHAWC Deltae=E.
// #############################################################
void GGJ::FindQe (double E, double& _Qr, double& _Qi, double dQ, double& err, int& iter)
{ 
  double nuSL = sqrt (1. - 4.*(EF + H));
  double fQ   = pow (fS, 1./3.) /fA;
  double fD   = pow (fS, nuSL/3.);
  double Eps = 1.e-10, x = fD*E/2.;
  int ND, IterMax = Maxiter;
  complex<double> Q = complex<double> (_Qr /fQ, _Qi /fQ), Deltae;

  FindRoote (x, Q, Deltae, err, iter, Eps, IterMax, dQ);

  _Qr = fQ*real(Q);
  _Qi = fQ*imag(Q);
}

// #############################################################
// Function to find TOMUHAWC growth-rate when TOMUHAWC Deltao=E.
// #############################################################
void GGJ::FindQo (double E, double& _Qr, double& _Qi, double dQ, double& err, int& iter)
{
  double nuSL = sqrt (1. - 4.*(EF + H));
  double fQ   = pow (fS, 1./3.) /fA;
  double fD   = pow (fS, nuSL/3.);
  double Eps = 1.e-8, x = fD*E/2.;
  int ND, IterMax = Maxiter;
  complex<double> Q = complex<double> (_Qr /fQ, _Qi /fQ), Deltao;
       
  FindRooto (x, Q, Deltao, err, iter, Eps, IterMax, dQ);

  _Qr = fQ*real(Q);
  _Qi = fQ*imag(Q);
}

// ###############################################################
// Function to calculate TOMUHAWC Deltae as function of TOMUHAWC Q
// ###############################################################
void GGJ::GetDeltae (double _Qr, double _Qi, double& Deltar, double& Deltai)
{ 
  double Qrsave = Qr;
  double Qisave = Qi;

  double nuSL = sqrt (1. - 4.*(EF + H));
  double fQ   = pow (fS, 1./3.) /fA;
  double fD   = pow (fS, nuSL/3.);

  complex<double> Deltae, Deltao;
  double epse, dele, epso, delo;

  Qr = _Qr /fQ;
  Qi = _Qi /fQ;
  Solv (Deltae, epse, dele, Deltao, epso, delo, 0);
  
  Deltar = 2.*real(Deltae) /fD;
  Deltai = 2.*imag(Deltae) /fD;

  Qr = Qrsave;
  Qi = Qisave;
}

// #################################################
// Function to scan TOMUHAWC Qr at fixed TOMUHAWC Qi
// and output TOMUHAWC Delta
// #################################################
void GGJ::ScanQr ()
{
  double Qrsave = Qr;
  double Qisave = Qi;

  double nuSL = sqrt (1. - 4.*(EF + H));
  double fQ   = pow (fS, 1./3.) /fA;
  double fD   = pow (fS, nuSL/3.);

  double _Qi, Qstart, Qend, epse, epso, dele, delo;
  int NQ;
  complex<double> Deltae, Deltao;
  
  printf ("\nQi Qstart Qend NQ ?? ");
  scanf ("%lf %lf %lf %d", &_Qi, &Qstart, &Qend, &NQ);
  Qi = _Qi /fQ;

  FILE* file = OpenFile ("Stage4/scanqr.out");
  for (int i = 0; i < NQ; i++)
    {
      Qr = Qstart + (Qend - Qstart)* double (i) / double (NQ - 1);
      Qr /= fQ;
      
      Solv (Deltae, epse, dele, Deltao, epso, delo, 1);
 
      printf ("Q = (%10.3e, %10.3e) Deltae = (%10.3e, %10.3e) eps = %10.3e del = %10.3e Deltao = (%10.3e, %10.3e) eps = %10.3e del=%10.3e\n",
	      fQ*Qr, fQ*Qi, 2.*real(Deltae)/fD, 2.*imag(Deltae)/fD, epse, dele,  
	      2.*real(Deltao)/fD, 2.*imag(Deltao)/fD, epso, delo);
      fprintf (file, "%e %e %e %e %e %e %e %e %e %e\n",
	       fQ*Qr, fQ*Qi, 2.*real(Deltae)/fD, 2.*imag(Deltae)/fD, epse, dele,
	       2.*real(Deltao)/fD, 2.*imag(Deltao)/fD, epso, delo);
      fflush (file);
    }
  fclose (file);

  Qr = Qrsave;
  Qi = Qisave;
}

// #################################################
// Function to scan TOMUHAWC Qi at fixed TOMUHAWC Qr
// and output TOMUHAWC Delta
// #################################################
void GGJ::ScanQi ()
{
  double Qrsave = Qr;
  double Qisave = Qi;

  double nuSL = sqrt (1. - 4.*(EF + H));
  double fQ   = pow (fS, 1./3.) /fA;
  double fD   = pow (fS, nuSL/3.);

  double _Qr, Qstart, Qend, epse, epso, dele, delo;
  int NQ;
  complex<double> Deltae, Deltao;
  
  printf ("\nQr Qstart Qend NQ ?? ");
  scanf ("%lf %lf %lf %d", &_Qr, &Qstart, &Qend, &NQ);
  Qr = _Qr /fQ;

  FILE* file = OpenFile ("Stage4/scanqi.out");
  for (int i = 0; i < NQ; i++)
    {
      Qi  = Qstart + (Qend - Qstart)* double (i) / double (NQ - 1);
      Qi /= fQ;
      
      Solv (Deltae, epse, dele, Deltao, epso, delo, 1);
 
      printf ("Q = (%10.3e, %10.3e) Deltae = (%10.3e, %10.3e) eps = %10.3e del = %10.3e Deltao = (%10.3e, %10.3e) eps = %10.3e del=%10.3e\n",
	      fQ*Qr, fQ*Qi, 2.*real(Deltae)/fD, 2.*imag(Deltae)/fD, epse, dele, 
	      2.*real(Deltao) /fD, 2.*imag(Deltao) /fD, epso, delo);
      fprintf (file, "%e %e %e %e %e %e %e %e %e %e\n",
	       fQ*Qr, fQ*Qi, 2.*real(Deltae)/fD, 2.*imag(Deltae)/fD, epse, dele,
	       2.*real(Deltao)/fD, 2.*imag(Deltao/fD), epso, delo);
      fflush (file);
    }
  fclose (file);

  Qr = Qrsave;
  Qi = Qisave;
}

// #######################################################
// Function to scan TOMUHAWC Deltae and return TOMUHAWC Q.
// #######################################################
void GGJ::ScanDeltae ()
{
  double Dstart, Dend, Eps = 1.e-10, dQ = 1.e-8, err, x;
  int ND, iter, IterMax = Maxiter;
  complex<double> Q = complex<double> (Qr, Qi), Deltae;
  
  double nuSL = sqrt (1. - 4.*(EF + H));
  double fQ   = pow (fS, 1./3.) /fA;
  double fD   = pow (fS, nuSL/3.);

  printf ("\nDstart Dend ND ?? ");
  scanf ("%lf %lf %d", &Dstart, &Dend, &ND);

  FILE* file = OpenFile ("Stage4/scande.out");
  for (int i = 0; i < ND; i++)
    {
      x  = Dstart + (Dend - Dstart)* double (i) / double (ND - 1);
      x *= fD/2.; 
       
      FindRoote (x, Q, Deltae, err, iter, Eps, IterMax, dQ);
      printf  ("Deltae = (%11.4e, %11.4e)  Q = (%14.7e, %14.7e)  err = %11.4e  iter = %3d\n",
	       2.*real(Deltae)/fD, 2.*imag(Deltae)/fD, fQ*real(Q), fQ*imag(Q), err, iter);
      fprintf (file, "%e %e %e %e %e %d\n",
	       fQ*real(Q), fQ*imag(Q), 2.*real(Deltae)/fD, 2.*imag(Deltae)/fD, err, iter);
      fflush  (file);
    }
  fclose (file);
}

// #######################################################
// Function to scan TOMUHAWC Deltae and return TOMUHAWC Q.
// #######################################################
void GGJ::ScanDeltao ()
{
  double Dstart, Dend, Eps = 1.e-6, dQ = 1.e-8, err, x;
  int ND, iter, IterMax = Maxiter;
  complex<double> Q = complex<double> (Qr, Qi), Deltao;
 
  double nuSL = sqrt (1. - 4.*(EF + H));
  double fQ   = pow (fS, 1./3.) /fA;
  double fD   = pow (fS, nuSL/3.);
 
  printf ("\nDstart Dend ND ?? ");
  scanf ("%lf %lf %d", &Dstart, &Dend, &ND);

  FILE* file = OpenFile ("Stage4/scando.out");
  for (int i = 0; i < ND; i++)
    {
      x = Dstart + (Dend - Dstart)* double (i) / double (ND - 1);
      x *= fD/2.; 
       
      FindRooto (x, Q, Deltao, err, iter, Eps, IterMax, dQ);
      printf  ("Deltao = (%11.4e, %11.4e)  Q = (%14.7e, %14.7e)  err = %11.4e  iter = %3d\n",
	       2.*real(Deltao)/fD, 2.*imag(Deltao)/fD, fQ*real(Q), fQ*imag(Q), err, iter);
      fprintf (file, "%e %e %e %e %e %d\n",
	       fQ*real(Q), fQ*imag(Q), 2.*real(Deltao)/fD, 2.*imag(Deltao)/fD, err, iter);
      fflush (file);
    }
  fclose (file);
}

// ################################################################
// Function to find root of Deltae = x.
// Search via Newton iteration until |Deltae-x| < Eps 
// or number of iterations exceeds IterMax.
// On entry Q contains first guess. On exit Q contains root 
// and iter contain number of iterations.
// Jacobian calculated via finite differencing with grid spacing dQ.
// ################################################################
void GGJ::FindRoote (double x, complex<double>& Q, complex<double>& Deltae, double& err, int& iter,
		     double Eps, int IterMax, double dQ)
{
  double Qrsave = Qr;
  double Qisave = Qi;

  complex<double> Deltao;
  double epse, dele, epso, delo, errold;

  iter = 0;
  Qr = real (Q);
  Qi = imag (Q);
  Solv (Deltae, epse, dele, Deltao, epso, delo, 0);
  err = hypot (Deltae-x);
  if (hypot (Deltae) > 1.)
    err /= hypot (Deltae);
  //printf ("Q = (%11.4e, %11.4e) Deltae = (%11.4e, %11.4e) err = (%11.4e, %11.4e) eps = %11.4e del = %11.4e\n", 
  //  Qr, Qi, real(Deltae), imag(Deltae), real(Deltae)-x, imag(Deltae), epse, dele);

  if (err  < Eps) 
    {
      Qr = Qrsave;
      Qi = Qisave;
      return;
    }

  do
    {
      complex<double> Deltae1, Deltae2, dDdQ;
      double J00, J01, J10, J11, det, dQr, dQi;
      
      iter += 1;
      Qr   += dQ;
      Solv (Deltae1, epse, dele, Deltao, epso, delo, 0);
      Qr   -=  2.*dQ;
      Solv (Deltae2, epse, dele, Deltao, epso, delo, 0);
      dDdQ  = (Deltae1 - Deltae2) /2./dQ;
      J00   = real (dDdQ);
      J10   = imag (dDdQ);
      Qr   += dQ;

      Qi   += dQ;
      Solv (Deltae1, epse, dele, Deltao, epso, delo, 0);
      Qi   -= 2.*dQ;
      Solv (Deltae2, epse, dele, Deltao, epso, delo, 0);
      dDdQ  = (Deltae1 - Deltae2) /2./dQ;
      J01   = real (dDdQ);
      J11   = imag (dDdQ);
      Qi   += dQ;

      det = J00*J11 - J01*J10;

      dQr = (- J11 * (real(Deltae)-x) + J01 * imag(Deltae)) /det;
      dQi = (  J10 * (real(Deltae)-x) - J00 * imag(Deltae)) /det;

      Qr += dQr;
      Qi += dQi;

      errold = err;
      Solv (Deltae, epse, dele, Deltao, epso, delo, 0);
      err = hypot (Deltae-x);
      if (hypot (Deltae) > 1.)
	err /= hypot (Deltae);

      //  printf ("Q = (%11.4e, %11.4e) Deltae = (%11.4e, %11.4e) err = (%11.4e, %11.4e) eps = %11.4e del = %11.4e\n", 
      //      Qr, Qi, real(Deltae), imag(Deltae), real(Deltae)-x, imag(Deltae), epse, dele);
    }
  while (iter < IterMax && err > Eps);
  
  Q  = complex<double> (Qr, Qi);
  //err = (epse > err) ? epse : err;
  //err = (dele > err) ? dele : err;

  Qr = Qrsave;
  Qi = Qisave;
}

// ################################################################
// Function to find root of Deltao = x.
// Search via Newton iteration until |Deltao-x| < Eps 
// or number of iterations exceeds IterMax.
// On entry Q contains first guess. On exit Q contains root 
// and iter contain number of iterations.
// Jacobian calculated via finite differencing with grid spacing dQ.
// ################################################################
void GGJ::FindRooto (double x, complex<double>& Q, complex<double>& Deltao, double& err, int& iter,
		     double Eps, int IterMax, double dQ)
{
  double Qrsave = Qr;
  double Qisave = Qi;

  complex<double> Deltae;
  double epse, dele, epso, delo, errold;

  iter = 0;
  Qr = real (Q);
  Qi = imag (Q);
  Solv (Deltae, epse, dele, Deltao, epso, delo, 0);
  err = hypot (Deltao-x);
  if (hypot (Deltao) > 1.)
    err /= hypot (Deltao);
  // printf ("Q = (%11.4e, %11.4e) Deltao = (%11.4e, %11.4e) err = (%11.4e, %11.4e) eps = %11.4e del = %11.4e\n", 
  //	  Qr, Qi, real(Deltao), imag(Deltao), real(Deltao)-x, imag(Deltao), epso, delo);

  if (err  < Eps) 
    {
      Qr = Qrsave;
      Qi = Qisave;
      return;
    }

  do
    {
      complex<double> Deltao1, Deltao2, dDdQ;
      double J00, J01, J10, J11, det, dQr, dQi;
      
      iter += 1;
      Qr   += dQ;
      Solv (Deltae, epse, dele, Deltao1, epso, delo, 0);
      Qr   -=  2.*dQ;
      Solv (Deltae, epse, dele, Deltao2, epso, delo, 0);
      dDdQ  = (Deltao1 - Deltao2) /2./dQ;
      J00   = real (dDdQ);
      J10   = imag (dDdQ);
      Qr   += dQ;

      Qi   += dQ;
      Solv (Deltae, epse, dele, Deltao1, epso, delo, 0);
      Qi   -= 2.*dQ;
      Solv (Deltae, epse, dele, Deltao2, epso, delo, 0);
      dDdQ  = (Deltao1 - Deltao2) /2./dQ;
      J01   = real (dDdQ);
      J11   = imag (dDdQ);
      Qi   += dQ;

      det = J00*J11 - J01*J10;

      dQr = (- J11 * (real(Deltao)-x) + J01 * imag(Deltao)) /det;
      dQi = (  J10 * (real(Deltao)-x) - J00 * imag(Deltao)) /det;

      Qr += dQr;
      Qi += dQi;

      errold = err;
      Solv (Deltae, epse, dele, Deltao, epso, delo, 0);
      err = hypot (Deltao-x);
      if (hypot (Deltao) > 1.)
	err /= hypot (Deltao);

      //printf ("Q = (%11.4e, %11.4e) Deltao = (%11.4e, %11.4e) err = (%11.4e, %11.4e) eps = %11.4e del = %11.4e\n", 
      //	      Qr, Qi, real(Deltao), imag(Deltao), real(Deltao)-x, imag(Deltao), epso, delo);
    }
  while (iter < IterMax && err > Eps);
  
  Q  = complex<double> (Qr, Qi);
  err = (epse > err) ? epse : err;
  err = (dele > err) ? dele : err;

  Qr = Qrsave;
  Qi = Qisave;
}

// #####################################
// Function to solve GGJ layer equations
// #####################################
void GGJ::Solv ()
{
  complex<double> Deltae, Deltao;
  double epse, dele, epso, delo;

  Solv (Deltae, epse, dele, Deltao, epso, delo, 1);

  printf ("Q = (%11.4e, %11.4e)  Deltae = (%11.4e, %11.4e)  eps = %11.4e  del = %11.4e\n", 
	  Qr, Qi, real(Deltae), imag(Deltae), epse, dele);
  printf ("Q = (%11.4e, %11.4e)  Deltao = (%11.4e, %11.4e)  eps = %11.4e  del = %11.4e\n", 
	  Qr, Qi, real(Deltao), imag(Deltao), epso, delo);
}

// #####################################
// Function to solve GGJ layer equations
// #####################################
void GGJ::Solv (complex<double>& Deltae, double& epse, double& dele,
		complex<double>& Deltao, double& epso, double& delo, int flag)
{
  // Set Q
  complex<double> Q = complex<double> (Qr, Qi);

  // Set L
  double L0;
  if (hypot(Q) < 1.)
    L0 = 6.*pow(hypot(Q), 0.25);
  else
    {
      double Z = 1.;
      if (fabs(KFG) > Z)
	Z = fabs(KFG);
      if (fabs(KH*H) > Z)
	Z = fabs(KH*H);

      L0 = 2.*hypot(Q)*sqrt(Z);
    }

  // Set M
  int M0 = 9;

  complex<double> Deltae1, Deltao1, Deltae2, Deltao2;
  CalcDelta (Q, L0+0.05*L0, M0,  1, Deltae1, 0);
  CalcDelta (Q, L0+0.05*L0, M0, -1, Deltao1, 0);
  CalcDelta (Q, L0-0.05*L0, M0,  1, Deltae2, 0);
  CalcDelta (Q, L0-0.05*L0, M0, -1, Deltao2, 0);
  Deltae = (Deltae1 + Deltae2)/2.;
  Deltao = (Deltao1 + Deltao2)/2.;
  epse = hypot ((Deltae1 - Deltae2)/0.1/L0);
  epso = hypot ((Deltao1 - Deltao2)/0.1/L0);
  CalcDelta (Q, L0+0.05*L0, M0+1,  1, Deltae2, flag);
  CalcDelta (Q, L0+0.05*L0, M0+1, -1, Deltao2, 0);
  dele = hypot ((Deltae1 - Deltae2)/1.);
  delo = hypot ((Deltao1 - Deltao2)/1.);
}

// ###########################
// Function to calculate Delta
// ###########################
void GGJ::CalcDelta (complex<double> Q, double L, int M, int sign, complex<double>& Delta, int flag)
{
  // Allocate memory
  gsl_matrix_complex *EE  = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *A   = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *B   = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *C   = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *S   = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *T   = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *U   = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *Up  = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *Upp = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *V   = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *I   = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *UV  = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *VU  = gsl_matrix_complex_alloc (3, 3); 
  gsl_matrix_complex *VUp = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *V2  = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *V2U = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *U2  = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *V3  = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *DD  = gsl_matrix_complex_alloc (3, 3);

  double *X = new double[N];

  complex<double> *psi3 = new complex<double>[M];
  complex<double> *xi3  = new complex<double>[M];
  complex<double> *v3   = new complex<double>[M];
  complex<double> *psi4 = new complex<double>[M];
  complex<double> *xi4  = new complex<double>[M];
  complex<double> *v4   = new complex<double>[M];
  complex<double> *psi5 = new complex<double>[M];
  complex<double> *xi5  = new complex<double>[M];
  complex<double> *v5   = new complex<double>[M];
  complex<double> *psi6 = new complex<double>[M];
  complex<double> *xi6  = new complex<double>[M];
  complex<double> *v6   = new complex<double>[M];

  complex<double> *T3   = new complex<double>[3];
  complex<double> *T3p  = new complex<double>[3];
  complex<double> *T4   = new complex<double>[3];
  complex<double> *T4p  = new complex<double>[3];
  complex<double> *T5   = new complex<double>[3];
  complex<double> *T5p  = new complex<double>[3];
  complex<double> *T6   = new complex<double>[3];
  complex<double> *T6p  = new complex<double>[3];

  complex<double> *a    = new complex<double>[3];
  complex<double> *b    = new complex<double>[3];
  complex<double> *c    = new complex<double>[3];

  complex<double> *Psi  = new complex<double>[N];
  complex<double> *Xi   = new complex<double>[N];
  complex<double> *Ups  = new complex<double>[N];

  // Set h
  double h = L / double(N-1);

  // Set up grid
  for (int i = 0; i < N; i++)
    X[i] = double(i) * h;

  // Initialize EE-matrix
  CalcABC (Q, X[0], h, A, B, C);
  if (sign < 0)
    {
      complex<double> sum;
      sum = 
	  ComplexDouble (gsl_matrix_complex_get (A, 0, 0))
	- ComplexDouble (gsl_matrix_complex_get (C, 0, 0));
      gsl_matrix_complex_set (S, 0, 0, GslComplex (sum));
      sum = 
	  ComplexDouble (gsl_matrix_complex_get (A, 1, 0))
	- ComplexDouble (gsl_matrix_complex_get (C, 1, 0));
      gsl_matrix_complex_set (S, 1, 0, GslComplex (sum));
      sum = 
	  ComplexDouble (gsl_matrix_complex_get (A, 2, 0))
	- ComplexDouble (gsl_matrix_complex_get (C, 2, 0));
      gsl_matrix_complex_set (S, 2, 0, GslComplex (sum));
      sum = 
	  ComplexDouble (gsl_matrix_complex_get (A, 0, 1))
	+ ComplexDouble (gsl_matrix_complex_get (C, 0, 1));
      gsl_matrix_complex_set (S, 0, 1, GslComplex (sum));
      sum = 
	  ComplexDouble (gsl_matrix_complex_get (A, 1, 1))
	+ ComplexDouble (gsl_matrix_complex_get (C, 1, 1));
      gsl_matrix_complex_set (S, 1, 1, GslComplex (sum));
      sum = 
	  ComplexDouble (gsl_matrix_complex_get (A, 2, 1))
	+ ComplexDouble (gsl_matrix_complex_get (C, 2, 1));
      gsl_matrix_complex_set (S, 2, 1, GslComplex (sum));
      sum = 
	  ComplexDouble (gsl_matrix_complex_get (A, 0, 2))
	+ ComplexDouble (gsl_matrix_complex_get (C, 0, 2));
      gsl_matrix_complex_set (S, 0, 2, GslComplex (sum));
      sum = 
	  ComplexDouble (gsl_matrix_complex_get (A, 1, 2))
	+ ComplexDouble (gsl_matrix_complex_get (C, 1, 2));
      gsl_matrix_complex_set (S, 1, 2, GslComplex (sum));
      sum = 
	  ComplexDouble (gsl_matrix_complex_get (A, 2, 2))
	+ ComplexDouble (gsl_matrix_complex_get (C, 2, 2));
      gsl_matrix_complex_set (S, 2, 2, GslComplex (sum));
     }
  else
    {
   complex<double> sum;
      sum = 
	  ComplexDouble (gsl_matrix_complex_get (A, 0, 0))
	+ ComplexDouble (gsl_matrix_complex_get (C, 0, 0));
      gsl_matrix_complex_set (S, 0, 0, GslComplex (sum));
      sum = 
	  ComplexDouble (gsl_matrix_complex_get (A, 1, 0))
	+ ComplexDouble (gsl_matrix_complex_get (C, 1, 0));
      gsl_matrix_complex_set (S, 1, 0, GslComplex (sum));
      sum = 
	  ComplexDouble (gsl_matrix_complex_get (A, 2, 0))
	+ ComplexDouble (gsl_matrix_complex_get (C, 2, 0));
      gsl_matrix_complex_set (S, 2, 0, GslComplex (sum));
      sum = 
	  ComplexDouble (gsl_matrix_complex_get (A, 0, 1))
	- ComplexDouble (gsl_matrix_complex_get (C, 0, 1));
      gsl_matrix_complex_set (S, 0, 1, GslComplex (sum));
      sum = 
	  ComplexDouble (gsl_matrix_complex_get (A, 1, 1))
	- ComplexDouble (gsl_matrix_complex_get (C, 1, 1));
      gsl_matrix_complex_set (S, 1, 1, GslComplex (sum));
      sum = 
	  ComplexDouble (gsl_matrix_complex_get (A, 2, 1))
	- ComplexDouble (gsl_matrix_complex_get (C, 2, 1));
      gsl_matrix_complex_set (S, 2, 1, GslComplex (sum));
      sum = 
	  ComplexDouble (gsl_matrix_complex_get (A, 0, 2))
	- ComplexDouble (gsl_matrix_complex_get (C, 0, 2));
      gsl_matrix_complex_set (S, 0, 2, GslComplex (sum));
      sum = 
	  ComplexDouble (gsl_matrix_complex_get (A, 1, 2))
	- ComplexDouble (gsl_matrix_complex_get (C, 1, 2));
      gsl_matrix_complex_set (S, 1, 2, GslComplex (sum));
      sum = 
	  ComplexDouble (gsl_matrix_complex_get (A, 2, 2))
	- ComplexDouble (gsl_matrix_complex_get (C, 2, 2));
      gsl_matrix_complex_set (S, 2, 2, GslComplex (sum));
    }

  InvMatrix (B, S, EE, complex<double> (-1., 0.));
 
  // Evolve EE-matrix
  for (int l = 1; l < N-1; l++)
    {
      CalcABC (Q, X[l], h, A, B, C);
   
      MultMatrix (C, EE, S);
      AddMatrix  (B, S,  T);
      InvMatrix  (T, A,  EE, complex<double> (-1., 0.));
    }

  // Calculate DD-matrix
  CalcUV (Q, L, U, Up, Upp, V, I, UV, VU, VUp, V2, V2U, U2, V3);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      {
	complex<double> rhs = 
	  +                   ComplexDouble (gsl_matrix_complex_get (I,  i, j)) 
	  - (h/2.)      *     ComplexDouble (gsl_matrix_complex_get (V,  i, j))
	  + (h*h/6.)    * (   ComplexDouble (gsl_matrix_complex_get (V2, i, j)) 
			   +  ComplexDouble (gsl_matrix_complex_get (U,  i, j)))
	  - (h*h*h/24.) * (2.*ComplexDouble (gsl_matrix_complex_get (Up, i, j)) 
			   +  ComplexDouble (gsl_matrix_complex_get (UV, i, j))
			   +  ComplexDouble (gsl_matrix_complex_get (VU, i, j)) 
			   +  ComplexDouble (gsl_matrix_complex_get (V3, i, j)));
	gsl_matrix_complex_set (S, i, j, GslComplex (rhs));

	complex<double> rhs1 = 
	  +                    ComplexDouble (gsl_matrix_complex_get (I,   i, j))
	  -                    ComplexDouble (gsl_matrix_complex_get (EE,  i, j))
	  + (h*h/2.)      *    ComplexDouble (gsl_matrix_complex_get (U,   i, j))
	  - (h*h*h/6.)    * (  ComplexDouble (gsl_matrix_complex_get (Up,  i, j))
			     + ComplexDouble (gsl_matrix_complex_get (VU,  i, j)))
	  + (h*h*h*h/24.) * (  ComplexDouble (gsl_matrix_complex_get (Upp, i, j))
			     + ComplexDouble (gsl_matrix_complex_get (VUp, i, j))
			     + ComplexDouble (gsl_matrix_complex_get (V2U, i, j))
			     + ComplexDouble (gsl_matrix_complex_get (U2,  i, j)));
	gsl_matrix_complex_set (T, i, j, GslComplex (rhs1));
      }
  InvMatrix (S, T, DD, complex<double> (1./h, 0.));

  double p3, p4; complex<double> s5, s6;
  complex<double> Eplus, Eminu;

  // Calculate asymptotic series
  CalcSeries (Q, M, p3, p4, s5, s6, psi3, xi3, v3, psi4, xi4, v4, psi5, xi5, v5, psi6, xi6, v6);

  // Calculate asymptotic solution at x=L
  CalcAsymptotic (Q, L, M, p3, p4, s5, s6, 
		  psi3, xi3, v3, psi4, xi4, v4, psi5, xi5, v5, psi6, xi6, v6,
		  T3, T3p, T4, T4p, T5, T5p, T6, T6p);

  // Match asymptotic series at x=L
  Match (DD, T3, T3p, T4, T4p, T5, T5p, T6, T6p, Delta, Eplus, Eminu);

  // Calculate layer solution
  Psi[N-1] = Delta*T3[0] + T4[0] + Eplus*T5[0] + Eminu*T6[0];
  Xi [N-1] = Delta*T3[1] + T4[1] + Eplus*T5[1] + Eminu*T6[1];
  Ups[N-1] = Delta*T3[2] + T4[2] + Eplus*T5[2] + Eminu*T6[2];
  LoadVector (a, Psi[N-1], Xi[N-1], Ups[N-1]);
  MultMatVec (EE, a, b);
  UnloadVec  (b, Psi[N-2], Xi[N-2], Ups[N-2]);
  for (int i = N-2; i > 0; i--)
    {
      CalcABC    (Q, X[i], h, A, B, C);
      LoadVector (a, Psi[i+1], Xi[i+1], Ups[i+1]);
      MultMatVec (A, a, b);
      LoadVector (a, Psi[i], Xi[i], Ups[i]);
      MultMatVec (B, a, c);
      AddVector  (b, c, a);
      InvMatVec  (C, a, b, -1.);
      UnloadVec  (b, Psi[i-1], Xi[i-1], Ups[i-1]);
    }
  if (flag)
    {
      FILE *file;
      if (sign > 0)
	file = OpenFile ("Stage4/layere.out");
      else 
	file = OpenFile ("Stage4/layero.out");
      for (int i = 0; i < N; i++)
	fprintf (file, "%e %e %e %e %e %e %e\n", X[i], real(Psi[i]), imag(Psi[i]),
		 real(Xi[i]), imag(Xi[i]), real(Ups[i]), imag(Ups[i]));
      fclose (file);
    }

  // Clean up
  delete[] X;    delete[] a;   delete[] b;  delete[] c;    
  delete[] psi3; delete[] xi3; delete[] v3; delete[] psi4; delete[] xi4; delete[] v4;
  delete[] psi5; delete[] xi5; delete[] v5; delete[] psi6; delete[] xi6; delete[] v6;
  delete[] T3;   delete[] T3p; delete[] T4; delete[] T4p;
  delete[] T5;   delete[] T5p; delete[] T6; delete[] T6p;
  delete[] Psi;  delete[] Xi;  delete[] Ups;
  gsl_matrix_complex_free (EE);  gsl_matrix_complex_free (A);  gsl_matrix_complex_free (B);  
  gsl_matrix_complex_free (C);   gsl_matrix_complex_free (S);  gsl_matrix_complex_free (T);   
  gsl_matrix_complex_free (I);   gsl_matrix_complex_free (UV); gsl_matrix_complex_free (VU); 
  gsl_matrix_complex_free (VUp); gsl_matrix_complex_free (V2); gsl_matrix_complex_free (V2U);
  gsl_matrix_complex_free (U2);  gsl_matrix_complex_free (V3); gsl_matrix_complex_free (DD);
}

// ###################################################################
// Function to match asymptotic solution to finite-difference solution
// ###################################################################
void GGJ::Match (gsl_matrix_complex *DD, 
		 complex<double> *T3, complex<double> *T3p,
		 complex<double> *T4, complex<double> *T4p,
		 complex<double> *T5, complex<double> *T5p,
		 complex<double> *T6, complex<double> *T6p,
		 complex<double>& Delta, complex<double>& Eplus, complex<double>& Eminu)
{
  complex<double> *a3   = new complex<double>[3];
  complex<double> *a4   = new complex<double>[3];
  complex<double> *a5   = new complex<double>[3];
  complex<double> *a6   = new complex<double>[3];
  complex<double> *a56  = new complex<double>[3];
  complex<double> *a36  = new complex<double>[3];
  complex<double> *a35  = new complex<double>[3];

  // Calculate a3, a4, a5, a6
  MultMatVec (DD, T3, a3);
  for (int i = 0; i < 3; i++)
    a3[i] -= T3p[i];  
  MultMatVec (DD, T4, a4);
  for (int i = 0; i < 3; i++)
    a4[i] -= T4p[i];
  MultMatVec (DD, T5, a5);
  for (int i = 0; i < 3; i++)
    a5[i] -= T5p[i];  
  MultMatVec (DD, T6, a6);
  for (int i = 0; i < 3; i++)
    a6[i] -= T6p[i];

  // Calculate Delta, Eplus, Eminu
  CrossProduct (a5, a6, a56);
  CrossProduct (a3, a5, a35);
  CrossProduct (a3, a6, a36);
  Delta = - DotProduct (a4, a56) / DotProduct (a3, a56);
  Eplus = - DotProduct (a4, a36) / DotProduct (a5, a36);
  Eminu = - DotProduct (a4, a35) / DotProduct (a6, a35);

  delete[] a3;   delete[] a4;  delete[] a5; delete[] a6;   
  delete[] a56;  delete[] a35; delete[] a36;
}

// ##########################################
// Function to calculate asymptotic solutions
// ##########################################
void GGJ::CalcAsymptotic (complex<double> Q, double x, int M,
			  double p3, double p4,  complex<double>  s5,  complex<double>  s6, 
			  complex<double> *psi3, complex<double> *xi3, complex<double> *v3,
			  complex<double> *psi4, complex<double> *xi4, complex<double> *v4,
			  complex<double> *psi5, complex<double> *xi5, complex<double> *v5,
			  complex<double> *psi6, complex<double> *xi6, complex<double> *v6,
			  complex<double> *T3,   complex<double> *T3p,
			  complex<double> *T4,   complex<double> *T4p,
			  complex<double> *T5,   complex<double> *T5p,
			  complex<double> *T6,   complex<double> *T6p)
{
  // Calculate T3 and T3p
  T3 [0] = T3 [1] = T3 [2] = complex<double> (0., 0.);
  T3p[0] = T3p[1] = T3p[2] = complex<double> (0., 0.);
  for (int j = 0; j < M; j++)
    {
      double jj = double (j);
      T3 [0] += psi3[j]*pow(x, p3-2.*jj+1.);
      T3 [1] += xi3 [j]*pow(x, p3-2.*jj);
      T3 [2] += v3  [j]*pow(x, p3-2.*jj);
      T3p[0] += psi3[j]*(p3-2.*jj+1.)*pow(x, p3-2.*jj);
      T3p[1] += xi3 [j]*(p3-2.*jj   )*pow(x, p3-2.*jj-1.);
      T3p[2] += v3  [j]*(p3-2.*jj   )*pow(x, p3-2.*jj-1.);
    }

  // Calculate T4 and T4p
  T4 [0] = T4 [1] = T4 [2] = complex<double> (0., 0.);
  T4p[0] = T4p[1] = T4p[2] = complex<double> (0., 0.);
  for (int j = 0; j < M; j++)
    {
      double jj = double (j);
      T4 [0] += psi4[j]*pow(x, p4-2.*jj+1.);
      T4 [1] += xi4 [j]*pow(x, p4-2.*jj);
      T4 [2] += v4  [j]*pow(x, p4-2.*jj);
      T4p[0] += psi4[j]*(p4-2.*jj+1.)*pow(x, p4-2.*jj);
      T4p[1] += xi4 [j]*(p4-2.*jj   )*pow(x, p4-2.*jj-1.);
      T4p[2] += v4  [j]*(p4-2.*jj   )*pow(x, p4-2.*jj-1.);
    }

  // Calculate T5 and T5p
  T5 [0] = T5 [1] = T5 [2] = complex<double> (0., 0.);
  T5p[0] = T5p[1] = T5p[2] = complex<double> (0., 0.);
  for (int j = 0; j < M; j++)
    {
      double jj = double (j);
      T5 [0] += exp(-x*x/2./sqrt(Q))*psi5[j]*pow(x, s5-2.*jj-1.);
      T5 [1] += exp(-x*x/2./sqrt(Q))*xi5 [j]*pow(x, s5-2.*jj);
      T5 [2] += exp(-x*x/2./sqrt(Q))*v5  [j]*pow(x, s5-2.*jj);
      T5p[0] += exp(-x*x/2./sqrt(Q))*psi5[j]*(s5-2.*jj-1.)*pow(x, s5-2.*jj-2.) - (x/sqrt(Q))*T5[0];
      T5p[1] += exp(-x*x/2./sqrt(Q))*xi5 [j]*(s5-2.*jj   )*pow(x, s5-2.*jj-1.) - (x/sqrt(Q))*T5[1];
      T5p[2] += exp(-x*x/2./sqrt(Q))*v5  [j]*(s5-2.*jj   )*pow(x, s5-2.*jj-1.) - (x/sqrt(Q))*T5[2];
    }

  // Calculate T6 and T6p
  T6 [0] = T6 [1] = T6 [2] = complex<double> (0., 0.);
  T6p[0] = T6p[1] = T6p[2] = complex<double> (0., 0.);
  for (int j = 0; j < M; j++)
    {
      double jj = double (j);
      T6 [0] += exp(-x*x/2./sqrt(Q))*psi6[j]*pow(x, s6-2.*jj-1.);
      T6 [1] += exp(-x*x/2./sqrt(Q))*xi6 [j]*pow(x, s6-2.*jj);
      T6 [2] += exp(-x*x/2./sqrt(Q))*v6  [j]*pow(x, s6-2.*jj);
      T6p[0] += exp(-x*x/2./sqrt(Q))*psi6[j]*(s6-2.*jj-1.)*pow(x, s6-2.*jj-2.) - (x/sqrt(Q))*T6[0];
      T6p[1] += exp(-x*x/2./sqrt(Q))*xi6 [j]*(s6-2.*jj   )*pow(x, s6-2.*jj-1.) - (x/sqrt(Q))*T6[1];
      T6p[2] += exp(-x*x/2./sqrt(Q))*v6  [j]*(s6-2.*jj   )*pow(x, s6-2.*jj-1.) - (x/sqrt(Q))*T6[2];
    }
}

// #######################################
// Function to calculate asymptotic series
// #######################################
void GGJ::CalcSeries (complex<double> Q, int M,
		      double& p3, double &p4, complex<double>& s5,  complex<double>& s6,
		      complex<double> *psi3,  complex<double> *xi3, complex<double> *v3,
		      complex<double> *psi4,  complex<double> *xi4, complex<double> *v4,
		      complex<double> *psi5,  complex<double> *xi5, complex<double> *v5,
		      complex<double> *psi6,  complex<double> *xi6, complex<double> *v6)
{		     
  gsl_matrix_complex *W   = gsl_matrix_complex_alloc (3, 3);
  complex<double>    *rhs = new complex<double>[3];
  complex<double>    *a   = new complex<double>[3];

  // Calculate third asymptotic series
  p3 = - 0.5 + sqrt (0.25 - EF - H);
  psi3[0] = complex<double> (1., 0.);
  xi3 [0] = complex<double> (1., 0.);
  v3  [0] = complex<double> (1., 0.);

  complex<double> delta = - Q*Q*KFG - Q*Q*KEG - Q*Q*KH*(p3+1.);
  psi3[1] = (1./2./(1.-2.*p3)) * (- Q*Q*p3*(p3-1.) - delta * (EF-H*(p3-2.)));
  xi3 [1] = psi3[1] + (H/Q)*(p3+1.) + EF/Q;
  v3  [1] = psi3[1] + delta;

  for (int j = 2; j < M; j++)
    {
      double jj = double (j);
      delta = 
	- Q*Q*KFG*v3[j-1]
	- Q*Q*KEG*xi3[j-1] - Q*Q*KH*(p3+3.-2.*jj)*psi3[j-1]
	+ Q*(p3+4.-2.*jj)*(p3+3.-2.*jj)*v3[j-2];
      psi3[j] = (1./2./jj/(2.*jj-1.-2.*p3)) 
	* ( - Q*Q*(p3+2.-2.*jj)*(p3+1.-2.*jj)*xi3[j-1]
	    - (EF-H*(p3-2.*jj))*delta);
      xi3 [j] = psi3[j] + (H/Q)*(p3+3.-2.*jj)*psi3[j-1]
	+ (EF/Q)*v3[j-1] + Q*(p3+4.-2.*jj)*(p3+3.-2.*jj)*xi3[j-2];
      v3  [j] = psi3[j] + delta;
    }

  // Calculate fourth asymptotic series
  p4 = - 0.5 - sqrt (0.25 - EF - H);
  psi4[0] = complex<double> (1., 0.);
  xi4 [0] = complex<double> (1., 0.);
  v4  [0] = complex<double> (1., 0.);

  delta   = - Q*Q*KFG - Q*Q*KEG - Q*Q*KH*(p4+1.);
  psi4[1] = (1./2./(1.-2.*p4)) * (- Q*Q*p4*(p4-1.) - delta * (EF-H*(p4-2.)));
  xi4 [1] = psi4[1] + (H/Q)*(p4+1.) + EF/Q;
  v4  [1] = psi4[1] + delta;

  for (int j = 2; j < M; j++)
    {
      double jj = double (j);
      delta = 
	- Q*Q*KFG*v4[j-1]
	- Q*Q*KEG*xi4[j-1] - Q*Q*KH*(p4+3.-2.*jj)*psi4[j-1]
	+ Q*(p4+4.-2.*jj)*(p4+3.-2.*jj)*v4[j-2];
      psi4[j] = (1./2./jj/(2.*jj-1.-2.*p4)) 
	* ( - Q*Q*(p4+2.-2.*jj)*(p4+1.-2.*jj)*xi4[j-1]
	    - (EF-H*(p4-2.*jj))*delta);
      xi4 [j] = psi4[j] + (H/Q)*(p4+3.-2.*jj)*psi4[j-1]
	+ (EF/Q)*v4[j-1] + Q*(p4+4.-2.*jj)*(p4+3.-2.*jj)*xi4[j-2];
      v4  [j] = psi4[j] + delta;
    }

  // Calculate fifth asymptotic series
  complex<double> fa5 = 
    Q*Q*Q*((KFG-1.)*(KFG-1.)+KH*H*(KH*H+2.*(KFG+1.)))
    + 4.*((-KEG-1.)*(EF+H*H)+H*H);
  s5  = - 0.5 - 0.25 * pow(Q,1.5) * (1.+KFG+KH*H) + 0.25*sqrt(fa5);
  psi5[0] = 1.;
  xi5 [0] = (- EF/Q + H*(Q-H/sqrt(Q))/sqrt(Q))      /Q/(EF+H*(2.*s5+1.));
  v5  [0] = (- sqrt(Q)*(2.*s5+1.) - Q*(Q-H/sqrt(Q)))/Q/(EF+H*(2.*s5+1.));

  for (int j = 1; j < M; j++)
    {
      double (jj) = double (j);
      
      gsl_complex z;

      z = GslComplex (1./Q);
      gsl_matrix_complex_set (W, 0, 0, z);

      z = GslComplex (Q);
      gsl_matrix_complex_set (W, 0, 1, z);

      z = GslComplex(H/sqrt(Q));
      gsl_matrix_complex_set (W, 0, 2, z);

      z = GslComplex(Q-H/sqrt(Q));
      gsl_matrix_complex_set (W, 1, 0, z);

      z = GslComplex (-pow(Q,1.5)*(2.*s5+1.-4.*jj));
      gsl_matrix_complex_set (W, 1, 1, z);
 
      z = GslComplex (EF);
      gsl_matrix_complex_set (W, 1, 2, z);

      z = GslComplex(1.+pow(Q,1.5)*KH);
      gsl_matrix_complex_set (W, 2, 0, z);

      z = GslComplex (-Q*Q*KEG);
      gsl_matrix_complex_set (W, 2, 1, z);

      z = GslComplex (-Q*Q*KFG-sqrt(Q)*(2.*s5+1.-4.*jj));
      gsl_matrix_complex_set (W, 2, 2, z);

      rhs[0] = 
	+ (Q+(2.*s5+3.-4.*jj)/sqrt(Q))*psi5[j-1]
	+ H*(s5+2.-2.*jj)*v5[j-1];
      if (j > 1)
	rhs[0] -= (s5+3.-2.*jj)*(s5+2.-2.*jj)*psi5[j-2];
      rhs[1] = 
	- H*(s5+1.-2.*jj)*psi5[j-1]
	- Q*Q*(s5+2.-2.*jj)*(s5+1.-2.*jj)*xi5[j-1];
      rhs[2] = 
	- Q*(s5+2.-2.*jj)*(s5+1.-2.*jj)*v5[j-1]
	+ Q*Q*KH*(s5+1.-2.*jj)*psi5[j-1];

      InvMatVec (W, rhs, a, complex<double>(1., 0.));
      psi5[j] = a[0];
      xi5 [j] = a[1];
      v5  [j] = a[2];
    }

  // Calculate sixth asymptotic series
  complex<double> fa6 = 
    Q*Q*Q*((KFG-1.)*(KFG-1.)+KH*H*(KH*H+2.*(KFG+1.)))
    + 4.*((-KEG-1.)*(EF+H*H)+H*H);
  s6  = - 0.5 - 0.25 * pow(Q,1.5) * (1.+KFG+KH*H) - 0.25*sqrt(fa6);
  psi6[0] = 1.;
  xi6 [0] = (- EF/Q + H*(Q-H/sqrt(Q))/sqrt(Q))      /Q/(EF+H*(2.*s6+1.));
  v6  [0] = (- sqrt(Q)*(2.*s6+1.) - Q*(Q-H/sqrt(Q)))/Q/(EF+H*(2.*s6+1.));

  for (int j = 1; j < M; j++)
    {
      double (jj) = double (j);
            
      gsl_complex z;

      z = GslComplex (1./Q);
      gsl_matrix_complex_set (W, 0, 0, z);

      z = GslComplex (Q);
      gsl_matrix_complex_set (W, 0, 1, z);

      z = GslComplex(H/sqrt(Q));
      gsl_matrix_complex_set (W, 0, 2, z);

      z = GslComplex(Q-H/sqrt(Q));
      gsl_matrix_complex_set (W, 1, 0, z);

      z = GslComplex (-pow(Q,1.5)*(2.*s6+1.-4.*jj));
      gsl_matrix_complex_set (W, 1, 1, z);
 
      z = GslComplex (EF);
      gsl_matrix_complex_set (W, 1, 2, z);

      z = GslComplex(1.+pow(Q,1.5)*KH);
      gsl_matrix_complex_set (W, 2, 0, z);

      z = GslComplex (-Q*Q*KEG);
      gsl_matrix_complex_set (W, 2, 1, z);

      z = GslComplex (-Q*Q*KFG-sqrt(Q)*(2.*s6+1.-4.*jj));
      gsl_matrix_complex_set (W, 2, 2, z);

      rhs[0] = 
	+ (Q+(2.*s6+3.-4.*jj)/sqrt(Q))*psi6[j-1]
	+ H*(s6+2.-2.*jj)*v6[j-1];
      if (j > 1)
	rhs[0] -= (s6+3.-2.*jj)*(s6+2.-2.*jj)*psi6[j-2];
      rhs[1] = 
	- H*(s6+1.-2.*jj)*psi6[j-1]
	- Q*Q*(s6+2.-2.*jj)*(s6+1.-2.*jj)*xi6[j-1];
      rhs[2] = 
	- Q*(s6+2.-2.*jj)*(s6+1.-2.*jj)*v6[j-1]
	+ Q*Q*KH*(s6+1.-2.*jj)*psi6[j-1];

      InvMatVec (W, rhs, a, complex<double> (1., 0.));
      psi6[j] = a[0];
      xi6 [j] = a[1];
      v6  [j] = a[2];
    }

  gsl_matrix_complex_free (W);   
  delete[] rhs; delete[] a;
}  

// ##########################################
// Function to calculate A, B, and C matrices
// ##########################################
void GGJ::CalcABC (complex<double> Q, double x, double h, 
		   gsl_matrix_complex *A, gsl_matrix_complex *B, gsl_matrix_complex *C)
{
  // Allocate memory
  gsl_matrix_complex *U   = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *Up  = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *Upp = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *V   = gsl_matrix_complex_alloc (3, 3);

  gsl_matrix_complex *I   = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *UV  = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *VU  = gsl_matrix_complex_alloc (3, 3); 
  gsl_matrix_complex *VUp = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *V2  = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *V2U = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *U2  = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *V3  = gsl_matrix_complex_alloc (3, 3);

  gsl_matrix_complex *M1  = gsl_matrix_complex_alloc (3, 3);
  gsl_matrix_complex *M2  = gsl_matrix_complex_alloc (3, 3);

  // Calculate U, Up, Upp, V, etc. matrices
  CalcUV (Q, x, U, Up, Upp, V, I, UV, VU, VUp, V2, V2U, U2, V3);
  
  // Set elements of M1, M2
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      {
	complex<double> rhs1 = 
	  + ComplexDouble (gsl_matrix_complex_get (Upp, i, j)) 
	  - ComplexDouble (gsl_matrix_complex_get (VUp, i, j)) 
	  - ComplexDouble (gsl_matrix_complex_get (V2U, i, j))
	  + ComplexDouble (gsl_matrix_complex_get (U2,  i, j));
	complex<double> rhs2 = 
	  + 2.*ComplexDouble (gsl_matrix_complex_get (Up, i, j))
	  +    ComplexDouble (gsl_matrix_complex_get (UV, i, j)) 
	  -    ComplexDouble (gsl_matrix_complex_get (VU, i, j))
	  -    ComplexDouble (gsl_matrix_complex_get (V3, i, j));
	gsl_matrix_complex_set (M1, i, j, GslComplex (rhs1));
	gsl_matrix_complex_set (M2, i, j, GslComplex (rhs2));
      }

  // Set elements of A, B, C
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      {
	complex<double> rhsa = 
	  +               ComplexDouble (gsl_matrix_complex_get (I,  i, j))
	  - (h/2.)      * ComplexDouble (gsl_matrix_complex_get (V,  i, j))
	  - (h*h*h/24.) * ComplexDouble (gsl_matrix_complex_get (M2, i, j));
	complex<double> rhsb = 
	  - 2.            * ComplexDouble (gsl_matrix_complex_get (I,  i, j))
	  - h*h           * ComplexDouble (gsl_matrix_complex_get (U,  i, j))
	  - (h*h*h*h/12.) * ComplexDouble (gsl_matrix_complex_get (M1, i, j));
	complex<double> rhsc = 
	  +               ComplexDouble (gsl_matrix_complex_get (I,  i, j))
	  + (h/2.)      * ComplexDouble (gsl_matrix_complex_get (V,  i, j))
	  + (h*h*h/24.) * ComplexDouble (gsl_matrix_complex_get (M2, i, j));
	  
	gsl_matrix_complex_set (A, i, j, GslComplex (rhsa));
	gsl_matrix_complex_set (B, i, j, GslComplex (rhsb));
	gsl_matrix_complex_set (C, i, j, GslComplex (rhsc));
      }

  // Clean up
  gsl_matrix_complex_free (U);  gsl_matrix_complex_free (Up);  gsl_matrix_complex_free (Upp);
  gsl_matrix_complex_free (V);  gsl_matrix_complex_free (I);   gsl_matrix_complex_free (UV); 
  gsl_matrix_complex_free (VU); gsl_matrix_complex_free (VUp); gsl_matrix_complex_free (V2U); 
  gsl_matrix_complex_free (U2); gsl_matrix_complex_free (V2);  gsl_matrix_complex_free (V3);
  gsl_matrix_complex_free (M1); gsl_matrix_complex_free (M2);
}

// ################################################
// Function to calculate U, Up, Upp, and V matrices
// ################################################
void GGJ::CalcUV (complex<double> Q, double x, 
		  gsl_matrix_complex *U,   gsl_matrix_complex *Up,  gsl_matrix_complex *Upp, 
		  gsl_matrix_complex *V,   gsl_matrix_complex *I,   gsl_matrix_complex *UV,  
		  gsl_matrix_complex *VU,  gsl_matrix_complex *VUp, gsl_matrix_complex *V2, 
		  gsl_matrix_complex *V2U, gsl_matrix_complex *U2,  gsl_matrix_complex *V3)
{
  // Initialize matrix elements to zero
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      {
	gsl_matrix_complex_set (U,   i, j, gsl_complex_rect(0., 0.));
	gsl_matrix_complex_set (Up,  i, j, gsl_complex_rect(0., 0.));
	gsl_matrix_complex_set (Upp, i, j, gsl_complex_rect(0., 0.));
	gsl_matrix_complex_set (V,   i, j, gsl_complex_rect(0., 0.));
	gsl_matrix_complex_set (I,   i, j, gsl_complex_rect(0., 0.));
      }

  // Set elements of U
  gsl_matrix_complex_set (U, 0, 0, GslComplex (  Q)            );
  gsl_matrix_complex_set (U, 0, 1, GslComplex (- x*Q)          );
  gsl_matrix_complex_set (U, 1, 0, GslComplex (- x/Q)          );
  gsl_matrix_complex_set (U, 1, 1, GslComplex (  x*x/Q)        );
  gsl_matrix_complex_set (U, 1, 2, GslComplex (- EF/Q/Q)       );
  gsl_matrix_complex_set (U, 2, 0, GslComplex (- x/Q)          );
  gsl_matrix_complex_set (U, 2, 1, GslComplex (  KEG*Q)        );
  gsl_matrix_complex_set (U, 2, 2, GslComplex (  x*x/Q + KFG*Q));

  // Set elements of Up
  gsl_matrix_complex_set (Up, 0, 1, GslComplex (- Q)     );
  gsl_matrix_complex_set (Up, 1, 0, GslComplex (- 1./Q)  );
  gsl_matrix_complex_set (Up, 1, 1, GslComplex (  2.*x/Q));
  gsl_matrix_complex_set (Up, 2, 0, GslComplex (- 1./Q)  );
  gsl_matrix_complex_set (Up, 2, 2, GslComplex (  2.*x/Q));

  // Set elements of Upp
  gsl_matrix_complex_set (Upp, 1, 1, GslComplex (2./Q));
  gsl_matrix_complex_set (Upp, 2, 2, GslComplex (2./Q));

  // Set elements of V
  gsl_matrix_complex_set (V, 0, 2, GslComplex (  H)    );
  gsl_matrix_complex_set (V, 1, 0, GslComplex (- H/Q/Q));
  gsl_matrix_complex_set (V, 2, 0, GslComplex (  KH*Q) );

  // Set elements of I
  gsl_matrix_complex_set (I, 0, 0, GslComplex (1.));
  gsl_matrix_complex_set (I, 1, 1, GslComplex (1.));
  gsl_matrix_complex_set (I, 2, 2, GslComplex (1.));

  // Set elements of UV, VU, etc
  MultMatrix (U,  V,  UV);
  MultMatrix (V,  U,  VU);
  MultMatrix (V,  Up, VUp);
  MultMatrix (V,  V,  V2);
  MultMatrix (V2, U,  V2U);
  MultMatrix (U,  U,  U2);
  MultMatrix (V2, V,  V3);
}

// #######################################
// Function to add 3x3 matrices: C = A + B
// #######################################
void GGJ::AddMatrix (gsl_matrix_complex *A, gsl_matrix_complex *B, gsl_matrix_complex *C)
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      {	
	complex<double> sum = 
	    ComplexDouble (gsl_matrix_complex_get (A, i, j)) 
	  + ComplexDouble (gsl_matrix_complex_get (B, i, j));
	gsl_matrix_complex_set (C, i, j, GslComplex (sum));
      }
}

// ######################################
// Function to add 3x1 vectors: c = a + a
// ######################################
void GGJ::AddVector (complex<double> *a, complex<double> *b, complex<double> *c)
{
  for (int i = 0; i < 3; i++)
    c[i] = a[i] + b[i];
}

// ##########################################
// Function to multiply 3x3 matrices: C = A*B
// ##########################################
void GGJ::MultMatrix (gsl_matrix_complex *A, gsl_matrix_complex *B, gsl_matrix_complex *C)
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      {
	complex<double> sum = complex<double> (0., 0.);
	for (int k = 0; k < 3; k++)
	  sum += 
	      ComplexDouble (gsl_matrix_complex_get (A, i, k)) 
	    * ComplexDouble (gsl_matrix_complex_get (B, k, j));
	gsl_matrix_complex_set (C, i, j, GslComplex (sum));
      }
}

// #######################################################
// Function to multiply 3x3 matrix and 3x1 vector: c = A*b
// #######################################################
void GGJ::MultMatVec (gsl_matrix_complex *A, complex<double> *b, complex<double> *c)
{
  for (int i = 0; i < 3; i++)
    {
      complex<double> sum = complex<double> (0., 0.);
      for (int k = 0; k < 3; k++)
	sum += ComplexDouble (gsl_matrix_complex_get (A, i, k)) * b[k];
      c[i] = sum;
    }
}

// #######################################################
// Function to invert 3x3 matrix equation: C = g*A^(-1)*B
// #######################################################
void GGJ::InvMatrix (gsl_matrix_complex *A, gsl_matrix_complex *B, gsl_matrix_complex *C, complex<double> g)
{
  complex<double> *b = new complex<double>[3];
  complex<double> *c = new complex<double>[3];

  for (int j = 0; j < 3; j++)
    {
      for (int i = 0; i < 3; i++)
	b[i] = ComplexDouble (gsl_matrix_complex_get (B, i, j));

      InvMatVec (A, b, c, g);

      for (int i = 0; i < 3; i++)
	gsl_matrix_complex_set (C, i, j, GslComplex (c[i]));
    }

  delete[] b;
  delete[] c;
}

// #############################################################
// Function to invert 3x3 matrix/vector equation: c = g*A^(-1)*b
// #############################################################
void GGJ::InvMatVec (gsl_matrix_complex *A, complex<double> *b, complex<double> *c, complex<double> g)
{
  gsl_matrix_complex  *LU = gsl_matrix_complex_alloc (3, 3);
  gsl_vector_complex  *bb = gsl_vector_complex_alloc (3);
  gsl_vector_complex  *cc = gsl_vector_complex_alloc (3);
  gsl_permutation     *p  = gsl_permutation_alloc    (3);
  int s; 

  for (int i = 0; i < 3; i++)
    gsl_vector_complex_set (bb, i, GslComplex (b[i]));

  gsl_matrix_complex_memcpy    (LU, A);
  gsl_linalg_complex_LU_decomp (LU, p, &s);
  gsl_linalg_complex_LU_solve  (LU, p, bb, cc);

  for (int i = 0; i < 3; i++)
    c[i] = g*ComplexDouble (gsl_vector_complex_get (cc, i));

  gsl_matrix_complex_free (LU);
  gsl_vector_complex_free (bb);
  gsl_vector_complex_free (cc);
  gsl_permutation_free    (p);
}

// ###########################################################
// Function to calculate cross product of 3x1 vectors: c = axb
// ###########################################################
void GGJ::CrossProduct (complex<double> *a, complex<double> *b, complex<double> *c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

// #########################################################
// Function to calculate dot product of 3x1 vectors: c = a.b
// #########################################################
complex<double> GGJ::DotProduct (complex<double> *a, complex<double> *b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// ###########################
// Function to load 3x1 vector
// ###########################
void GGJ::LoadVector (complex<double> *a, complex<double> a0, complex<double> a1, complex<double> a2)
{
  a[0] = a0;
  a[1] = a1;
  a[2] = a2;
}

// #############################
// Function to unload 3x1 vector
// #############################
void GGJ::UnloadVec (complex<double> *a, complex<double>& a0, complex<double>& a1, complex<double>& a2)
{
  a0 = a[0];
  a1 = a[1];
  a2 = a[2];
}

// ####################################################
// Function to convert complex<double> into gsl_complex
// ####################################################
gsl_complex GGJ::GslComplex (complex<double> z)
{
  return gsl_complex_rect (real(z), imag(z));
}

// ###########################################
// Function to convert double into gsl_complex
// ###########################################
gsl_complex GGJ::GslComplex (double x)
{
  return gsl_complex_rect (x, 0.);
}

// ####################################################
// Function to convert gsl_complex into complex<double>
// ####################################################
complex<double> GGJ::ComplexDouble (gsl_complex z)
{
  return complex<double> (GSL_REAL(z), GSL_IMAG(z));
}

// #############################################
// Function to return modulus of complex<double>
// #############################################
double GGJ::hypot (complex<double> z)
{
  return sqrt(real(z)*real(z) + imag(z)*imag(z));
}

// ################################
// Target function for zero finding
// ################################
double GGJ::Feval (double x)
{
  double Deltar, Deltai;
  GetDeltae (0., x, Deltar, Deltai);
  
  return Deltai;
}

// ###################################################################
//  Routine to find approximate root of F(x) = 0 using Ridder's method. 
//  Search takes place in interval (0., 1.). 
//  Interval is chopped into nint equal segments.
// ###################################################################
double GGJ::RootFind (double x1, double x2)
{
  double F1, F2 = 0., root;

  // Chop search interval into nint segments  
  for (int i = 0; i < nint; i++)
    {
      double x1_seg = x1 + (x2 - x1) * double (i)     / double (nint);
      double x2_seg = x1 + (x2 - x1) * double (i + 1) / double (nint);
      
      if (i == 0) 
	F1 = Feval (x1_seg);
      else 
	F1 = F2;
      F2 = Feval (x2_seg);
      //printf ("%e %e %e %e\n", x1_seg, F1, x2_seg, F2);
      
      // Call Ridder's method for segment containing zero
      if (F1 * F2 < 0.)
	{
	  Ridder (x1_seg, x2_seg, F1, F2, root);
	  break;
	}
    }

  return root;
}

// #############################################
// Ridder's method for finding root of F(x) = 0.
// #############################################
void GGJ::Ridder (double x1, double x2, double F1, double F2, double& x)
{
  // Iteration loop  
  x = x2; double xold, Fx; int iter = 0;
  do 
    {              
      // Calculate F(x3), where x3 is midpoint of current interval 
      double x3 = (x1 + x2) / 2.;    
      double F3 = Feval (x3);
      
      // Iterate x using Ridder's method 
      xold = x;           
      x = x3 - (x3 - x1) * (F2 - F1) * F3 /
	(sqrt (F3 * F3 - F1 * F2) * fabs (F2 - F1));
      Fx = Feval (x);
       
      // Make new value of x upper/lower bound of refined search 
      // interval, as appropriate 
      if (Fx * F1 < 0.) 
	{  
	  x2 = x;           
	  F2 = Fx; 
	}
      else 
	{
	  x1 = x;
	  F1 = Fx; 
	}
      //printf ("%d %e %e\n", iter, x, Fx);
      iter++;
    } 
  // Iterate until absolute change in x falls below Eta
  while (fabs (x - xold) > Eta && fabs(Fx) > Eta && iter < Maxiter); 
}

// #####################################
// Function to open new file for writing
// #####################################
FILE* GGJ::OpenFile (char* filename)
{
  FILE *file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("GGJ::OpenFile: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

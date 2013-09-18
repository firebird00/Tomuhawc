// Profile.cpp

#include "Flux.h"

// ######################################
// Function to evaluate profile functions
// ######################################
void Flux::Profile ()
{
  double *QX  = new double [J];
  double *P1  = new double [J];
  double *P2  = new double [J];
  double *P3  = new double [J];
  double *P4  = new double [J];
  double *P5  = new double [J];
  double *P6  = new double [J];
  double *P1P = new double [J];
  double *P2P = new double [J];
  double *P3P = new double [J];

  for (int j = 0; j < J; j++)
    {
      P1 [j] = epsa*epsa*rP[j]*rP[j];
      P2 [j] = GPP[j];
      P3 [j] = QGP[j] * PPP[j];
      P5 [j] = P1[j] /QGP[j]/QGP[j];
      P6 [j] = Gamma * PP[j];
      P1P[j] = 2.*epsa*epsa*rP[j];
    }

  for (int j = 0; j < J; j++)
    {
      double qpval, qgpval, p2pval, p3pval;

      qpval  = Interpolate (J, rP, QP,  rP[j], 1);
      qgpval = Interpolate (J, rP, QGP, rP[j], 1);
      p2pval = Interpolate (J, rP, P2,  rP[j], 1);
      p3pval = Interpolate (J, rP, P3,  rP[j], 1);
      
      QX [j] = qpval;
      P4 [j] = - rP[j] * qgpval / QGP[j];
      P2P[j] = p2pval;
      P3P[j] = p3pval;
    }	

  FILE *file = OpenFile ("Stage2/profile.out");
  for (int j = 0; j < J; j++)
    fprintf (file, "%e %e %e %e %e %e %e %e %e %e %e %e\n", 
	     rP[j], QP[j], QX[j], P1[j], P2[j], P3[j], P4[j], P5[j], P6[j], P1P[j], P2P[j], P3P[j]);
  double one = 1.;
  fprintf (file, "%e %e %e %e %e %e %e %e %e %e %e %e\n", 
	   one, 
	   QA,
	   Extrapolate1 (J, rP, QX,  one), 
	   Extrapolate1 (J, rP, P1,  one), 
	   Extrapolate1 (J, rP, P2,  one), 
	   Extrapolate1 (J, rP, P3,  one), 
	   Extrapolate1 (J, rP, P4,  one), 
	   Extrapolate1 (J, rP, P5,  one), 
	   Extrapolate1 (J, rP, P6,  one),
	   Extrapolate1 (J, rP, P1P, one),
	   Extrapolate1 (J, rP, P2P, one),
	   Extrapolate1 (J, rP, P3P, one));
  fclose (file);

  delete[] QX; 
  delete[] P1;  delete[] P2;  delete[] P3; 
  delete[] P4;  delete[] P5;  delete[] P6; 
  delete[] P1P; delete[] P2P; delete[] P3P;
}


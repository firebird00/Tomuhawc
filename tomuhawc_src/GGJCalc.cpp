// GGJCalc.cpp

#include "Tomuhawc.h"

// #####################################################################
// Function to calculate Glasser, Green, Johnson layer theory parameters
//
// R .. flux surface label
//
// #####################################################################
void Thawc::GGJCalc (double R, double& EF, double& HH, double& ET, double& FT, 
		     double& HT, double& FA, double& FR)
{
  // Get profile data
  double Q   = gsl_spline_eval (Sq,     R, Aq);
  double QP  = gsl_spline_eval (Sqp,    R, Aqp);
  double P1  = gsl_spline_eval (Sprof1, R, Aprof1);
  double P3  = gsl_spline_eval (Sprof3, R, Aprof3);
  double P4  = gsl_spline_eval (Sprof4, R, Aprof4);
  double P5  = gsl_spline_eval (Sprof5, R, Aprof5);
  double P6  = gsl_spline_eval (Sprof6, R, Aprof6);

  // Get metric data
  double MET1, MET1P, MET3, MET4, MET5, MET8, MET9;
  MET1  = gsl_spline_eval (SM1[0],  R, AM1[0] );
  MET3  = gsl_spline_eval (SM3[0],  R, AM3[0] );
  MET4  = gsl_spline_eval (SM4[0],  R, AM4[0] );
  MET5  = gsl_spline_eval (SM5[0],  R, AM5[0] );  
  MET8  = gsl_spline_eval (SM8,     R, AM8    );
  MET9  = gsl_spline_eval (SM9,     R, AM9    );
  MET1P = gsl_spline_eval (SM1P[0], R, AM1P[0]);

  // Calculate common quantities
  double S = R*QP/Q;

  // Calculate GGJ parameters
  double facx = 1.   + MET8*P1/Q/Q;
  double facc = MET3 +      P1/Q/Q;
  double facy = MET9 * facx - MET1*MET1;

  FR = facx /facc;
  FA = facy /P5/FR;

  double GG = facc * MET1 * FR/FA/P6; 

  EF = (P3/S/S)*(facc * (P4*MET1 + P3*MET5 - R*MET1P) - P3*MET4*MET4 + S*MET1/FR);
  HH = (P3/S/S)*(S*MET4 - S*MET1/FR);

  ET = - GG + (facc * (P4*MET1 - R*MET1P) + S*MET1/FR) *FR/FA/P3/P5;
  FT =   GG + (facc*MET5 - MET4*MET4)                  *FR/FA   /P5;
  HT = (S*MET4 - S*MET1/FR)                            *FR/FA/P3/P5;
}


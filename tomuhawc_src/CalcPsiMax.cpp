// CalcYmax.cpp

#include "Tomuhawc.h"

// #################################
// Function to calculate maximum Psi
// #################################
double Thawc::CalcPsiMax (double Y[])
{
  double ymax = 0.;
  for (int i = 0; i < dim+vac; i++)
    for (int k = 0; k < dim; k++)
      ymax = (fabs(Y[dim1*i+k]) > ymax) ? fabs(Y[dim1*i+k]) : ymax;

  return log10(ymax);
}

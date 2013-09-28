// Segment1.cpp

#include "Tomuhawc.h"

// ####################################################################################
// Function to integrate single solution vector to r = rx
// via adaptive step method
//
// r                   .. flux-surface label
// rx                  .. integrate solution vector until r = rx
// Y[dimN]             .. solution vector
//
//   i = 0,dim-1         .. independent solutions launched from axis
//   i = dim,dim+vac-1   .. small solutions launched from rational surfaces
//
//      Y[dim1*i     +0] .. y[dim1*i     +dim-1]  - psi
//      Y[dim1*i+dim +0] .. y[dim1*i+dim +dim-1]  - Z
// flag                .. truncation error flag
// interactive         .. set if program in iteractive mode
// _adap               .. name of adaptive step monitor file
// _soln               .. name of solution file
//
// side                .. number of sideband harmonics
// vac                 .. total number of rational surfaces (including vacuum surfaces)
// num                 .. number of rational surfaces in plasma
// dim  = vac + 2*side .. dimension of coupling matrices
// dim1 = 2*dim        .. dimension of single solution vector 
//
// flag = 0            .. truncation error is absolute
// flag = 1            .. truncation error is relative
// flag = 2            .. truncation error is mixed
//
// NOTE: Core integration always uses absolute truncation error
//
// meth = 0            .. 4th order (classical) Runge-Kutta method
// meth = 1            .. Embedded Runge-Kutta-Fehlberg (4,5) method
// meth = 2            .. Embedded Runge-Kutta Cash-Karp (4,5) method
// meth = 3            .. Embedded Runge-Kutta Prince-Dormand (8,9) method
// meth = 4            .. Implicit 4th order Runge-Kutta at Gaussian points
// meth = 5            .. M=2 implicit Gear method
//
// ####################################################################################
void Thawc::Segment1 (double &r, double rx, double Y[], 
		      int flag, int interactive,
		      char *_adap, char *_soln)
{
  gsl_odeiv_system sys = {pRhs, NULL, dimN, this};
  const gsl_odeiv_step_type *T;
  gsl_odeiv_step    *s;
  gsl_odeiv_control *c;
  gsl_odeiv_evolve  *e;
  
  if (meth == 0)
    T = gsl_odeiv_step_rk4;
  else if (meth == 1)
    T = gsl_odeiv_step_rkf45;
  else if (meth == 2)
    T = gsl_odeiv_step_rkck;
  else if (meth == 3)
    T = gsl_odeiv_step_rk8pd;
  else if (meth == 4)
    T = gsl_odeiv_step_rk4imp;
  else if (meth == 5)
    T = gsl_odeiv_step_gear2;
  
  s = gsl_odeiv_step_alloc   (T, dimN);
  e = gsl_odeiv_evolve_alloc (dimN);
  
  if (flag == 0)
    c = gsl_odeiv_control_y_new (acc, 0.);
  else if (flag == 1)
    c = gsl_odeiv_control_y_new (0., acc);
  else if (flag == 2)
    c = gsl_odeiv_control_y_new (acc/2., acc/2.);
  
  FILE *adap, *soln; 
  if (interactive) 
    {
      adap = OpenFilea (_adap);
      soln = OpenFilea (_soln);
    }
  double h = h0;
  int count = 1, status; 
  while (r < rx)
    {
      status = GSL_SUCCESS;
      status = gsl_odeiv_evolve_apply (e, c, s, &sys, &r, rx, &h, Y);
      step++;
      
      if (status != GSL_SUCCESS)
	{
	  printf ("Segment1: status != GSL_SUCCESS\n");
	  exit (1);
	}
      
      if (count%skip == 0)
	if (interactive) 
	  {
	    fprintf (adap, "%e %e %e\n", r, log10(h), CalcPsiMax(Y));

	    fprintf (soln, "%e", r);
	    for (int i = 0; i < dimN; i++)
	      fprintf (soln, " %20.15e", Y[i]);
	    fprintf (soln, "\n");
	    fflush (adap); fflush (soln);
	  }
     
      count++;
    }
  if (interactive)
    {
      fclose (adap); 
      fclose (soln);
    }

  gsl_odeiv_evolve_free  (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free    (s);
}

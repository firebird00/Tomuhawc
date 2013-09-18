// CalcrP.cpp

#include "Flux.h"

// ##################################
// Function to calculate r(P) profile
// ##################################
void Flux::CalcrP ()
{
  gsl_odeiv_system           sys2 = {pRhs2, NULL, 1, this};
  const gsl_odeiv_step_type *T    = gsl_odeiv_step_rk8pd;
  gsl_odeiv_step            *s    = gsl_odeiv_step_alloc (T, 1);
  gsl_odeiv_control         *c    = gsl_odeiv_control_y_new (acc/2., acc/2.);
  gsl_odeiv_evolve          *e    = gsl_odeiv_evolve_alloc (1);
  double                    *y    = new double [1]; 
  double                     r, h;

  // Calculate epsa
  r    = 0.;
  h    = h0;
  y[0] = 0.;
  while (r < P[0])
    {
      int status = gsl_odeiv_evolve_apply (e, c, s, &sys2, &r, P[0], &h, y);
      
      if (status != GSL_SUCCESS)
	{
	  printf ("Profile: status != GSL_SUCCESS\n");
	  exit (1);
	}
    }
  epsa = sqrt (y[0]);

  // Calculate epsb
  r    = psib;
  h    = h0;
  y[0] = 0.;
  while (r < P[0])
    {
      int status = gsl_odeiv_evolve_apply (e, c, s, &sys2, &r, P[0], &h, y);
      
      if (status != GSL_SUCCESS)
	{
	  printf ("Profile: status != GSL_SUCCESS\n");
	  exit (1);
	}
    }
  epsb = sqrt (y[0]);

  // Calculate r[Psi]
  for (int j = 1; j < J; j++)
    {
      r    = P[j];
      h    = h0;
      y[0] = 0.;
      
      while (r < P[0])
	{
	  int status = gsl_odeiv_evolve_apply (e, c, s, &sys2, &r, P[0], &h, y);
	  
	  if (status != GSL_SUCCESS)
	    {
	      printf ("Profile: status != GSL_SUCCESS\n");
	      exit (1);
	    }
	}
      rP[j] = sqrt (y[0]) /epsa;
    }

  gsl_odeiv_evolve_free  (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free    (s);
  delete[]                y;
}

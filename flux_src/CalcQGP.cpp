// CalcQGP.cpp

#include "Flux.h"

// #######################################
// Function to calculate q(P)/g(P) profile
// #######################################
void Flux::CalcQGP ()
{
  gsl_odeiv_system           sys1 = {pRhs1, NULL, 3, this};
  const gsl_odeiv_step_type *T    = gsl_odeiv_step_rk8pd;
  gsl_odeiv_step            *s    = gsl_odeiv_step_alloc (T, 3);
  gsl_odeiv_control         *c    = gsl_odeiv_control_y_new (acc/2., acc/2.);
  gsl_odeiv_evolve          *e    = gsl_odeiv_evolve_alloc (3);
  double                    *y    = new double [3]; 
  double                     r, h, y1old, y2old;
  
  for (int j = 1; j < J; j++)
    {
      r    = 0.;
      h    = h0;
      y[0] = RP[j];
      y[1] = Z[jc];
      y[2] = 0.;
      
      do
	{
	  y1old = y[1];
	  y2old = y[2];
	  int status = gsl_odeiv_evolve_apply (e, c, s, &sys1, &r, 100., &h, y);

	  if (status != GSL_SUCCESS)
	    {
	      printf ("Profile: status != GSL_SUCCESS\n");
	      exit (1);
	    }
	}
      while (y[1] > 0.);

      QGP[j] = (y1old*y[2] - y[1]*y2old) / (y1old - y[1]); 
      QP [j] = QGP[j]*GP[j];
     }

  gsl_odeiv_evolve_free  (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free    (s);
  delete[]                y;
}

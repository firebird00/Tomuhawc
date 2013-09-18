// Calcst.cpp

#include "Flux.h"

// ####################################
// Function to calculate straight angle
// ####################################
void Flux::Calctst ()
{
  gsl_odeiv_system           sys3 = {pRhs3, NULL, 2, this};
  const gsl_odeiv_step_type *T    = gsl_odeiv_step_rk8pd;
  gsl_odeiv_step            *s    = gsl_odeiv_step_alloc (T, 2);
  gsl_odeiv_control         *c    = gsl_odeiv_control_y_new (acc/2., acc/2.);
  gsl_odeiv_evolve          *e    = gsl_odeiv_evolve_alloc (2);
  double                    *y    = new double [2]; 
  double                     r, h, psir, psiz, rt, zt;
  
  for (int j = 1; j < J; j++)
    {
      r    = 0.;
      h    = h0;
      y[0] = RP[j];
      y[1] = Z[jc];
      qgp  = QGP[j] * fabs(psic);
	      
      gsl_matrix_set (Rst, j, 0, (y[0]-Raxis) /rP[j]);
      gsl_matrix_set (Zst, j, 0, y[1] /rP[j]);
      psir = GetPsiR (y[0], y[1]);
      psiz = GetPsiZ (y[0], y[1]);
      rt   = - qgp * y[0] * psiz /rP[j];
      zt   = + qgp * y[0] * psir /rP[j];
      gsl_matrix_set (Rt, j, 0, rt);
      gsl_matrix_set (Zt, j, 0, zt);

      for (int k = 1; k < K; k++)
	{
	  while (r < th[k])
	    {
	      int status = gsl_odeiv_evolve_apply (e, c, s, &sys3, &r, th[k], &h, y);
	      
	      if (status != GSL_SUCCESS)
		{
		  printf ("Profile: status != GSL_SUCCESS\n");
		  exit (1);
		}
	    }
	  gsl_matrix_set (Rst, j, k, (y[0]-Raxis) /rP[j]);
	  gsl_matrix_set (Zst, j, k, y[1] /rP[j]);

	  psir = GetPsiR (y[0], y[1]);
	  psiz = GetPsiZ (y[0], y[1]);
	  rt   = - qgp * y[0] * psiz /rP[j];
	  zt   = + qgp * y[0] * psir /rP[j];
	  gsl_matrix_set (Rt, j, k, rt);
	  gsl_matrix_set (Zt, j, k, zt);
	}
     }

  gsl_odeiv_evolve_free  (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free    (s);
  delete[]                y;
}

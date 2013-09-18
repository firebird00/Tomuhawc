// InterpolatePsi.cpp

#include "Flux.h"

// ##############################################
// Function to interpolate Psi on uniform 2D grid
// order = 0: Psi
// order = 1: Psi_x
// order = 2: Psi_y
// order = 3: Psi_xx
// order = 4: Psi_yy
// ##############################################
double Flux::InterpolatePsi (double RR, double ZZ, int order)
{
  double i  = double (II-1) * (RR - R[0]) /(R[II-1] - R[0]);
  int    i0 = int (i); 
  if (i0 < 0 || i0 > II-1)
    {
      printf ("Interpolation error: II = %3d i0 = %3d r = %11.4e R[0] = %11.4e R[II-1] = %11.4e\n", 
	      II, i0, RR, R[0], R[II-1]);
      exit (1);
    }
  if (i - double (i0) > 0.5)
    i0 += 1;
  if (i0 == 0)
    i0 += 1;
  if (i0 == II-1)
    i0 -= 1;

  double j  = double (JJ-1) * (ZZ - Z[0]) /(Z[JJ-1] - Z[0]);
  int    j0 = int (j); 
  if (j0 < 0 || j0 > JJ-1)
    {
      printf ("Interpolation error: JJ = %3d j0 = %3d z = %11.4e Z[0] = %11.4e Z[JJ-1] = %11.4e\n", 
	      JJ, j0, ZZ, Z[0], Z[JJ-1]);
      exit (1);
    }
  if (j - double (j0) > 0.5)
    j0 += 1;
  if (j0 == 0)
    j0 += 1;
  if (j0 == JJ-1)
    j0 -= 1;

  double dR = (R[1] - R[0]);
  double dZ = (Z[1] - Z[0]);
  double x  = (RR - R[i0]) /dR;
  double z  = (ZZ - Z[j0]) /dZ;

  double val;
  if (order == 0)
    {
      double xm = + 0.5 * x * (x-1.);
      double x0 = - (x+1.) * (x-1.);
      double x1 = + 0.5 * x * (x+1.);
      double zm = + 0.5 * z * (z-1.);
      double z0 = - (z+1.) * (z-1.);
      double z1 = + 0.5 * z * (z+1.);

      val =
	+ xm * (+ zm * gsl_matrix_get (Psi, i0-1, j0-1)
		+ z0 * gsl_matrix_get (Psi, i0-1, j0  )
		+ z1 * gsl_matrix_get (Psi, i0-1, j0+1))
	+ x0 * (+ zm * gsl_matrix_get (Psi, i0  , j0-1)
		+ z0 * gsl_matrix_get (Psi, i0  , j0  )
		+ z1 * gsl_matrix_get (Psi, i0  , j0+1))
	+ x1 * (+ zm * gsl_matrix_get (Psi, i0+1, j0-1)
		+ z0 * gsl_matrix_get (Psi, i0+1, j0  )
		+ z1 * gsl_matrix_get (Psi, i0+1, j0+1));
    }
  else if (order == 1)
    {
      double xm = + x - 0.5;
      double x0 = - 2. * x;
      double x1 = + x + 0.5;
      double zm = + 0.5 * z * (z-1.);
      double z0 = - (z+1.) * (z-1.);
      double z1 = + 0.5 * z * (z+1.);
 
      val = 
	+ xm * (+ zm * gsl_matrix_get (Psi, i0-1, j0-1)
		+ z0 * gsl_matrix_get (Psi, i0-1, j0  )
		+ z1 * gsl_matrix_get (Psi, i0-1, j0+1))
	+ x0 * (+ zm * gsl_matrix_get (Psi, i0  , j0-1)
		+ z0 * gsl_matrix_get (Psi, i0  , j0  )
		+ z1 * gsl_matrix_get (Psi, i0  , j0+1))
	+ x1 * (+ zm * gsl_matrix_get (Psi, i0+1, j0-1)
		+ z0 * gsl_matrix_get (Psi, i0+1, j0  )
		+ z1 * gsl_matrix_get (Psi, i0+1, j0+1));

      val /= dR;
    }
  else if (order == 2)
    {
      double xm = + 0.5 * x * (x-1.);
      double x0 = - (x+1.) * (x-1.);
      double x1 = + 0.5 * x * (x+1.);
      double zm = + z - 0.5;
      double z0 = - 2. * z;
      double z1 = + z + 0.5;

      val =
	+ xm * (+ zm * gsl_matrix_get (Psi, i0-1, j0-1)
		+ z0 * gsl_matrix_get (Psi, i0-1, j0  )
		+ z1 * gsl_matrix_get (Psi, i0-1, j0+1))
	+ x0 * (+ zm * gsl_matrix_get (Psi, i0  , j0-1)
		+ z0 * gsl_matrix_get (Psi, i0  , j0  )
		+ z1 * gsl_matrix_get (Psi, i0  , j0+1))
	+ x1 * (+ zm * gsl_matrix_get (Psi, i0+1, j0-1)
		+ z0 * gsl_matrix_get (Psi, i0+1, j0  )
		+ z1 * gsl_matrix_get (Psi, i0+1, j0+1));

      val /= dZ;
    }
  else if (order == 3)
    {
      double xm = + 1.;
      double x0 = - 2.;
      double x1 = + 1.;
      double zm = + 0.5 * z * (z-1.);
      double z0 = - (z+1.) * (z-1.);
      double z1 = + 0.5 * z * (z+1.);
 
      val = 
	+ xm * (+ zm * gsl_matrix_get (Psi, i0-1, j0-1)
		+ z0 * gsl_matrix_get (Psi, i0-1, j0  )
		+ z1 * gsl_matrix_get (Psi, i0-1, j0+1))
	+ x0 * (+ zm * gsl_matrix_get (Psi, i0  , j0-1)
		+ z0 * gsl_matrix_get (Psi, i0  , j0  )
		+ z1 * gsl_matrix_get (Psi, i0  , j0+1))
	+ x1 * (+ zm * gsl_matrix_get (Psi, i0+1, j0-1)
		+ z0 * gsl_matrix_get (Psi, i0+1, j0  )
		+ z1 * gsl_matrix_get (Psi, i0+1, j0+1));
      
      val /= (dR*dR);
    }
  else if (order == 4)
    {
      double xm = + 0.5 * x * (x-1.);
      double x0 = - (x+1.) * (x-1.);
      double x1 = + 0.5 * x * (x+1.);
      double zm = + 1.;
      double z0 = - 2.;
      double z1 = + 1.;

      val =
	+ xm * (+ zm * gsl_matrix_get (Psi, i0-1, j0-1)
		+ z0 * gsl_matrix_get (Psi, i0-1, j0  )
		+ z1 * gsl_matrix_get (Psi, i0-1, j0+1))
	+ x0 * (+ zm * gsl_matrix_get (Psi, i0  , j0-1)
		+ z0 * gsl_matrix_get (Psi, i0  , j0  )
		+ z1 * gsl_matrix_get (Psi, i0  , j0+1))
	+ x1 * (+ zm * gsl_matrix_get (Psi, i0+1, j0-1)
		+ z0 * gsl_matrix_get (Psi, i0+1, j0  )
		+ z1 * gsl_matrix_get (Psi, i0+1, j0+1));

      val /= (dZ*dZ);
    }

  return val;
}

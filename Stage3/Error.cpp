// Program to calculate optimum wall displacement to exert torque on given rational surface

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>

double Extrapolate (int I, double *X, gsl_matrix *Y, double x, int j);
double Interpolate (int I, double *X, gsl_matrix *Y, double x, int j);

int main ()
{
  // Determine J and K
  int J, K;
  FILE *file = fopen ("../Stage2/jk.out", "r");
  if (file == NULL)
    {
      printf ("Error opening ../Stage2/jk.out\n");
      exit (1);
    }	
  if (fscanf (file, "%d %d", &J, &K) != 2)
    {
      printf ("Error reading ../Stage2/jk.out\n");
      exit (1);
    }
  fclose (file);

  // Read r and t arrays
  double *r = new double[J];
  double *t = new double[K];
  file = fopen ("../Stage2/rx.out", "r");
  if (file == NULL)
    {
      printf ("Error opening ../Stage2/rx.out\n");
      exit (1);
    }	
  for (int j = 0; j < J; j++)
    if (fscanf (file, "%lf", &r[j]) != 1)
      {
	printf ("Error reading ../Stage2/rx.out\n");
	exit (1);
      }	
  fclose (file);
  file = fopen ("../Stage2/tt.out", "r");
  if (file == NULL)
    {
      printf ("Error opening ../Stage2/tt.out\n");
      exit (1);
    }	
  for (int k = 0; k < K; k++)
    if (fscanf (file, "%lf", &t[k]) != 1)
      {
	printf ("Error reading ../Stage2/tt.out\n");
	exit (1);
      }	
  fclose (file);

  // Read R and Z arrays
  gsl_matrix *R = gsl_matrix_alloc (J, K);
  gsl_matrix *Z = gsl_matrix_alloc (J, K);
  file = fopen ("../Stage2/Rt.out", "r");
  if (file == NULL)
    {
      printf ("Error opening ../Stage2/Rt.out\n");
      exit (1);
    }
  for (int j = 0; j < J; j++)
    for (int k = 0; k < K; k++)
      {
	double val;
	if (fscanf (file, "%lf", &val) != 1)
	  {
	    printf ("Error reading ../Stage2/Rt.out\n");
	    exit (1);
	  }
	gsl_matrix_set (R, j, k, val);
      }
  fclose (file);
  file = fopen ("../Stage2/Zt.out", "r");
  if (file == NULL)
    {
      printf ("Error opening ../Stage2/Zt.out\n");
      exit (1);
    }
  for (int j = 0; j < J; j++)
    for (int k = 0; k < K; k++)
      {
	double val;
	if (fscanf (file, "%lf", &val) != 1)
	  {
	    printf ("Error reading ../Stage2/Zt.out\n");
	    exit (1);
	  }
	gsl_matrix_set (Z, j, k, val);
      }
  fclose (file);

  // Read in dim and num
  int dim, num, vac;
  file = fopen ("info.out", "r");
  if (file == NULL)
    {
      printf ("Error opening info.out\n");
      exit (1);
    }
  if (fscanf (file, "%d %d %d\n", &dim, &num, &vac) != 3)
    {
      printf ("Error reading info.out\n");
      exit (1);
    }
  fclose (file);
  int dim1 = dim*num;

  // Determine parity
  int tear;
  printf ("\nTearing/twisting parity (1/0) ? ");
  scanf ("%d", &tear);

  // Determine rational surface
  int rat;
  printf ("Rational surface number (1,2,3, etc) ? ");
  scanf ("%d", &rat);
  if (rat < 1 || rat > num)
    {
      printf ("Error: Invalid rational surface number\n");
      return 1;
    }
  rat = rat - 1;

  // Read in displacement 
  double *xim = new double [dim1];
  if (tear)
    {
      file = fopen ("xie.out", "r");
      if (file == NULL)
	{
	  printf ("Error opening xie.out\n");
	  exit (1);
	}
    }
  else
    {	
      file = fopen ("xio.out", "r");  
      if (file == NULL)
	{
	  printf ("Error opening xio.out\n");
	  exit (1);
	} 
    }	 
  for (int i = 0; i < dim1; i++)
    if (fscanf (file, "%lf", &xim[i]) != 1)
      {
	printf ("Error reading displacement.\n");
	exit (1);
      }
  
  fclose (file);

  // Read in mode numbers
  int *m = new int [dim1];
  file = fopen ("modes.out", "r");
  if (file == NULL)
    {
      printf ("Error opening modes.out\n");
      exit (1);
    } 
  for (int i = 0; i < dim1; i++)
    if (fscanf (file, "%d", &m[i]) != 1)
      {
	printf ("Error reading modes.out\n");
	exit (1);
      } 
  fclose (file);

  // Calculate normalization constant
  double amp;
  printf ("Amplitude ? ");
  scanf ("%lf", &amp);
  printf ("\n");

  double norm = 0.;
  for (int i = 0; i < dim; i++)
    norm += xim[rat*dim+i]*xim[rat*dim+i];
  norm = sqrt(norm)/amp;

  // Read rational surface data
  double *rs = new double[vac];
  double dum;
  file = fopen ("ratsur.out", "r");
  if (file == NULL)
    {
      printf ("Error opening ratsur.out\n");
      exit (1);
    } 
  for (int i = 0; i < vac; i++)
    if (fscanf (file, "%lf %lf", &rs[i], &dum) != 2)
      {
	printf ("Error reading ratsur.out\n");
	exit (1);
      }	
  fclose (file);

  // Plot boundary displacement
  file = fopen ("error.out", "w");
  for (int k = 0; k < K; k++)
    {
      double rb = 1.;
      double Rb = Extrapolate (J, r, R, rb,      k);
      double Zb = Extrapolate (J, r, Z, rb,      k);
      double Rx = Interpolate (J, r, R, rs[rat], k);
      double Zx = Interpolate (J, r, Z, rs[rat], k);

      double th = t[k];

      double xic = 0., xis = 0.;
      for (int i = 0; i < dim; i++)
	{
	  xic += xim[rat*dim+i] * cos(double(m[rat*dim+i]) * th);
	  xis += xim[rat*dim+i] * sin(double(m[rat*dim+i]) * th);
	}
      xic /= norm;
      xis /= norm;
      double xi = sqrt(xic*xic + xis*xis);

      double Rp  = Extrapolate (J, r, R, 1.+xi,  k);
      double Zp  = Extrapolate (J, r, Z, 1.+xi,  k);
      double Rm  = Extrapolate (J, r, R, 1.-xi,  k);
      double Zm  = Extrapolate (J, r, Z, 1.-xi,  k);
      double Rcp = Extrapolate (J, r, R, 1.+xic, k);
      double Zcp = Extrapolate (J, r, Z, 1.+xic, k);
      double Rcm = Extrapolate (J, r, R, 1.-xic, k);
      double Zcm = Extrapolate (J, r, Z, 1.-xic, k);
      double Rsp = Extrapolate (J, r, R, 1.+xis, k);
      double Zsp = Extrapolate (J, r, Z, 1.+xis, k);
      double Rsm = Extrapolate (J, r, R, 1.-xis, k);
      double Zsm = Extrapolate (J, r, Z, 1.-xis, k);

      fprintf (file, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", 
	       Rb, Zb, Rx, Zx, Rp, Zp, Rm, Zm, Rcp, Zcp, Rcm, Zcm, Rsp, Zsp, Rsm, Zsm);
    }
  for (int k = 1; k < K; k++)
    {
      double rb = 1.;
      double Rb = Extrapolate (J, r, R, rb,      K-1-k);
      double Zb = Extrapolate (J, r, Z, rb,      K-1-k);
      double Rx = Interpolate (J, r, R, rs[rat], K-1-k);
      double Zx = Interpolate (J, r, Z, rs[rat], K-1-k);

      double th = t[K-1-k];

      double xic = 0., xis = 0.;
      for (int i = 0; i < dim; i++)
	{
	  xic += xim[rat*dim+i] * cos(double(m[rat*dim+i]) * th);
	  xis += xim[rat*dim+i] * sin(double(m[rat*dim+i]) * th);
	}
      xic /= norm;
      xis /= norm;
      double xi = sqrt(xic*xic + xis*xis);

      double Rp  = Extrapolate (J, r, R, 1.+xi,  K-1-k);
      double Zp  = Extrapolate (J, r, Z, 1.+xi,  K-1-k);
      double Rm  = Extrapolate (J, r, R, 1.-xi,  K-1-k);
      double Zm  = Extrapolate (J, r, Z, 1.-xi,  K-1-k);
      double Rcp = Extrapolate (J, r, R, 1.+xic, K-1-k);
      double Zcp = Extrapolate (J, r, Z, 1.+xic, K-1-k);
      double Rcm = Extrapolate (J, r, R, 1.-xic, K-1-k);
      double Zcm = Extrapolate (J, r, Z, 1.-xic, K-1-k);
      double Rsp = Extrapolate (J, r, R, 1.-xis, K-1-k);
      double Zsp = Extrapolate (J, r, Z, 1.-xis, K-1-k);
      double Rsm = Extrapolate (J, r, R, 1.+xis, K-1-k);
      double Zsm = Extrapolate (J, r, Z, 1.+xis, K-1-k);


      fprintf (file, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", 
	       Rb, -Zb, Rx, -Zx, Rp, -Zp, Rm, -Zm, Rcp, -Zcp, Rcm, -Zcm, Rsp, -Zsp, Rsm, -Zsm);
    }
  fclose (file);

  // Clean up
  delete[] r; delete[] t;
  gsl_matrix_free (R); gsl_matrix_free (Z);
  delete[] xim; delete[] m; delete[] rs;
}

double Extrapolate (int I, double *X, gsl_matrix *Y, double x, int j)
{
  int    i0 = I-2;
  double sm = (x-X[i0  ]) * (x-X[i0+1]) /(X[i0-1]-X[i0  ]) /(X[i0-1]-X[i0+1]);
  double s0 = (x-X[i0-1]) * (x-X[i0+1]) /(X[i0  ]-X[i0-1]) /(X[i0  ]-X[i0+1]);
  double s1 = (x-X[i0-1]) * (x-X[i0  ]) /(X[i0+1]-X[i0-1]) /(X[i0+1]-X[i0  ]);
      
  return sm * gsl_matrix_get (Y, i0-1, j) + s0 * gsl_matrix_get (Y, i0, j) + s1 * gsl_matrix_get (Y, i0+1, j);
}

double Interpolate (int I, double *X, gsl_matrix *Y, double x, int j)
{  
  int i0 = 0;
  for (int i = 1; i < I; i++)
    if (x > X[i])
      i0 = i;
  if (i0 < 0 || i0 > I-1)
    {
      printf ("Interpolation error: I = %3d i0 = %3d x = %11.4e X[0] = %11.4e X[I-1] = %11.4e\n",
	      I, i0, x, X[0], X[I-1]);
      exit (1);
    }
  if (x - X[i0] > 0.5 *(X[i0+1] - X[i0]))
    i0 += 1;
  if (i0 == 0)
    i0 += 1;
  if (i0 == I-1)
    i0 -= 1;

  double val;
  double sm = (x-X[i0  ]) * (x-X[i0+1]) /(X[i0-1]-X[i0  ]) /(X[i0-1]-X[i0+1]);
  double s0 = (x-X[i0-1]) * (x-X[i0+1]) /(X[i0  ]-X[i0-1]) /(X[i0  ]-X[i0+1]);
  double s1 = (x-X[i0-1]) * (x-X[i0  ]) /(X[i0+1]-X[i0-1]) /(X[i0+1]-X[i0  ]);
  
  val = sm * gsl_matrix_get (Y, i0-1, j) + s0 * gsl_matrix_get (Y, i0, j) 
    + s1 * gsl_matrix_get (Y, i0+1, j);
  
  return val;
}

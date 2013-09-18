// Equilibrium.cpp

#include "Tomuhawc.h"

// ###########################################
// Function to read in equilibrium metric data
//
// sw = 0: read in metric data 
// sw = 1: clean up 
//
// interactive .. set if in interactive mode
//
// ###########################################
void Thawc::Equilibrium (int sw, int interactive)
{
  if (sw == 0)
    {
      // Read in grid data
      FILE *file = OpenFiler ("Stage2/grid.out");
      int    err = fscanf (file, "%d %d", &nint, &diag);
      if (err != 2)
	{
	  printf ("Error: Cannot read grid.out\n");
	  exit (1);
	}
      fclose (file);
      printf ("nint = %4d  diag = %4d\n",
	      nint, diag);
      if (nint < 2 || diag < 0)
	{
	  printf ("Error: Bad data values for nint and/or diag\n");
	  exit (1);
	}
      if (side > diag)
	{
	  printf ("Error: side larger than diag\n");
	  exit (1);
	}

      // Read in Gamma, abar, a, rb
      file = OpenFiler ("Stage2/info.out");
      err  = fscanf (file, "%lf %lf %lf %lf", &Gamma, &abar, &a, &rb);
      if (err != 4)
	{
	  printf ("Error: Cannot read info.out\n");
	  exit (1);
	}
      fclose (file);
      printf ("abar = %11.4e  a   = %11.4e\n", abar, a);

      // Read in profile data
      double *r   = new double [nint];
      double *q   = new double [nint];
      double *qp  = new double [nint];
      double *p1  = new double [nint];
      double *p2  = new double [nint];
      double *p3  = new double [nint];
      double *p4  = new double [nint];
      double *p5  = new double [nint];
      double *p6  = new double [nint];
      double *p1p = new double [nint];
      double *p2p = new double [nint];
      double *p3p = new double [nint];

      file = OpenFiler ("Stage2/profile.out");
      for (int i = 0; i < nint; i++)
	{
	  err = fscanf (file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
			&r[i], &q[i], &qp[i], &p1[i], &p2[i], &p3[i], 
			&p4[i], &p5[i], &p6[i], &p1p[i], &p2p[i], &p3p[i]);
	  if (err != 12)
	    {
	      printf ("Error: Cannot read profile.out\n");
	      exit (1);
	    }
	}
      fclose (file);
      ra = r[nint-1]; 
      qa = q[nint-1];
      printf ("r0/a = %11.4e  b/a = %11.4e  a/a = %11.4e\n", r0, rb, ra);

      // Fit profile data
      Aq      = gsl_interp_accel_alloc ();
      Aqp     = gsl_interp_accel_alloc ();
      Aprof1  = gsl_interp_accel_alloc ();
      Aprof2  = gsl_interp_accel_alloc ();
      Aprof3  = gsl_interp_accel_alloc ();
      Aprof4  = gsl_interp_accel_alloc ();
      Aprof5  = gsl_interp_accel_alloc ();
      Aprof6  = gsl_interp_accel_alloc ();
      Aprof1p = gsl_interp_accel_alloc ();
      Aprof2p = gsl_interp_accel_alloc ();
      Aprof3p = gsl_interp_accel_alloc ();

      Sq      = gsl_spline_alloc (gsl_interp_cspline, nint);
      Sqp     = gsl_spline_alloc (gsl_interp_cspline, nint);
      Sprof1  = gsl_spline_alloc (gsl_interp_cspline, nint);
      Sprof2  = gsl_spline_alloc (gsl_interp_cspline, nint);
      Sprof3  = gsl_spline_alloc (gsl_interp_cspline, nint);
      Sprof4  = gsl_spline_alloc (gsl_interp_cspline, nint);
      Sprof5  = gsl_spline_alloc (gsl_interp_cspline, nint);
      Sprof6  = gsl_spline_alloc (gsl_interp_cspline, nint);
      Sprof1p = gsl_spline_alloc (gsl_interp_cspline, nint);
      Sprof2p = gsl_spline_alloc (gsl_interp_cspline, nint);
      Sprof3p = gsl_spline_alloc (gsl_interp_cspline, nint);

      gsl_spline_init (Sq,      r, q,   nint);
      gsl_spline_init (Sqp,     r, qp,  nint);
      gsl_spline_init (Sprof1,  r, p1,  nint);
      gsl_spline_init (Sprof2,  r, p2,  nint);
      gsl_spline_init (Sprof3,  r, p3,  nint);
      gsl_spline_init (Sprof4,  r, p4,  nint);
      gsl_spline_init (Sprof5,  r, p5,  nint);
      gsl_spline_init (Sprof6,  r, p6,  nint);
      gsl_spline_init (Sprof1p, r, p1p, nint);
      gsl_spline_init (Sprof2p, r, p2p, nint);
      gsl_spline_init (Sprof3p, r, p3p, nint);

      // Read in metric functions
      double **M1  = new double* [diag+1];
      double **M2  = new double* [diag+1];
      double **M3  = new double* [diag+1];
      double **M4  = new double* [diag+1];
      double **M5  = new double* [diag+1];
      double **M6  = new double* [diag+1];
      double **M7  = new double* [diag+1];
      double **M1P = new double* [diag+1];

      for (int j = 0; j <= diag; j++) 
	{	
	  M1 [j] = new double [nint];
	  M2 [j] = new double [nint];
	  M3 [j] = new double [nint];
	  M4 [j] = new double [nint];
	  M5 [j] = new double [nint];
	  M6 [j] = new double [nint];
	  M7 [j] = new double [nint];
	  M1P[j] = new double [nint];
	}
      double *M8  = new double [nint];
      double *M9  = new double [nint];
      double *M3P = new double [nint];   
      double *M4P = new double [nint];
      
      file = OpenFiler ("Stage2/M1.out");
      for (int i = 0; i < nint; i++)
	for (int j = 0; j <= diag; j++)
	  {
	    err = fscanf (file, "%lf", &M1[j][i]);
		    
	    if (err == EOF)
	      {
		printf ("Error: Cannot read M1.out\n");
		exit (1);
	      }
	  }	
      fclose (file);
      file = OpenFiler ("Stage2/M2.out");
      for (int i = 0; i < nint; i++)
	for (int j = 0; j <= diag; j++)
	  {
	    err = fscanf (file, "%lf", &M2[j][i]);
		    
	    if (err == EOF)
	      {
		printf ("Error: Cannot read M2.out\n");
		exit (1);
	      }
	  }	
      fclose (file);
      file = OpenFiler ("Stage2/M3.out");
      for (int i = 0; i < nint; i++)
	for (int j = 0; j <= diag; j++)
	  {
	    err = fscanf (file, "%lf", &M3[j][i]);
		    
	    if (err == EOF)
	      {
		printf ("Error: Cannot read M3.out\n");
		exit (1);
	      }
	  }	
      fclose (file);
      file = OpenFiler ("Stage2/M4.out");
      for (int i = 0; i < nint; i++)
	for (int j = 0; j <= diag; j++)
	  {
	    err = fscanf (file, "%lf", &M4[j][i]);
		    
	    if (err == EOF)
	      {
		printf ("Error: Cannot read M4.out\n");
		exit (1);
	      }
	  }	
      fclose (file);
      file = OpenFiler ("Stage2/M5.out");
      for (int i = 0; i < nint; i++)
	for (int j = 0; j <= diag; j++)
	  {
	    err = fscanf (file, "%lf", &M5[j][i]);
		    
	    if (err == EOF)
	      {
		printf ("Error: Cannot read M5.out\n");
		exit (1);
	      }
	  }	
      fclose (file);
      file = OpenFiler ("Stage2/M6.out");
      for (int i = 0; i < nint; i++)
	for (int j = 0; j <= diag; j++)
	  {
	    err = fscanf (file, "%lf", &M6[j][i]);
		    
	    if (err == EOF)
	      {
		printf ("Error: Cannot read M6.out\n");
		exit (1);
	      }
	  }	
      fclose (file);
      file = OpenFiler ("Stage2/M7.out");
      for (int i = 0; i < nint; i++)
	for (int j = 0; j <= diag; j++)
	  {
	    err = fscanf (file, "%lf", &M7[j][i]);
		    
	    if (err == EOF)
	      {
		printf ("Error: Cannot read M7.out\n");
		exit (1);
	      }
	  }	
      fclose (file);
      file = OpenFiler ("Stage2/M8.out");
      for (int i = 0; i < nint; i++)
	{
	  err = fscanf (file, "%lf", &M8[i]);
	  
	  if (err == EOF)
	    {
	      printf ("Error: Cannot read M8.out\n");
	      exit (1);
	    }
	}	
      fclose (file);
      file = OpenFiler ("Stage2/M9.out");
       for (int i = 0; i < nint; i++)
	{
	  err = fscanf (file, "%lf", &M9[i]);
	  
	  if (err == EOF)
	    {
	      printf ("Error: Cannot read M9.out\n");
	      exit (1);
	    }
	}	
      fclose (file);
      file = OpenFiler ("Stage2/M1P.out");
      for (int i = 0; i < nint; i++)
	for (int j = 0; j <= diag; j++)
	  {
	    err = fscanf (file, "%lf", &M1P[j][i]);
		    
	    if (err == EOF)
	      {
		printf ("Error: Cannot read M1P.out\n");
		exit (1);
	      }
	  }	
      fclose (file);
          file = OpenFiler ("Stage2/M3P.out");
      for (int i = 0; i < nint; i++)
	{
	  err = fscanf (file, "%lf", &M3P[i]);
	  
	  if (err == EOF)
	    {
	      printf ("Error: Cannot read M3P.out\n");
	      exit (1);
	    }
	}	
      fclose (file);
      file = OpenFiler ("Stage2/M4P.out");
       for (int i = 0; i < nint; i++)
	{
	  err = fscanf (file, "%lf", &M4P[i]);
	  
	  if (err == EOF)
	    {
	      printf ("Error: Cannot read M4P.out\n");
	      exit (1);
	    }
	}	
      fclose (file);
      
      // Fit metric functions
      AM1  = new gsl_interp_accel* [diag+1];
      AM2  = new gsl_interp_accel* [diag+1];
      AM3  = new gsl_interp_accel* [diag+1];
      AM4  = new gsl_interp_accel* [diag+1];
      AM5  = new gsl_interp_accel* [diag+1];
      AM6  = new gsl_interp_accel* [diag+1];
      AM7  = new gsl_interp_accel* [diag+1];
      AM1P = new gsl_interp_accel* [diag+1];

      SM1  = new gsl_spline* [diag+1];
      SM2  = new gsl_spline* [diag+1];
      SM3  = new gsl_spline* [diag+1];
      SM4  = new gsl_spline* [diag+1];
      SM5  = new gsl_spline* [diag+1];
      SM6  = new gsl_spline* [diag+1];
      SM7  = new gsl_spline* [diag+1];
      SM1P = new gsl_spline* [diag+1];

      for (int j = 0; j <= diag; j++) 
	{	
	  AM1 [j] = gsl_interp_accel_alloc ();
	  AM2 [j] = gsl_interp_accel_alloc ();
	  AM3 [j] = gsl_interp_accel_alloc ();
	  AM4 [j] = gsl_interp_accel_alloc ();
	  AM5 [j] = gsl_interp_accel_alloc ();
	  AM6 [j] = gsl_interp_accel_alloc ();
	  AM7 [j] = gsl_interp_accel_alloc ();
	  AM1P[j] = gsl_interp_accel_alloc ();

	  SM1 [j] = gsl_spline_alloc (gsl_interp_cspline, nint);
	  SM2 [j] = gsl_spline_alloc (gsl_interp_cspline, nint);
	  SM3 [j] = gsl_spline_alloc (gsl_interp_cspline, nint);
	  SM4 [j] = gsl_spline_alloc (gsl_interp_cspline, nint);
	  SM5 [j] = gsl_spline_alloc (gsl_interp_cspline, nint);
	  SM6 [j] = gsl_spline_alloc (gsl_interp_cspline, nint);
	  SM7 [j] = gsl_spline_alloc (gsl_interp_cspline, nint);
	  SM1P[j] = gsl_spline_alloc (gsl_interp_cspline, nint);
	}
      AM8  = gsl_interp_accel_alloc ();
      AM9  = gsl_interp_accel_alloc ();
      AM3P = gsl_interp_accel_alloc ();
      AM4P = gsl_interp_accel_alloc ();

      SM8  = gsl_spline_alloc (gsl_interp_cspline, nint);
      SM9  = gsl_spline_alloc (gsl_interp_cspline, nint);
      SM3P = gsl_spline_alloc (gsl_interp_cspline, nint);
      SM4P = gsl_spline_alloc (gsl_interp_cspline, nint);

      for (int j = 0; j <= diag; j++)
	{
	  gsl_spline_init (SM1 [j], r, M1 [j], nint);
	  gsl_spline_init (SM2 [j], r, M2 [j], nint);
	  gsl_spline_init (SM3 [j], r, M3 [j], nint);
	  gsl_spline_init (SM4 [j], r, M4 [j], nint);
	  gsl_spline_init (SM5 [j], r, M5 [j], nint);
	  gsl_spline_init (SM6 [j], r, M6 [j], nint);
	  gsl_spline_init (SM7 [j], r, M7 [j], nint);
	  gsl_spline_init (SM1P[j], r, M1P[j], nint);
	}   
      gsl_spline_init (SM8,  r, M8,  nint);
      gsl_spline_init (SM9,  r, M9,  nint);
      gsl_spline_init (SM3P, r, M3P, nint);
      gsl_spline_init (SM4P, r, M4P, nint);

      // Read in vacuum data
      Ptor0 = gsl_matrix_alloc (diag+1, diag+1);
      Qtor0 = gsl_matrix_alloc (diag+1, diag+1);
      Ptor1 = gsl_matrix_alloc (diag+1, diag+1);
      Qtor1 = gsl_matrix_alloc (diag+1, diag+1);
      Pm    = new double [diag+1];
      Qm    = new double [diag+1];
    
      double val;
      file = OpenFiler ("Stage2/P.out");
      for (int i = 0; i <= diag; i++)
	for (int j = 0; j <= diag; j++)
	  {
	    err = fscanf (file, "%lf", &val);
		    
	    if (err == EOF)
	      {
		printf ("Error: Cannot read P.out\n");
		exit (1);
	      }

	    gsl_matrix_set (Ptor0, i, j, val);
	  }	
      fclose (file);

      file = OpenFiler ("Stage2/Q.out");
      for (int i = 0; i <= diag; i++)
	for (int j = 0; j <= diag; j++)
	  {
	    err = fscanf (file, "%lf", &val);
		    
	    if (err == EOF)
	      {
		printf ("Error: Cannot read Q.out\n");
		exit (1);
	      }

	    gsl_matrix_set (Qtor0, i, j, val);
	  }	
      fclose (file);

      file = OpenFiler ("Stage2/PP.out");
      for (int i = 0; i <= diag; i++)
	for (int j = 0; j <= diag; j++)
	  {
	    err = fscanf (file, "%lf", &val);
		    
	    if (err == EOF)
	      {
		printf ("Error: Cannot read PP.out\n");
		exit (1);
	      }

	    gsl_matrix_set (Ptor1, i, j, val);
	  }	
      fclose (file);

      file = OpenFiler ("Stage2/QQ.out");
      for (int i = 0; i <= diag; i++)
	for (int j = 0; j <= diag; j++)
	  {
	    err = fscanf (file, "%lf", &val);
		    
	    if (err == EOF)
	      {
		printf ("Error: Cannot read QQ.out\n");
		exit (1);
	      }

	    gsl_matrix_set (Qtor1, i, j, val);
	  }	
      fclose (file);

      file = OpenFiler ("Stage2/PQm.out");
      for (int i = 0; i <= diag; i++)
	if (fscanf (file, "%lf %lf", &Pm[i], &Qm[i]) != 2)
	  {
	    printf ("Error: Cannot read PQm.out\n");
	    exit (1);
	  }
      fclose (file);
        
      // +++++++++++++++++++++++++++++++++++++
      // Calculate GGJ layer theory parameters
      // +++++++++++++++++++++++++++++++++++++
      if (interactive) file = OpenFile ("Stage3/GGJ.out");
      for (int i = 1; i < nint; i++)
	{
	  double EF, HH, ET, FT, HT, FA, FR;
	  if (r[i] < 1.)
	    GGJCalc (r[i], EF, HH, ET, FT, HT, FA, FR);
	  if (interactive) 
	    fprintf (file, "%e %e %e %e %e %e %e %e\n", r[i], EF, HH, ET, FT, HT, FA, FR);
	}
      if (interactive) fclose (file);

      // Clean up
      delete[] r;   delete[] q;   delete[] qp; 
      delete[] p1;  delete[] p2;  delete[] p3; 
      delete[] p4;  delete[] p5;  delete[] p6;
      delete[] p1p; delete[] p2p; delete[] p3p;
      for (int j = 0; j <= diag; j++) 
	{	
	  delete[] M1 [j];
	  delete[] M2 [j];
	  delete[] M3 [j];
	  delete[] M4 [j];
	  delete[] M5 [j];
	  delete[] M6 [j];
	  delete[] M7 [j];
	  delete[] M1P[j];
	}
      delete[] M1;  delete[] M2;  delete[] M3; 
      delete[] M4;  delete[] M5;  delete[] M6; 
      delete[] M7;  delete[] M8;  delete[] M9;
      delete[] M1P; delete[] M3P; delete[] M4P;
    }
  else 
    {
      // Deallocate memory
      gsl_interp_accel_free (Aq);
      gsl_interp_accel_free (Aqp);
      gsl_interp_accel_free (Aprof1);
      gsl_interp_accel_free (Aprof2);
      gsl_interp_accel_free (Aprof3);
      gsl_interp_accel_free (Aprof4);     
      gsl_interp_accel_free (Aprof5);
      gsl_interp_accel_free (Aprof6);
      gsl_interp_accel_free (Aprof1p);
      gsl_interp_accel_free (Aprof2p);
      gsl_interp_accel_free (Aprof3p);

      gsl_spline_free (Sq);
      gsl_spline_free (Sqp);
      gsl_spline_free (Sprof1);
      gsl_spline_free (Sprof2);
      gsl_spline_free (Sprof3);
      gsl_spline_free (Sprof4);
      gsl_spline_free (Sprof5);
      gsl_spline_free (Sprof6);
      gsl_spline_free (Sprof1p);
      gsl_spline_free (Sprof2p);
      gsl_spline_free (Sprof3p);

      for (int j = 0; j <= diag; j++) 
	{	
	  gsl_interp_accel_free (AM1 [j]);
	  gsl_interp_accel_free (AM2 [j]);
	  gsl_interp_accel_free (AM3 [j]);
	  gsl_interp_accel_free (AM4 [j]);
	  gsl_interp_accel_free (AM5 [j]);
	  gsl_interp_accel_free (AM6 [j]);
	  gsl_interp_accel_free (AM7 [j]);
	  gsl_interp_accel_free (AM1P[j]);
	  
	  gsl_spline_free (SM1 [j]);
	  gsl_spline_free (SM2 [j]);
	  gsl_spline_free (SM3 [j]);
	  gsl_spline_free (SM4 [j]);
	  gsl_spline_free (SM5 [j]);
	  gsl_spline_free (SM6 [j]);
	  gsl_spline_free (SM7 [j]);
	  gsl_spline_free (SM1P[j]);
	}
      delete[] AM1;  delete[] AM2;  delete[] AM3;  delete[] AM4;  delete[] AM5;  delete[] AM6; 
      delete[] AM7;  delete[] SM1;  delete[] SM2;  delete[] SM3;  delete[] SM4;  delete[] SM5; 
      delete[] SM6;  delete[] SM7;  delete[] AM8;  delete[] AM9;  delete[] SM8;  delete[] SM9;
      delete[] AM1P; delete[] AM3P; delete[] AM4P; delete[] SM1P; delete[] SM3P; delete[] SM4P;

      gsl_matrix_free (Ptor0); gsl_matrix_free (Qtor0); 
      gsl_matrix_free (Ptor1); gsl_matrix_free (Qtor1);
      delete[] Pm; delete[] Qm;
    }
}


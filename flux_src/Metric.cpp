// Metric.cpp

#include "Flux.h"

// #######################################
// Function to evaluate metric information
// #######################################
void Flux::Metric ()
{
  FILE *file1  = OpenFile ("Stage2/M1.out");
  FILE *file2  = OpenFile ("Stage2/M2.out");
  FILE *file3  = OpenFile ("Stage2/M3.out");
  FILE *file4  = OpenFile ("Stage2/M4.out");
  FILE *file5  = OpenFile ("Stage2/M5.out");
  FILE *file6  = OpenFile ("Stage2/M6.out");
  FILE *file7  = OpenFile ("Stage2/M7.out");
  FILE *file8  = OpenFile ("Stage2/M8.out"); 
  FILE *file9  = OpenFile ("Stage2/M9.out");
  FILE *file10 = OpenFile ("Stage2/M1P.out");
  FILE *file11 = OpenFile ("Stage2/M3P.out");
  FILE *file12 = OpenFile ("Stage2/M4P.out");

  double *M1   = new double [K];
  double *M2   = new double [K];
  double *M3   = new double [K];
  double *M4   = new double [K];
  double *M5   = new double [K];
  double *M6   = new double [K];
  double *M7   = new double [K];
  double *M8   = new double [K];
  double *M9   = new double [K]; 
  gsl_matrix *FM1  = gsl_matrix_alloc (J, K/2);
  gsl_matrix *FM2  = gsl_matrix_alloc (J, K/2);
  gsl_matrix *FM3  = gsl_matrix_alloc (J, K/2);
  gsl_matrix *FM4  = gsl_matrix_alloc (J, K/2);
  gsl_matrix *FM5  = gsl_matrix_alloc (J, K/2);
  gsl_matrix *FM6  = gsl_matrix_alloc (J, K/2);
  gsl_matrix *FM7  = gsl_matrix_alloc (J, K/2);
  gsl_matrix *FM8  = gsl_matrix_alloc (J, K/2);
  gsl_matrix *FM9  = gsl_matrix_alloc (J, K/2);
  gsl_matrix *FM1P = gsl_matrix_alloc (J, K/2);
  gsl_matrix *FM3P = gsl_matrix_alloc (J, K/2);
  gsl_matrix *FM4P = gsl_matrix_alloc (J, K/2);

  for (int j = 0; j < J; j++)
    {
      for (int k = 0; k < K; k++)
	{
	  double r2   = Raxis*Raxis + rP[j] * gsl_matrix_get (R2, j, k);
	  double igr2 = gsl_matrix_get (igradr2, j, k) /epsa/epsa;
	  double grt  = gsl_matrix_get (gradrt,  j, k);

	  M1 [k] = r2;
	  M2 [k] = igr2/r2;
	  M3 [k] = igr2;
	  M4 [k] = r2*igr2;
	  M5 [k] = r2*r2*igr2;
	  M6 [k] = grt;
	  M7 [k] = r2*grt;
	  M8 [k] = 1./igr2;
	  M9 [k] = r2*r2;
	}
      
      FTCos (j, M1,  FM1);
      FTCos (j, M2,  FM2);
      FTCos (j, M3,  FM3);
      FTCos (j, M4,  FM4);
      FTCos (j, M5,  FM5);
      FTSin (j, M6,  FM6);
      FTSin (j, M7,  FM7);
      FTCos (j, M8,  FM8);
      FTCos (j, M9,  FM9);

      for (int m = 0; m < K/2; m++)
	{
	  fprintf (file1 , "%17.10e ", gsl_matrix_get (FM1,  j, m));
	  fprintf (file2 , "%17.10e ", gsl_matrix_get (FM2,  j, m));
	  fprintf (file3 , "%17.10e ", gsl_matrix_get (FM3,  j, m));
	  fprintf (file4 , "%17.10e ", gsl_matrix_get (FM4,  j, m));
	  fprintf (file5 , "%17.10e ", gsl_matrix_get (FM5,  j, m));
	  fprintf (file6 , "%17.10e ", gsl_matrix_get (FM6,  j, m));
	  fprintf (file7 , "%17.10e ", gsl_matrix_get (FM7,  j, m));
	}
      fprintf (file8, "%17.10e\n", gsl_matrix_get (FM8, j, 0));
      fprintf (file9, "%17.10e\n", gsl_matrix_get (FM9, j, 0));
      fprintf (file1, "\n"); fprintf (file2, "\n"); fprintf (file3, "\n"); 
      fprintf (file4, "\n"); fprintf (file5, "\n"); fprintf (file6, "\n");
      fprintf (file7, "\n"); 
    }

  for (int j = 0; j < J; j++)
    for (int m = 0; m < K/2; m++)
      {
	gsl_matrix_set (FM1P, j, m, Interpolate1 (J, rP, FM1, rP[j], m, 1));
	gsl_matrix_set (FM3P, j, m, Interpolate1 (J, rP, FM3, rP[j], m, 1));
	gsl_matrix_set (FM4P, j, m, Interpolate1 (J, rP, FM4, rP[j], m, 1));
      }	

  double zero = 0.;
  for (int j = 0; j < J; j++)
    {
      fprintf (file11, "%17.10e\n", gsl_matrix_get (FM3P, j, 0));
      fprintf (file12, "%17.10e\n", gsl_matrix_get (FM4P, j, 0));

      for (int m = 0; m < K/2; m++)
	fprintf (file10, "%17.10e ", gsl_matrix_get (FM1P, j, m));
      fprintf (file10, "\n");
    }

  double rx = 1.;
  for (int m = 0; m < K/2; m++)
    {
      fprintf (file1,  "%17.10e ", Extrapolate2 (J, rP, FM1,  rx, m));
      fprintf (file2,  "%17.10e ", Extrapolate2 (J, rP, FM2,  rx, m));
      fprintf (file3,  "%17.10e ", Extrapolate2 (J, rP, FM3,  rx, m));
      fprintf (file4,  "%17.10e ", Extrapolate2 (J, rP, FM4,  rx, m));
      fprintf (file5,  "%17.10e ", Extrapolate2 (J, rP, FM5,  rx, m));
      fprintf (file6,  "%17.10e ", Extrapolate2 (J, rP, FM6,  rx, m));
      fprintf (file7,  "%17.10e ", Extrapolate2 (J, rP, FM7,  rx, m));
      fprintf (file10, "%17.10e ", Extrapolate2 (J, rP, FM1P, rx, m));
    }
  fprintf (file1, "\n");  fprintf (file2,  "\n"); fprintf (file3, "\n"); 
  fprintf (file4, "\n");  fprintf (file5,  "\n"); fprintf (file6, "\n");
  fprintf (file7, "\n");  fprintf (file10, "\n");  
  fprintf (file8,  "%17.10e\n", Extrapolate2 (J, rP, FM8,  rx, 0));
  fprintf (file9,  "%17.10e\n", Extrapolate2 (J, rP, FM9,  rx, 0));
  fprintf (file11, "%17.10e\n", Extrapolate2 (J, rP, FM3P, rx, 0));
  fprintf (file12, "%17.10e\n", Extrapolate2 (J, rP, FM4P, rx, 0));
  fclose (file1);  fclose (file2);  fclose (file3);  
  fclose (file4);  fclose (file5);  fclose (file6);  
  fclose (file7);  fclose (file8);  fclose (file9);  
  fclose (file10); fclose (file11); fclose (file12);
 
  delete[] M1;  delete[] M2;  delete[] M3;  
  delete[] M4;  delete[] M5;  delete[] M6;  
  delete[] M7;  delete[] M8;  delete[] M9;

  gsl_matrix_free (FM1);  gsl_matrix_free (FM2);  gsl_matrix_free (FM3); 
  gsl_matrix_free (FM4);  gsl_matrix_free (FM5);  gsl_matrix_free (FM6); 
  gsl_matrix_free (FM7);  gsl_matrix_free (FM8);  gsl_matrix_free (FM9); 
  gsl_matrix_free (FM1P); gsl_matrix_free (FM3P); gsl_matrix_free (FM4P);
}


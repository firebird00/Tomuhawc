// Log.cpp
// Various routines for logging information 

#include "Tomuhawc.h"

// #################################
// Function to log coupling matrices
// #################################
void Thawc::LogMatrices ()
{
  FILE *file = OpenFile ("Stage3/Bmat.out");
  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
	fprintf(file, "%e ", -fabs(gsl_matrix_get (bmat, i, j)));
      fprintf (file, "\n");
    }
  fclose (file);
  file = OpenFile ("Stage3/Cmat.out");
  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
	fprintf(file, "%e ", -fabs(gsl_matrix_get (cmat, i, j)));
      fprintf (file, "\n");
    }
  fclose (file);
  file = OpenFile ("Stage3/Dmat.out");
  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
	fprintf(file, "%e ", -fabs(gsl_matrix_get (dmat, i, j)));
      fprintf (file, "\n");
    }
  fclose (file);
  file = OpenFile ("Stage3/Emat.out");
  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
	fprintf(file, "%e ", -fabs(gsl_matrix_get (emat, i, j)));
      fprintf (file, "\n");
    }
  fclose (file);
  file = OpenFile ("Stage3/Lmat.out");
  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
	fprintf(file, "%e ", -fabs(gsl_matrix_get (Lmat, i, j)));
      fprintf (file, "\n");
    }
  fclose (file);
  file = OpenFile ("Stage3/Mmat.out");
  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
	fprintf(file, "%e ", -fabs(gsl_matrix_get (Mmat, i, j)));
      fprintf (file, "\n");
    }
  fclose (file);
  file = OpenFile ("Stage3/Nmat.out");
  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
	fprintf(file, "%e ", -fabs(gsl_matrix_get (Nmat, i, j)));
      fprintf (file, "\n");
    }
  fclose (file);
  file = OpenFile ("Stage3/Pmat.out");
  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
	fprintf(file, "%e ", -fabs(gsl_matrix_get (Pmat, i, j)));
      fprintf (file, "\n");
    }
  fclose (file);
}

// #########################################
// Function to log multiple solution vectors
// #########################################
void Thawc::LogVector (double r, gsl_matrix *YY)
{
  FILE *logfile = OpenFilea ("Stage3/logfile");
  fprintf (logfile, "Solution vectors at r = %20.13e\n", r);
  fprintf (logfile, "m_cen  :");
  for (int i = 0; i < dim; i++)
    fprintf (logfile, " %+3d        ", mpol[i]);
  for (int i = 0; i < vac; i++)
    fprintf (logfile, " %+3d        ", Mres[i]);
  fprintf (logfile, "\n");
  for (int k = 0; k < dim; k++)
    {
      fprintf (logfile, "m = %+3d: ", mpol[k]);
      for (int i = 0; i < dim; i++)
	{	
	  double y =  gsl_matrix_get (YY, k, i);
	  y = (fabs(y) < Eta) ? 0. : y;
	  fprintf (logfile, "%11.4e ", y);	
	} 
      for (int i = 0; i < vac; i++)
	{
	  double y =  gsl_matrix_get (YY, k, dim+i);
	  y = (fabs(y) < Eta) ? 0. : y;
	  fprintf (logfile, "%11.4e ", y);	
	}
      fprintf (logfile, "\n");
    }
  for (int k = dim; k < 2*dim; k++)
    {
      fprintf (logfile, "m = %+3d: ", mpol[k-dim]);
      for (int i = 0; i < dim; i++)
	{	
	  double y =  gsl_matrix_get (YY, k, i);
	  y = (fabs(y) < Eta) ? 0. : y;
	  fprintf (logfile, "%11.4e ", y);	
	}
      for (int i = 0; i < vac; i++)
	{
	  double y =  gsl_matrix_get (YY, k, dim+i);
	  y = (fabs(y) < Eta) ? 0. : y;
	  fprintf (logfile, "%11.4e ", y);	
	}
        fprintf (logfile, "\n");
    }
  fclose (logfile);
}

// #########################
// Function to log Y1-matrix
// #########################
void Thawc::LogY1 (gsl_matrix *Y1)
{
  FILE *logfile = OpenFilea ("Stage3/logfile");
  fprintf (logfile, "Edge boundary condition check:\n");
  fprintf (logfile, "m_res  :");
  for (int l = 0; l < vac; l++)
    fprintf (logfile, " %+3d        ", Mres[l]);
  fprintf (logfile, "\n");
  for (int j = 0; j < dim; j++)
    {
      fprintf (logfile, "m = %+3d: ", mpol[j]);
      for (int l = 0; l < vac; l++)
	{
	  double y = gsl_matrix_get (Y1, j, l);
	  y = (fabs(y) < Eta) ? 0. : y;
	  fprintf (logfile, "%11.4e ", y);
	}
      fprintf (logfile, "\n");
    }
  fclose (logfile);
}

// ##########################
// Function to log Fee-matrix
// ##########################
void Thawc::LogFee (gsl_matrix *Fee)
{
  FILE *logfile = OpenFilea ("Stage3/logfile");
  gsl_matrix *FeeT = gsl_matrix_alloc (vac, vac);
  gsl_matrix_transpose_memcpy (FeeT, Fee);
  gsl_matrix_sub (FeeT, Fee); 
  gsl_matrix_scale (FeeT, 0.5);

  fprintf (logfile, "\nFee-matrix:\n");
  fprintf (logfile, "mp     :");
  for (int l = 0; l < vac; l++)
    fprintf (logfile, " %+3d        ", Mres[l]);
  fprintf (logfile, "\n");
  for (int j = 0; j < vac; j++)
    {
      fprintf (logfile, "m = %+3d: ", Mres[j]);
      for (int l = 0; l < vac; l++)
	{
	  double y = gsl_matrix_get (Fee, j, l);
	  y = (fabs(y) < Eta) ? 0. : y;
	  fprintf (logfile, "%11.4e ", y);
	}
      fprintf (logfile, "\n");
    }
  fprintf (logfile, "\nFee-matrix residuals:\n");
  fprintf (logfile, "mp     :");
  for (int l = 0; l < vac; l++)
    fprintf (logfile, " %+3d        ", Mres[l]);
  fprintf (logfile, "\n");
  for (int j = 0; j < vac; j++)
    {
      fprintf (logfile, "m = %+3d: ", Mres[j]);
      for (int l = 0; l < vac; l++)
	{
	  double y = gsl_matrix_get (FeeT, j, l);
	  y = (fabs(y) < Eta) ? 0. : fabs(y);
	  fprintf (logfile, "%11.4e ", y);
	}
      fprintf (logfile, "\n");
    }

  gsl_matrix_free (FeeT);
  fclose (logfile); 

  logfile = OpenFile ("Stage3/Fee.out");
  for (int l = 0; l < vac; l++)
    for (int j = 0; j < vac; j++)
      fprintf (logfile, "%e ", gsl_matrix_get (Fee, l, j));
  fprintf (logfile, "\n");
  fclose (logfile);
}

// #########################
// Function to log Ee-matrix
// #########################
void Thawc::LogEe (gsl_matrix *Ee, int interactive)
{
  gsl_matrix *EeT = gsl_matrix_alloc (num, num);
  gsl_matrix_transpose_memcpy (EeT, Ee);
  gsl_matrix_sub              (EeT, Ee); 
  gsl_matrix_scale            (EeT, 0.5);
  
  printf ("\nEe-matrix:\n");
  printf ("mp     :");
  for (int l = 0; l < num; l++)
    printf (" %+3d        ", Mres[l]);
  printf ("\n");
  for (int j = 0; j < num; j++)
    {
      printf ("m = %+3d: ", Mres[j]);
      for (int l = 0; l < num; l++)
	{
	  double y = gsl_matrix_get (Ee, j, l);
	  y = (fabs(y) < Eta) ? 0. : y;
	  printf ("%11.4e ", y);
	}
      printf ("\n");
    }
  printf ("\n");
  
  printf ("Ee-matrix residuals:\n");
  printf ("mp     :");
  for (int l = 0; l < num; l++)
    printf (" %+3d        ", Mres[l]);
  printf ("\n");
  for (int j = 0; j < num; j++)
    {
      printf ("m = %+3d: ", Mres[j]);
      for (int l = 0; l < num; l++)
	{
	  double fac = (fabs (gsl_matrix_get (Ee, j, l)) > 1.)
	    ? fabs (gsl_matrix_get (Ee, j, l)) : 1.;
	  double y = gsl_matrix_get (EeT, j, l) /fac;
	  y = (fabs(y) < Eta) ? 0. : fabs(y);
	  printf ("%11.4e ", y);
	}
      printf ("\n");
    }
  printf ("\n");

  FILE* logfile;
  if (interactive)
    {	
      logfile = OpenFilea ("Stage3/logfile");
      
      fprintf (logfile, "\nEe-matrix:\n");
      fprintf (logfile, "mp     :");
      for (int l = 0; l < num; l++)
	fprintf (logfile, " %+3d        ", Mres[l]);
      fprintf (logfile, "\n");
      for (int j = 0; j < num; j++)
	{
	  fprintf (logfile, "m = %+3d: ", Mres[j]);
	  for (int l = 0; l < num; l++)
	    {
	      double y = gsl_matrix_get (Ee, j, l);
	      y = (fabs(y) < Eta) ? 0. : y;
	      fprintf (logfile, "%11.4e ", y);
	    }
	  fprintf (logfile, "\n");
	} 
      
      fprintf (logfile, "\nEee-matrix residuals:\n");
      fprintf (logfile, "mp     :");
      for (int l = 0; l < num; l++)
	fprintf (logfile, " %+3d        ", Mres[l]);
      fprintf (logfile, "\n");
      for (int j = 0; j < num; j++)
	{
	  fprintf (logfile, "m = %+3d: ", Mres[j]);
	  for (int l = 0; l < num; l++)
	    {
	      double fac = (fabs (gsl_matrix_get (Ee, j, l)) > 1.) 
		? fabs (gsl_matrix_get (Ee, j, l)) : 1.;
	      double y = gsl_matrix_get (EeT, j, l) /fac;
	      y = (fabs(y) < Eta) ? 0. : fabs(y);
	      fprintf (logfile, "%11.4e ", y);
	    }
	  fprintf (logfile, "\n");
	}
      fclose (logfile);
    }	
  gsl_matrix_free (EeT);

  logfile = OpenFile ("Stage3/Ee.out");
  for (int l = 0; l < num; l++)
    for (int j = 0; j < num; j++)
      fprintf (logfile, "%e ", gsl_matrix_get (Ee, l, j));
  fprintf (logfile, "\n");
  fclose (logfile);

  logfile = OpenFile ("Stage3/Delta.out");
  fprintf (logfile, "%d", num);
  for (int l = 0; l < num; l++)
    fprintf (logfile, " %d %e", Mres[l], gsl_matrix_get (Ee, l, l));
  fprintf (logfile, "\n");
  fclose (logfile);
}

// ##########################
// Function to log Foe-matrix
// ##########################
void Thawc::LogFoe (gsl_matrix *Foe)
{
  FILE *logfile = OpenFilea ("Stage3/logfile");

  fprintf (logfile, "\nFoe-matrix:\n");
  fprintf (logfile, "mp     :");
  for (int l = 0; l < vac; l++)
    fprintf (logfile, " %+3d        ", Mres[l]);
  fprintf (logfile, "\n");
  for (int j = 0; j < vac; j++)
    {
      fprintf (logfile, "m = %+3d: ", Mres[j]);
      for (int l = 0; l < vac; l++)
	{
	  double y = gsl_matrix_get (Foe, j, l);
	  y = (fabs(y) < Eta) ? 0. : y;
	  fprintf (logfile, "%11.4e ", y);
	}
      fprintf (logfile, "\n");
    }
  fclose (logfile);

  logfile = OpenFile ("Stage3/Foe.out");
  for (int l = 0; l < vac; l++)
    for (int j = 0; j < vac; j++)
      fprintf (logfile, "%e ", gsl_matrix_get (Foe, l, j));
  fprintf (logfile, "\n");
  fclose (logfile);
}

// #########################
// Function to log Gp-matrix
// #########################
void Thawc::LogGp (gsl_matrix *Gp, int interactive)
{
  printf ("Gp-matrix:\n");
  printf ("mp     :");
  for (int l = 0; l < num; l++)
    printf (" %+3d        ", Mres[l]);
  printf ("\n");
  for (int j = 0; j < num; j++)
    {
      printf ("m = %+3d: ", Mres[j]);
      for (int l = 0; l < num; l++)
	{
	  double y = gsl_matrix_get (Gp, j, l);
	  y = (fabs(y) < Eta) ? 0. : y;
	  printf ("%11.4e ", y);
	}
      printf ("\n");
    }
  printf ("\n");

  if (interactive)
    {
      FILE *logfile = OpenFilea ("Stage3/logfile");
      fprintf (logfile, "\nGp-matrix:\n");
      fprintf (logfile, "mp     :");
      for (int l = 0; l < num; l++)
	fprintf (logfile, " %+3d        ", Mres[l]);
      fprintf (logfile, "\n");
      for (int j = 0; j < num; j++)
	{
	  fprintf (logfile, "m = %+3d: ", Mres[j]);
	  for (int l = 0; l < num; l++)
	    {
	      double y = gsl_matrix_get (Gp, j, l);
	      y = (fabs(y) < Eta) ? 0. : y;
	      fprintf (logfile, "%11.4e ", y);
	    }
	  fprintf (logfile, "\n");
	}
      fclose (logfile);
    }
}
// ##########################
// Function to log Foo-matrix
// ##########################
void Thawc::LogFoo (gsl_matrix *Foo)
{
  FILE *logfile = OpenFilea ("Stage3/logfile");
  gsl_matrix *FooT = gsl_matrix_alloc (vac, vac);
  gsl_matrix_transpose_memcpy (FooT, Foo);
  gsl_matrix_sub              (FooT, Foo); 
  gsl_matrix_scale            (FooT, 0.5);

  fprintf (logfile, "\nFoo-matrix:\n");
  fprintf (logfile, "mp     :");
  for (int l = 0; l < vac; l++)
    fprintf (logfile, " %+3d        ", Mres[l]);
  fprintf (logfile, "\n");
  for (int j = 0; j < vac; j++)
    {
      fprintf (logfile, "m = %+3d: ", Mres[j]);
      for (int l = 0; l < vac; l++)
	{
	  double y = gsl_matrix_get (Foo, j, l);
	  y = (fabs(y) < Eta) ? 0. : y;
	  fprintf (logfile, "%11.4e ", y);
	}
      fprintf (logfile, "\n");
    }
  fprintf (logfile, "\nFoo-matrix residuals:\n");
  fprintf (logfile, "mp     :");
  for (int l = 0; l < num; l++)
    fprintf (logfile, " %+3d        ", Mres[l]);
  fprintf (logfile, "\n");
  for (int j = 0; j < num; j++)
    {
      fprintf (logfile, "m = %+3d: ", Mres[j]);
      for (int l = 0; l < num; l++)
	{
	  double y = gsl_matrix_get (FooT, j, l);
	  y = (fabs(y) < Eta) ? 0. : fabs(y);
	  fprintf (logfile, "%11.4e ", y);
	}
      fprintf (logfile, "\n");
    }
  fclose (logfile); 
  gsl_matrix_free (FooT);

  logfile = OpenFile ("Stage3/Foo.out");
  for (int l = 0; l < vac; l++)
    for (int j = 0; j < vac; j++)
      fprintf (logfile, "%e ", gsl_matrix_get (Foo, l, j));
  fprintf (logfile, "\n");
  fclose (logfile);
}

// #########################
// Function to log Eo-matrix
// #########################
void Thawc::LogEo (gsl_matrix *Eo, int interactive)
{
  gsl_matrix *EoT = gsl_matrix_alloc (num, num);
  gsl_matrix_transpose_memcpy (EoT, Eo);
  gsl_matrix_sub              (EoT, Eo); 
  gsl_matrix_scale            (EoT, 0.5);

  printf ("Eo-matrix:\n");
  printf ("mp     :");
  for (int l = 0; l < num; l++)
    printf (" %+3d        ", Mres[l]);
  printf ("\n");
  for (int j = 0; j < num; j++)
    {
      printf ("m = %+3d: ", Mres[j]);
      for (int l = 0; l < num; l++)
	{
	  double y = gsl_matrix_get (Eo, j, l);
	  y = (fabs(y) < Eta) ? 0. : y;
	  printf ("%11.4e ", y);
	}
      printf ("\n");
    }
  printf ("\n");

  printf ("Eo-matrix residuals:\n");
  printf ("mp     :");
  for (int l = 0; l < num; l++)
    printf (" %+3d        ", Mres[l]);
  printf ("\n");
  for (int j = 0; j < num; j++)
    {
      printf ("m = %+3d: ", Mres[j]);
      for (int l = 0; l < num; l++)
	{
	  double fac = (fabs (gsl_matrix_get (Eo, j, l)) > 1.) ? fabs (gsl_matrix_get (Eo, j, l)) : 1.;
	  double y = gsl_matrix_get (EoT, j, l) /fac;
	  y = (fabs(y) < Eta) ? 0. : fabs(y);
	  printf ("%11.4e ", y);
	}
      printf ("\n");
    }
  printf ("\n");

  FILE *logfile;
  if (interactive)
    {
      logfile = OpenFilea ("Stage3/logfile");

      fprintf (logfile, "\nEo-matrix:\n");
      fprintf (logfile, "mp     :");
      for (int l = 0; l < num; l++)
	fprintf (logfile, " %+3d        ", Mres[l]);
      fprintf (logfile, "\n");
      for (int j = 0; j < num; j++)
	{
	  fprintf (logfile, "m = %+3d: ", Mres[j]);
	  for (int l = 0; l < num; l++)
	    {
	      double y = gsl_matrix_get (Eo, j, l);
	      y = (fabs(y) < Eta) ? 0. : y;
	      fprintf (logfile, "%11.4e ", y);
	    }
	  fprintf (logfile, "\n");
	}
      
      fprintf (logfile, "\nEo-matrix residuals:\n");
      fprintf (logfile, "mp     :");
      for (int l = 0; l < num; l++)
	fprintf (logfile, " %+3d        ", Mres[l]);
      fprintf (logfile, "\n");
      for (int j = 0; j < num; j++)
	{
	  fprintf (logfile, "m = %+3d: ", Mres[j]);
	  for (int l = 0; l < num; l++)
	    {
	      double fac = (fabs (gsl_matrix_get (Eo, j, l)) > 1.)
		? fabs (gsl_matrix_get (Eo, j, l)) : 1.;
	      double y = gsl_matrix_get (EoT, j, l) /fac;
	      y = (fabs(y) < Eta) ? 0. : fabs(y);
	      fprintf (logfile, "%11.4e ", y);
	    }
	  fprintf (logfile, "\n");
	}
      fclose (logfile);
    }
  gsl_matrix_free (EoT);

  logfile = OpenFile ("Stage3/Eo.out");
  for (int l = 0; l < num; l++)
    for (int j = 0; j < num; j++)
      fprintf (logfile, "%e ", gsl_matrix_get (Eo, l, j));
  fprintf (logfile, "\n");
  fclose (logfile);
}

// ##########################
// Function to log Feo-matrix
// ##########################
void Thawc::LogFeo (gsl_matrix *Feo)
{
  FILE *logfile = OpenFilea ("Stage3/logfile");

  fprintf (logfile, "\nFeo-matrix:\n");
  fprintf (logfile, "mp     :");
  for (int l = 0; l < vac; l++)
    fprintf (logfile, " %+3d        ", Mres[l]);
  fprintf (logfile, "\n");
  for (int j = 0; j < vac; j++)
    {
      fprintf (logfile, "m = %+3d: ", Mres[j]);
      for (int l = 0; l < vac; l++)
	{
	  double y = gsl_matrix_get (Feo, j, l);
	  y = (fabs(y) < Eta) ? 0. : y;
	  fprintf (logfile, "%11.4e ", y);
	}
      fprintf (logfile, "\n");
    }
  fclose (logfile);

  logfile = OpenFile ("Stage3/Feo.out");
  for (int l = 0; l < vac; l++)
    for (int j = 0; j < vac; j++)
      fprintf (logfile, "%e ", gsl_matrix_get (Feo, l, j));
  fprintf (logfile, "\n");
  fclose (logfile);
}

// ########################
// Function to log G-matrix
// ########################
void Thawc::LogG (gsl_matrix *G, gsl_matrix *Gp, int interactive)
{
  gsl_matrix *GT = gsl_matrix_alloc (num, num);
  gsl_matrix_transpose_memcpy (GT, Gp);
  gsl_matrix_sub              (GT, G); 
  gsl_matrix_scale            (GT, 0.5);

  printf ("G-matrix:\n");
  printf ("mp     :");
  for (int l = 0; l < num; l++)
    printf (" %+3d        ", Mres[l]);
  printf ("\n");
  for (int j = 0; j < num; j++)
    {
      printf ("m = %+3d: ", Mres[j]);
      for (int l = 0; l < num; l++)
	{
	  double y = gsl_matrix_get (G, j, l);
	  y = (fabs(y) < Eta) ? 0. : y;
	  printf ("%11.4e ", y);
	}
      printf ("\n");
    }
  printf ("\n");

  printf ("G-matrix residuals:\n");
  printf ("mp     :");
  for (int l = 0; l < num; l++)
    printf (" %+3d        ", Mres[l]);
  printf ("\n");
  for (int j = 0; j < num; j++)
    {
      printf ("m = %+3d: ", Mres[j]);
      for (int l = 0; l < num; l++)
	{
	  double fac = (fabs (gsl_matrix_get (G, j, l)) > 1.)
	    ? fabs (gsl_matrix_get (G, j, l)) : 1.;
	  double y = gsl_matrix_get (GT, j, l) /fac;
	  y = (fabs(y) < Eta) ? 0. : y;
	  printf ("%11.4e ", y);
	}
      printf ("\n");
    }
  printf ("\n");

  FILE *logfile;
  if (interactive)
    {	
      logfile = OpenFilea ("Stage3/logfile");

      fprintf (logfile, "\nG-matrix:\n");
      fprintf (logfile, "mp     :");
      for (int l = 0; l < num; l++)
	fprintf (logfile, " %+3d        ", Mres[l]);
      fprintf (logfile, "\n");
      for (int j = 0; j < num; j++)
	{
	  fprintf (logfile, "m = %+3d: ", Mres[j]);
	  for (int l = 0; l < num; l++)
	    {
	      double y = gsl_matrix_get (G, j, l);
	      y = (fabs(y) < Eta) ? 0. : y;
	      fprintf (logfile, "%11.4e ", y);
	    }
	  fprintf (logfile, "\n");
	}
      
      fprintf (logfile, "\nG-matrix residuals:\n");
      fprintf (logfile, "mp     :");
      for (int l = 0; l < num; l++)
	fprintf (logfile, " %+3d        ", Mres[l]);
      fprintf (logfile, "\n");
      for (int j = 0; j < num; j++)
	{
	  fprintf (logfile, "m = %+3d: ", Mres[j]);
	  for (int l = 0; l < num; l++)
	    {
	      double fac = (fabs (gsl_matrix_get (G, j, l)) > 1.)
		? fabs (gsl_matrix_get (G, j, l)) : 1.;
	      double y = gsl_matrix_get (GT, j, l) /fac;
	      y = (fabs(y) < Eta) ? 0. : y;
	      fprintf (logfile, "%11.4e ", y);
	    }
	  fprintf (logfile, "\n");
	}
      fclose (logfile);
    }
  gsl_matrix_free (GT);

  logfile = OpenFile ("Stage3/G.out");
  for (int l = 0; l < num; l++)
    for (int j = 0; j < num; j++)
      fprintf (logfile, "%e ", gsl_matrix_get (G, l, j));
  fprintf (logfile, "\n");
  fclose (logfile);

  logfile = OpenFile ("Stage3/Gp.out");
  for (int l = 0; l < num; l++)
    for (int j = 0; j < num; j++)
      fprintf (logfile, "%e ", gsl_matrix_get (Gp, l, j));
  fprintf (logfile, "\n");
  fclose (logfile);
}

// #####################################
// Function to open new file for writing
// #####################################
FILE* Thawc::OpenFile (char *filename)
{
  FILE *file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("Thawc::OpenFile: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// ##########################################
// Function to open existing file for writing
// ##########################################
FILE* Thawc::OpenFilea (char *filename)
{ 
  FILE *file = fopen (filename, "a");
  if (file == NULL) 
    {
      printf ("Thawc::OpenFilea: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// ##########################################
// Function to open existing file for reading
// ##########################################
FILE* Thawc::OpenFiler (char *filename)
{
  FILE *file = fopen (filename, "r");
  if (file == NULL) 
    {
      printf ("Thawc::OpenFiler: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

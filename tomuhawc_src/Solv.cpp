// Solv.cpp

#include "Tomuhawc.h"

// #########################################################################
// Function to solve stability problem
//
// interactive ... 1/0 if function called in interactive/noninteractive mode
//
// #########################################################################
int Thawc::Solv (int interactive)
{
  // %%%%%%%%%%%%%%%%%%%%%%
  // Initialize calculation
  // %%%%%%%%%%%%%%%%%%%%%%

  // Initialize output files
  char *_adape = "Stage3/adape.out";
  char *_adapo = "Stage3/adapo.out";
  char *_solne = "Stage3/solne.out";
  char *_solno = "Stage3/solno.out";
  char *_edgee = "Stage3/Yea.out";
  char *_edgeo = "Stage3/Yoa.out";

  system ("rm Stage3/*.out");
 
  // ++++++++++++++++++++++++
  // Read in equilibrium data
  // ++++++++++++++++++++++++
  printf ("\nRead in equilibrium data:\n");
  Equilibrium (0, interactive);

  // ++++++++++++++++++++++++++
  // Log calculation parameters
  // ++++++++++++++++++++++++++
  if (interactive) Status (0);

  // ++++++++++++++++++++++
  // Find resonant surfaces
  // ++++++++++++++++++++++
  printf ("Find rational surfaces:\n");
  Resonant (0, interactive);

  // Abort calculation if no resonant surfaces found
  if (vac == 0)
    {
      printf ("No resonant surfaces found. Abort calculation.\n");
      if (interactive)
	{
	  FILE *logfile = OpenFilea ("Stage3/logfile");
	  fprintf (logfile, "No resonant surfaces found. Abort calculation.\n");
	  fclose (logfile);
	}
      Equilibrium (1, interactive);
      Resonant    (1, interactive);
      return 1;
    } 

  // Abort calculation if resonant surfaces in core
  if (Rres[0] < r1*rb)
    {
      printf ("Resonant surfaces in core:\n");
      if (interactive)
	{
	  FILE *logfile = OpenFilea ("Stage3/logfile");
	  fprintf (logfile, "Resonant surfaces in core.\n");
	  fclose (logfile);
	}
      Equilibrium (1, interactive);
      Resonant    (1, interactive);
      return 1;
    }
  
  // Abort calculation if Mercier indices not in range 0. < DI < 1.
  for (int j = 0; j < vac; j++)
    {
      double DI = 0.25 - EFres[j] - HHres[j];
      if (DI <= 0. || DI >= 1.)
	{
	  printf ("Mercier indices out of range. Abort calculation.\n");
	  if (interactive)
	    {	
	      FILE *logfile = OpenFilea ("Stage3/logfile");
	      fprintf (logfile, "Mercier indices out of range. Abort calculation.\n");
	      fclose (logfile);
	      Equilibrium (1, interactive);
	      Resonant    (1, interactive);
	    }	
	  return 1;
	}
    }

  // +++++++++++++++
  // Allocate memory
  // +++++++++++++++
  dim  = vac + 2*side;
  dim1 = 2*dim;
  dimN = dim1*(dim+vac);

  mpol = new int[dim];
  for (int i = 0; i < dim; i++)
    mpol[i] = Mres[0] - side + off + i;

  bmat = gsl_matrix_alloc (dim,  dim);
  cmat = gsl_matrix_alloc (dim,  dim);
  dmat = gsl_matrix_alloc (dim,  dim);
  emat = gsl_matrix_alloc (dim,  dim);
  Lmat = gsl_matrix_alloc (dim,  dim);
  Mmat = gsl_matrix_alloc (dim,  dim);
  Nmat = gsl_matrix_alloc (dim,  dim);
  Pmat = gsl_matrix_alloc (dim,  dim);
  Vmat = gsl_matrix_alloc (dim1, dim1);
  //Vacuum ();
  if (interactive)
    {
      Couple (rb);
      LogMatrices ();
    }
  
  gsl_matrix *YY   = gsl_matrix_calloc (dim1, dim+vac);
  gsl_matrix *Psi  = gsl_matrix_calloc (vac,  dim+vac);
  gsl_matrix *dPsi = gsl_matrix_calloc (vac,  dim+vac);
  gsl_matrix *x    = gsl_matrix_calloc (dim,  vac);
  gsl_matrix *Fee  = gsl_matrix_calloc (vac,  vac);
  gsl_matrix *Foe  = gsl_matrix_calloc (vac,  vac);
  gsl_matrix *Foo  = gsl_matrix_calloc (vac,  vac);
  gsl_matrix *Feo  = gsl_matrix_calloc (vac,  vac);
  gsl_matrix *Ee   = gsl_matrix_calloc (num,  num);
  gsl_matrix *Eo   = gsl_matrix_calloc (num,  num);
  gsl_matrix *Gp   = gsl_matrix_calloc (num,  num);
  gsl_matrix *G    = gsl_matrix_calloc (num,  num);

  // ++++++++++++++++++++++++++++++
  // Output calculation information
  // ++++++++++++++++++++++++++++++
  FILE *info = OpenFile ("Stage3/info.out");
  fprintf (info, "%d %d %d %d %e\n", dim, num, vac, off, rb);
  fclose (info);

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Integrate tearing parity solution vectors from magnetic axis to wall
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Initialize tearing partity solution vectors at magnetic axis
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  double r = r0; step = 0; func = 0; nflag = 0;
  printf ("Initialize tearing parity solution vectors at magnetic axis:\n");

  Launch    (r, YY);
  if (interactive) LogVector (r, YY);
  
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Integrate tearing parity solution vectors from magnetic axis to edge of plasma core 
  // while performing nfix logarithmically spaced fixups.
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  int nf = (nfix == 0) ? 0 : int ((r1*rb) * double(nfix));
  Segment_Fixup (r, r1*rb, nf, YY, Psi, dPsi, 1, 0, 0, interactive, _adape, _solne, _edgee);
  
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Integrate tearing parity solution vectors across rational surfaces
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  printf ("Integrate tearing parity solution vectors across rational surfaces:\n");

  double re;
  for (int j = 0; j < vac; j++)
    {
      // Integrate tearing parity solution vectors to just before jth rational surface
      // while performing nfix equally spaced fixups
      re = Rres[j] - del;
      nf = (nfix == 0) ? 0 : int ((re - r) * double (nfix));
      Segment_Fixup (r, re, nf, YY, Psi, dPsi, 0, flg, 0, interactive, _adape, _solne, _edgee);
 
      // Jump tearing parity solution vectors across jth rational surface
      double q = gsl_spline_eval (Sq,  Rres[j], Aq);
      printf ("Rational surface: r/b = %20.13e  m_res = %4d  s = %11.4e  DI = %11.4e  res = %11.4e\n",
	      Rres[j]/rb, Mres[j], Sres[j], 0.25-EFres[j]-HHres[j], fabs(q - Qres[j])); 
      if (interactive)
	{
	  FILE *logfile = OpenFilea ("Stage3/logfile");
	  fprintf (logfile, "Rational surface: r/b = %20.13e  m_res = %4d  s = %11.4e  DI = %11.4e  res = %11.4e\n",
		   Rres[j]/rb, Mres[j], Sres[j], 0.25-EFres[j]-HHres[j], fabs(q - Qres[j])); 
	  fclose (logfile);
	}
      Jump (r, j, YY, Psi, dPsi, 1, interactive);
      if (interactive) LogVector (r, YY);
    }
 
  // +++++++++++++++++++++++++++++++++++++++++++++++++
  // Integrate tearing parity solution vectors to wall
  // while performing nfix equally spaced fixups
  // +++++++++++++++++++++++++++++++++++++++++++++++++
  printf ("Integrate tearing parity solution vectors to wall:\n");
 
  nf = (nfix == 0) ? 0 : int ((ra - r) * double (nfix));
  Segment_Fixup (r, ra, nf, YY, Psi, dPsi, 0, flg, 1, interactive, _adape, _solne, _edgee);
  printf ("Wall: r/b = %11.4e\n", ra/rb);
  printf ("step = %10d  func = %10d\n", step, func);
  FILE *logfile = OpenFilea ("Stage3/logfile");
  fprintf (logfile, "step = %10d  func = %10d\n", step, func);
  fclose (logfile);

  // +++++++++++++++++++++++++++++++++
  // Apply boundary conditions at wall
  // +++++++++++++++++++++++++++++++++
  Boundary (YY, x, interactive);

  // ++++++++++++++++++++
  // Calculate Fee-matrix
  // ++++++++++++++++++++
  CalcFee (Psi, x, Fee, YY);

  // ++++++++++++++++++++
  // Calculate Foe-matrix
  // ++++++++++++++++++++
  CalcFoe (dPsi, x, Foe);

  // ++++++++++++++++++++++++++++
  // Calculate Ee and Gp matrices
  // ++++++++++++++++++++++++++++
  CalcEe (Fee, Foe, Ee, Gp);

  // ++++++++++++++++++++++++
  // Calculate eigenfunctions
  // ++++++++++++++++++++++++
  CalcEigFeea (x);
  CalcEigEea  (Ee);
  if (interactive)
    {
      CalcEigFee (x);
      CalcEigEe  (Ee);
    }

  if (twist && !nflag)
    {
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      // Integrate twisting parity solution vectors from magnetic axis to wall
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // Initialize twisting partity solution vectors at magnetic axis
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      r = r0; step = 0; func = 0;
      printf ("Initialize twisting parity solution vectors at magnetic axis:\n");
      
      Launch (r, YY);
      if (interactive) LogVector (r, YY);
      
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // Integrate twisting parity solution vectors from magnetic axis to edge of plasma core 
      // while performing nfix logarimaically spaced fixups.
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      nf = (nfix == 0) ? 0 : int ((r1*rb) * double(nfix));
      Segment_Fixup (r, r1*rb, nf, YY, Psi, dPsi, 1, 0, 0, interactive, _adapo, _solno, _edgeo);
      
      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // Integrate twisting parity solution vectors across rational surfaces
      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      printf ("Integrate twisting parity solution vectors across rational surfaces:\n");
      
      for (int j = 0; j < vac; j++)
	{
	  // Integrate twisting parity solution vectors to just before jth rational surface
	  // while performing nfix equally spaced fixups.
	  re = Rres[j] - del;
	  nf = (nfix == 0) ? 0 : int ((re - r) * double (nfix));
	  Segment_Fixup (r, re, nf, YY, Psi, dPsi, 0, flg, 0, interactive, _adapo, _solno, _edgeo);
	  
	  // Jump tearing parity solution vectors across jth rational surface
	  double q = gsl_spline_eval (Sq,  Rres[j], Aq);
	  printf ("Rational surface: r/b = %20.13e  m_res = %4d  s = %11.4e  DI = %11.4e  DR = %11.4e  res = %11.4e\n",
		  Rres[j]/rb, Mres[j], Sres[j], 0.25-EFres[j]-HHres[j],
		  6.35*(EFres[j]+HHres[j]*HHres[j])/(0.5+sqrt(0.25-EFres[j]-HHres[j])-HHres[j]), fabs(q - Qres[j])); 
	  if (interactive)
	    {
	      FILE *logfile = OpenFilea ("Stage3/logfile");
	      fprintf (logfile, "Rational surface: r/b = %20.13e  m_res = %4d  s = %11.4e  DI = %11.4e  DR = %11.4e\n",
		       Rres[j]/rb, Mres[j], Sres[j], 0.25-EFres[j]-HHres[j], 
		       6.35*(EFres[j]+HHres[j]*HHres[j])/(0.5+sqrt(0.25-EFres[j]-HHres[j])-HHres[j])); 
	      fclose (logfile);
	    }
	  Jump (r, j, YY, Psi, dPsi, 0, interactive);
	  if (interactive) LogVector (r, YY);
	}
      
      // ++++++++++++++++++++++++++++++++++++++++++++++++++
      // Integrate twisting parity solution vectors to wall
      // while performing nfix equally spaced fixups.
      // ++++++++++++++++++++++++++++++++++++++++++++++++++
      printf ("Integrate twisting parity solution vectors to wall:\n");
      
      nf = (nfix == 0) ? 0 : int ((ra - r) * double (nfix));
      Segment_Fixup (r, ra, nf, YY, Psi, dPsi, 0, flg, 1, interactive, _adapo, _solno, _edgeo);
       printf ("Wall: r/b = %11.4e\n", ra/rb);
      printf ("step = %10d  func = %10d\n", step, func);
      logfile = OpenFilea ("Stage3/logfile");
      fprintf (logfile, "step = %10d  func = %10d\n", step, func);
      fclose (logfile);
      
      // +++++++++++++++++++++++++++++++++
      // Apply boundary conditions at wall
      // +++++++++++++++++++++++++++++++++
      Boundary (YY, x, interactive);
      
      // ++++++++++++++++++++
      // Calculate Foo-matrix
      // ++++++++++++++++++++
      CalcFoo (Psi, x, Foo, YY);
      
      // ++++++++++++++++++++
      // Calculate Feo-matrix
      // ++++++++++++++++++++
      CalcFeo (dPsi, x, Feo);

      // +++++++++++++++++++++++++++
      // Calculate Eo and G matrices
      // +++++++++++++++++++++++++++
      CalcEo (Foo, Feo, Eo, G);

      // ++++++++++++++++++++++++
      // Calculate eigenfunctions
      // ++++++++++++++++++++++++
      CalcEigFooa (x);
      CalcEigEoa  (Eo);
      if (interactive)
	{
	  CalcEigFoo (x);
	  CalcEigEo  (Eo);
	}	
    }

  // %%%%%%%%%%%%%%%%%%%%%
  // Output stability data
  // %%%%%%%%%%%%%%%%%%%%%
  LogEe (Ee, interactive);
  if (!nflag)
    LogGp (Gp, interactive);
  if (twist && !nflag)
    {
      LogEo (Eo,    interactive);
      LogG  (G, Gp, interactive); 
    }
  if (interactive)
    {
      LogFee (Fee);
      if (!nflag) 
	LogFoe (Foe);
      if (twist && !nflag)
	{
	  LogFoo (Foo);
	  LogFeo (Feo);
	}
    }

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Output poloidal mode numbers
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FILE *mode = OpenFile ("Stage3/modes.out");
  for (int j = 0; j < vac; j++)
    for (int k = 0; k < dim; k++)
      fprintf (mode, "%d ", mpol[k]);
  fprintf (mode, "\n"); 
  fclose (mode); 
 
  // %%%%%%%%%%%%%%%%%%%%
  // Clean up calculation
  // %%%%%%%%%%%%%%%%%%%%
  gsl_matrix_free (YY);   gsl_matrix_free (Psi);  gsl_matrix_free (dPsi);
  gsl_matrix_free (x);    gsl_matrix_free (Fee);  gsl_matrix_free (Feo);
  gsl_matrix_free (Foe);  gsl_matrix_free (Foo);  gsl_matrix_free (G);
  gsl_matrix_free (Ee);   gsl_matrix_free (Eo);   gsl_matrix_free (Gp);  
  gsl_matrix_free (bmat); gsl_matrix_free (cmat);
  gsl_matrix_free (dmat); gsl_matrix_free (emat);
  gsl_matrix_free (Lmat); gsl_matrix_free (Mmat);
  gsl_matrix_free (Nmat); gsl_matrix_free (Pmat);
  gsl_matrix_free (Vmat);
  
  delete[] mpol;

  Equilibrium (1, interactive);
  Resonant    (1, interactive);

  return 0;
}




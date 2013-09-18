// main.cpp

// Program to call CHEASE recursively to obtain equilibrium with given 
// qc = Q(PSIC), qb = Q(PSIB), and betaN
//
// Equilibrium specified in profile.in
//
// Increment data specified in profile.inc
// 
// Plasma boundary:
//
// R = 1. + ASPCT cos[t + TRIANG sin(t)]
//
// Z =      ASPCT ELONG sin(t)    0 <= t <= 2pi  
//
// Profiles:
//
// PSI - polidal flux
//
// s = (1 - PSI/PSIC)^1/2
// 
// PSI = PSIC, s = 0  at magnetic axis
// PSI = PSIB, s = sb at plasma boundary (0 < sb <= 1)
// PSI = 0,    s = 1  at wall
//
// dP/dPSI   = P0 [1 - (s/sb)^2 MUP]^NUP   0 <= s <= sb
//
// dP/dPSI   = 0                           sb < s <= 1
//
// T dT/dPSI = T0 [1 - (s/sb)^2 MU ]^NU    0 <= s <= sb
//
// T dT/dPSI = 0                           sb < s <= 1
//
// CHEASE determines T0 automatically to fix qc. 
// Program iterates P0 and NU to fix betaN and qb, respectively.
//
// NBPTS - number of points specified on boundary.
// NPTS  - number of points specified in profiles.
// NS    - number of radial points in CHEASE calculation
// NT    - number of angular points in CHEASE calculation
// NRBOX - number of R points in CHEASE output
// NZBOX - number of Z points in CHEASE output
//
// Output in folder Stage1/ (see README)
//
// If called without argument then program outputs usage message and exits.
// If called with argument  0 then program runs single PEST equilibrium calculation and exits.
// If called with argument -1 then program iterates once without stability calculation and exits.
// If called with argument -2 then program writes CHEASE namelist file and exits.
// If called with argument +N then program iterates N times with stability calculation and exits.

// Header files
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Function prototypes 
void WriteNamelist (int NIDEAL, double QC, int NS, int NT, int NRBOX, int NZBOX);
void WriteProfile  (int NBPTS, double ASPCT, double ELONG, double TRIANG, 
		    int NPTS, double SB, double P0, double MUP, double NUP, double MU, double NU);
FILE* OpenFile  (char *filename);
FILE* OpenFiler (char *filename);

// Equilibrium search accuracy control parameters
double epsb  = 2.e-5;   // Maximum tolerable difference between requested and actual BETAN
double epsq  = 2.e-6;   // Maximum tolerable difference between requested and actual QB
int    itmax = 20;      // Maximum allowed number of iterations

// #############
// Main function
// #############
int main (int argc, char *argv[])
{
  // ########################################
  // Determine manner in which program called
  // ########################################
  int control = 0;
  if (argc > 1)
    control = atoi(argv[1]);
  else
    {
      printf ("\nUsage:\n");
      printf ("profile    - print this message and exit\n");
      printf ("profile  0 - perform single PEST equilibrium calculation and exit\n");
      printf ("profile -1 - perform single equilibrium iteration and exit\n");
      printf ("profile -2 - write CHEASE namelist file and exit\n");
      printf ("profile +N - perform N equilibrium iterations/stability calculations and exit\n\n");
      return 0;
    }	
  int N;
  if (control > 0) 
    N = control;
  else
    N = 1;

  // ################### 
  // Read increment file
  // ###################
  char *c = new char[100];
  double DASPCT, DELONG, DTRIANG, DQC, DQB, DBETAN;
  FILE *file = OpenFiler ("profile.inc");
  if (fscanf (file, "%s %lf", c, &DASPCT) != 2)
    {
      printf ("Error reading profile.inc (line 1)\n");
      exit (1);
    }
  if (fscanf (file, "%s %lf", c, &DELONG) != 2)
    {
      printf ("Error reading profile.inc (line 2)\n");
      exit (1);
    }
  if (fscanf (file, "%s %lf", c, &DTRIANG) != 2)
    {
      printf ("Error reading profile.inc (line 3)\n");
      exit (1);
    }
  if (fscanf (file, "%s %lf", c, &DQC) != 2)
    {
      printf ("Error reading profile.inc (line 4)\n");
      exit (1);
    }
  if (fscanf (file, "%s %lf", c, &DQB) != 2)
    {
      printf ("Error reading profile.inc (line 5)\n");
      exit (1);
    }
  if (fscanf (file, "%s %lf", c, &DBETAN) != 2)
    {
      printf ("Error reading profile.inc (line 6)\n");
      exit (1);
    }
  fclose (file);
  delete[] c;
  printf ("\n***************\nProgram PROFILE\n***************\n");
  printf ("\nIncrement Parameters (from profile.inc):\nDASPCT = %11.4e  DELONG = %11.4e  DTRIANG = %11.4e\n",
	  DASPCT, DELONG, DTRIANG);
  printf ("DQC    = %11.4e  DQB    = %11.4e  DBETAN  = %11.4e\n",
	  DQC, DQB, DBETAN);
  
  // ###############
  // Read input file
  // ###############
  int iter;
  for (int i = 0; i < N; i++)
    {	
      iter = 0;
      char  *c = new char[100];
      int    NBPTS, NPTS, NS, NT, NRBOX, NZBOX;
      double SB, ASPCT, ELONG, TRIANG, QC, QB, BETAN, MUP, MU, NUP, NU, P0;
      
      // Read input file
      FILE *file = OpenFiler ("profile.in");
      if (fscanf (file, "%s %d", c, &NBPTS) != 2)
	{
	  printf ("Error reading profile.in (line 1)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %d", c, &NPTS) != 2)
	{
	  printf ("Error reading profile.in (line 2)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %d", c, &NS) != 2)
	{
	  printf ("Error reading profile.in (line 3)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %d", c, &NT) != 2)
	{
	  printf ("Error reading profile.in (line 4)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %d", c, &NRBOX) != 2)
	{
	  printf ("Error reading profile.in (line 5)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %d", c, &NZBOX) != 2)
	{
	  printf ("Error reading profile.in (line 6)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &SB) != 2)
	{
	  printf ("Error reading profile.in (line 7)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &ASPCT) != 2)
	{
	  printf ("Error reading profile.in (line 8)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &ELONG) != 2)
	{
	  printf ("Error reading profile.in (line 9)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &TRIANG) != 2)
	{
	  printf ("Error reading profile.in (line 10)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &QC) != 2)
	{
	  printf ("Error reading profile.in (line 11)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &QB) != 2)
	{
	  printf ("Error reading profile.in (line 12)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &BETAN) != 2)
	{
	  printf ("Error reading profile.in (line 13)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &MUP) != 2)
	{
	  printf ("Error reading profile.in (line 14)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &MU) != 2)
	{
	  printf ("Error reading profile.in (line 15)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &NUP) != 2)
	{
	  printf ("Error reading profile.in (line 16)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &NU) != 2)
	{
	  printf ("Error reading profile.in (line 17)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &P0) != 2)
	{
	  printf ("Error reading profile.in (line 18)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &epsb) != 2)
	{
	  printf ("Error reading profile.in (line 19)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &epsq) != 2)
	{
	  printf ("Error reading profile.in (line 20)\n");
	  exit (1);
	}
      fclose (file);
      delete[] c;
      
      printf ("\nInput Parameters (from profile.in):\nNBPTS =  %3d         NPTS  =  %3d\n",
	      NBPTS, NPTS);
      printf ("NRBOX = %4d         NZBOX = %4d         NS    =  %3d         NT     =  %3d\n",
	      NRBOX, NZBOX, NS, NT);
      printf ("SB    = %11.4e  ASPCT = %11.4e  ELONG = %11.4e  TRIANG = %11.4e\n",
	      SB, ASPCT, ELONG, TRIANG);
      printf ("QC    = %11.4e  QB    = %11.4e  BETAN = %11.4e\n",
	      QC, QB, BETAN);
      printf ("MUP   = %11.4e  MU    = %11.4e  NUP   = %11.4e\n",
	      MUP, MU, NUP);
      printf ("NU    = %11.4e  P0    = %11.4e\n",
	      NU, P0);
      printf ("epsb  = %11.4e  epsq  = %11.4e\n\n",
	      epsb, epsq);
      
      // ###############
      // Begin main loop
      // ###############
      double psic, Ip, qc, qa, qb, betap, li, beta, betan, dum;

      // ............
      // Sanity check
      // ............
      if (P0 > 0.)
	{
	  printf ("Error P0 > 0.\n");
	  exit (1);
	}
      
      // ...............................
      // Check for zero beta calculation
      // ...............................
      int finitebeta = (BETAN < epsb) ? 0 : 1;
      if (!finitebeta)
	P0 = 0.;
      if (finitebeta && fabs (P0) < epsb)
	P0 = - epsb;

      // ...................
      // Initial calculation
      // ...................
      WriteNamelist (6, QC, NS, NT, NRBOX, NZBOX);
      WriteProfile  (NBPTS, ASPCT, ELONG, TRIANG, NPTS, SB, P0, MUP, NUP, MU, NU);
      if (control == 0)
	{
	  WriteNamelist (3, QC, NS, NT, NRBOX, NZBOX);
	  system ("./chease");
	  return 0;
	}
      if (control == -2)
	{
	  WriteNamelist (6, QC, NS, NT, NRBOX, NZBOX);
	  return 0;
	}
      if (system ("./chease >& /dev/null"))
	{
	  printf ("Error calling CHEASE\n");
	  exit (1);
	}
      if (system ("./flux 1 >& /dev/null"))
	{
	  printf ("Error calling flux\n");
	  exit (1);
	}
      
      file = OpenFiler ("Stage2/flux.out");
      fscanf (file, "%lf %lf %lf %lf %lf %lf %lf",
	      &Ip, &qc, &qa, &betap, &li, &beta, &betan);
      fclose (file);
      file = OpenFiler ("Stage2/Edge.out");
      fscanf (file, "%lf %lf %lf", &dum, &dum, &qb);
      fclose (file);
      printf ("CHEASE/FLUX: nu = %17.10e p0 = %17.10e qc = %17.10e qb = %17.10e qa = %17.10e betaN = %17.10e\n",
	      NU, P0, qc, qb, qa, betan);

      // ..........
      // Iterate P0
      // ..........
      if (fabs(BETAN-betan) > epsb && finitebeta)
	do
	  {
	    P0 = P0 * BETAN /betan;
	    
	    WriteNamelist (6, QC, NS, NT, NRBOX, NZBOX);
	    WriteProfile  (NBPTS, ASPCT, ELONG, TRIANG, NPTS, SB, P0, MUP, NUP, MU, NU);
	    
	    if (system ("./chease >& /dev/null"))
	      {
		printf ("Error calling CHEASE\n");
		exit (1);
	      }
	    if (system ("./flux 1 >& /dev/null"))
	      {
		printf ("Error calling flux\n");
		exit (1);
	      }
	    iter++;
	    
	    file = OpenFiler ("Stage2/flux.out");
	    fscanf (file, "%lf %lf %lf %lf %lf %lf %lf",
		    &Ip, &qc, &qa, &betap, &li, &beta, &betan);
	    fclose (file);
	    file = OpenFiler ("Stage2/Edge.out");
	    fscanf (file, "%lf %lf %lf", &dum, &dum, &qb);
	    fclose (file);
	    printf ("CHEASE/FLUX: nu = %17.10e p0 = %17.10e qc = %17.10e qb = %17.10e qa = %17.10e betaN = %17.10e\n",
		    NU, P0, qc, qb, qa, betan);
	  } 
	while (fabs(BETAN-betan) > epsb);	 
      
      // ..........
      // Iterate NU
      // ..........
      if (fabs(QB-qb) > epsq)
	do
	  {
	    NU = NU * QB /qb;
	    
	    WriteNamelist (6, QC, NS, NT, NRBOX, NZBOX);
	    WriteProfile  (NBPTS, ASPCT, ELONG, TRIANG, NPTS, SB, P0, MUP, NUP, MU, NU);
	    
	    if (system ("./chease >& /dev/null"))
	      {
		printf ("Error calling CHEASE\n");
		exit (1);
	      }
	    if (system ("./flux 1 >& /dev/null"))
	      {
		printf ("Error calling flux\n");
		exit (1);
	      }
	    iter++;
	    
	    file = OpenFiler ("Stage2/flux.out");
	    fscanf (file, "%lf %lf %lf %lf %lf %lf %lf",
		    &Ip, &qc, &qa, &betap, &li, &beta, &betan);
	    fclose (file);
	    file = OpenFiler ("Stage2/Edge.out");
	    fscanf (file, "%lf %lf %lf", &dum, &dum, &qb);
	    fclose (file);
	    printf ("CHEASE/FLUX: nu = %17.10e p0 = %17.10e qc = %17.10e qb = %17.10e qa = %17.10e betaN = %17.10e\n",
		    NU, P0, qc, qb, qa, betan);
	    
	    // ..........
	    // Iterate P0
	    // ..........
	    if (finitebeta)
	      do
		{
		  P0 = P0 * BETAN /betan;
		  
		  WriteNamelist (6, QC, NS, NT, NRBOX, NZBOX);
		  WriteProfile  (NBPTS, ASPCT, ELONG, TRIANG, NPTS, SB, P0, MUP, NUP, MU, NU);
		  
		  if (system ("./chease >& /dev/null"))
		    {
		      printf ("Error calling CHEASE\n");
		      exit (1);
		    }
		  if (system ("./flux 1 >& /dev/null"))
		    {
		      printf ("Error calling flux\n");
		      exit (1);
		    }
		  iter++;
		  
		  file = OpenFiler ("Stage2/flux.out");
		  fscanf (file, "%lf %lf %lf %lf %lf %lf %lf",
			  &Ip, &qc, &qa, &betap, &li, &beta, &betan);
		  fclose (file);
		  file = OpenFiler ("Stage2/Edge.out");
		  fscanf (file, "%lf %lf %lf", &dum, &dum, &qb);
		  fclose (file);
		  printf ("CHEASE/FLUX: nu = %17.10e p0 = %17.10e qc = %17.10e qb = %17.10e qa = %17.10e betaN = %17.10e\n",
			  NU, P0, qc, qb, qa, betan);
		} 
	      while (fabs(BETAN-betan) > epsb && iter < itmax);	 
	  } 
	while (fabs(QB-qb) > epsq && iter < itmax);
      
      // Output equilibrium data
      file = OpenFile ("Stage1/profile.out");
      fprintf (file, "%d %d %e %e %e %e %e %e %e %e %e %e %e %e\n",
	       NBPTS, NPTS, ASPCT, ELONG, TRIANG, QC, QB, BETAN, MUP, MU, NUP, NU, P0, SB);
      fclose (file);
      // #############
      // End main loop
      // #############
      
      // ################################
      // Calculate flux coordinate system
      // ################################
      system ("./flux 1");
      
      // #################################################
      // Perform stability calculation if control positive
      // #################################################
      if (control >  0) 
	{
	  // Calculate stability matrix
	  system ("thawc 1");
	  
	  // Calculate layer data
	  system ("ggj 1");
      
	  // Archive calculation data
	  system ("cat Stage1/profile.out  >> profile.out.1");
	  system ("cat Stage2/flux.out     >> flux.out.1");
	  system ("cat Stage2/jk.out       >> jk.out.1");
	  system ("cat Stage2/rx.out       >> rx.out.1");
	  system ("cat Stage2/tt.out       >> tt.out.1");
	  system ("cat Stage2/Rt.out       >> Rt.out.1");
	  system ("cat Stage2/Zt.out       >> Zt.out.1");
	  system ("cat Stage3/info.out     >> info.out.1");
	  system ("cat Stage3/xie.out      >> xie.out.1");
	  system ("cat Stage3/modes.out    >> modes.out.1");
	  system ("cat Stage3/ratsur.out   >> ratsur.out.1");
	  system ("cat Stage3/Delta.out    >> Delta.out.1");
	  system ("cat Stage3/rational.out >> rational.out.1");
	  system ("cat Stage4/ggj.out      >> ggj.out.1");
	}
      
      // ##################
      // Rewrite profile.in
      // ##################
      file = OpenFile ("profile.in");
      fprintf (file, "NBPTS  %3d\n",     NBPTS);
      fprintf (file, "NPTS   %3d\n",     NPTS);
      fprintf (file, "NS     %3d\n",     NS);
      fprintf (file, "NT     %3d\n",     NT);
      fprintf (file, "NRBOX  %3d\n",     NRBOX);
      fprintf (file, "NZBOX  %3d\n",     NZBOX);
      fprintf (file, "SB     %11.4e\n",  SB);
      fprintf (file, "ASPCT  %11.4e\n",  ASPCT  + DASPCT);
      fprintf (file, "ELONG  %11.4e\n",  ELONG  + DELONG);
      fprintf (file, "TRIANG %11.4e\n",  TRIANG + DTRIANG);
      fprintf (file, "QC     %11.4e\n",  QC     + DQC);
      fprintf (file, "QB     %11.4e\n",  QB     + DQB);
      fprintf (file, "BETAN  %11.4e\n",  BETAN  + DBETAN);
      fprintf (file, "MUP    %11.4e\n",  MUP);
      fprintf (file, "MU     %11.4e\n",  MU);
      fprintf (file, "NUP    %11.4e\n",  NUP);
      fprintf (file, "NU     %17.10e\n", NU);
      fprintf (file, "P0     %17.10e\n", P0);
      fprintf (file, "EPSB   %17.10e\n", epsb);
      fprintf (file, "EPQQ   %17.10e\n", epsq);
      fclose  (file);
    }
}

// #####################################
// Function to open new file for writing
// #####################################
FILE* OpenFile (char *filename)
{
  FILE *file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

// ##########################################
// Function to open existing file for reading
// ##########################################
FILE* OpenFiler (char *filename)
{
  FILE* file = fopen (filename, "r");
  if (file == NULL) 
    {
      printf ("Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

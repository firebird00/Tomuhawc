// Stage1.cpp

#include "Flux.h"

// #########################################################
// Function to input EQDSK data and output intermediate data
// #########################################################
void Flux::Stage1 (int flag)
{
  // ...............
  // Read EQDSK file
  // ...............

  char  *FILENAME = "EQDSK_COCOS_02.OUT";
  int    INRBOX, INZBOX;
  double RBOXLEN, ZBOXLEN, RBOXLFT;
  double R0, B0;
  double RAXIS, ZAXIS, PSIAXIS;
  double xPSI[N], xg[N], xP[N], xggp[N], xPp[N], xq[N];
  double psi[N][N];
  int    NBPTS;
  double RBPTS[N], ZBPTS[N];
  double IP, QC, BETAP, LI;
  double xR[N], xZ[N];
  double mu0 = 4.*M_PI*1.e-7;

  int err = EQDSK_Read (FILENAME, 
			INRBOX, INZBOX, 
			RBOXLEN, ZBOXLEN, RBOXLFT, 
			R0, B0,
			RAXIS, ZAXIS, PSIAXIS, 
			xg, xP, xggp, xPp, xq, 
			psi,
			NBPTS, RBPTS, ZBPTS,
			IP, QC, QA, BETAP, LI);

  if (err) 
    {
      printf ("\nError reading EQDSK file: err = %3d\n\n", err);
      exit (1);
    }

  double BETAN = 40.*M_PI*BETA/IP/IASPCT;
  printf ("\n************\nProgram FLUX\n************\n");
  printf ("\nReading EQDSK file: %s\n\n", 
	  FILENAME);
  printf ("INRBOX  =  %d        INZBOX  =  %d        NBPTS   =  %d\n", 
	  INRBOX, INZBOX, NBPTS);
  printf ("RBOXLEN = %11.4e ZBOXLEN = %11.4e RBOXLFT = %11.4e\n",
	  RBOXLEN, ZBOXLEN, RBOXLFT);
  printf ("R0EXP   = %11.4e B0EXP   = %11.4e\n",
	  R0, B0);
  printf ("RAXIS   = %11.4e ZAXIS   = %11.4e PSIAXIS = %11.4e\n",
	  RAXIS, ZAXIS, PSIAXIS);
  printf ("IP      = %11.4e QC      = %11.4e QA      = %11.4e\n",
	  IP, QC, QA);
  printf ("BETAP   = %11.4e LI      = %11.4e BETA    = %11.4e IASPCT = %11.4e\n",
	  BETAP, LI, BETA, IASPCT);
  printf ("BETAN   = %11.4e\n\n",
	  BETAN);

  FILE *file = fopen ("Stage2/flux.out", "w");
  fprintf (file, "%17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e\n", IP, QC, QA, BETAP, LI, BETA, BETAN);
  fclose (file);

  // ................
  // if flag = 0 exit
  // ................
  if (flag == 0) 
    return;

  // ..................
  // Output stage1 data
  // ..................
  double xpsic = PSIAXIS /(R0*R0*B0);
  for (int i = 0; i < INRBOX; i++)
    {
      xR [i]  = RBOXLFT + RBOXLEN * double (i) /double (INRBOX-1);
      xPSI[i] = xpsic * double ((INRBOX - 1 - i)) /double (INRBOX-1); 
    }	
  for (int j = 0; j < INZBOX; j++)
      xZ[j] = - ZBOXLEN/2. + ZBOXLEN * double (j) /double (INZBOX-1);

  file = OpenFile ("Stage1/Points.out");
  fprintf (file, "%d %d %d\n", INRBOX, INZBOX, NBPTS);	
  fclose (file);

  file = OpenFile ("Stage1/Box.out");
  fprintf (file, "%17.10e %17.10e %17.10e %17.10e\n", xR[0] /R0, xZ[0] /R0, xR[INRBOX-1] /R0, xZ[INZBOX-1] /R0);
  fclose (file);

  file = OpenFile ("Stage1/R.out");
  for (int i = 0; i < INRBOX; i++)
    fprintf (file, "%17.10e\n", xR[i] /R0);
  fclose (file);

  file = OpenFile ("Stage1/Z.out");
  for (int j = 0; j < INZBOX; j++)
    fprintf (file, "%17.10e\n", xZ[j] /R0);
  fclose (file);

  file = OpenFile ("Stage1/Psi.out");
  for (int i = 0; i < INRBOX; i++)
    {
      for (int j = 0; j < INZBOX; j++)
	fprintf (file, "%17.10e ", psi[i][j] /PSIAXIS);
      fprintf (file, "\n");
    }
  fclose (file);

  file = OpenFile ("Stage1/Axis.out");
  fprintf (file, "%17.10e %17.10e %17.10e\n", RAXIS /R0, ZAXIS /R0, xpsic);
  fclose (file);

  file = OpenFile ("Stage1/Profiles.out");
  for (int i = 0; i < INRBOX; i++)
    fprintf (file, "%17.10e %17.10e %17.10e %17.10e %17.10e %17.10e\n", 
	     xPSI[i], xg[i] /(R0*B0), xP[i] /(B0*B0/mu0), 
	     (R0*R0*B0) * xggp[i] /(R0*R0*B0*B0), (R0*R0*B0) * xPp[i] /(B0*B0/mu0), xq[i]);
  fclose (file);

  file = OpenFile ("Stage1/Boundary.out");
  for (int i = 0; i < NBPTS; i++)
    fprintf (file, "%17.10e %17.10e\n", RBPTS[i] /R0, ZBPTS[i] /R0);
  fclose (file);
}

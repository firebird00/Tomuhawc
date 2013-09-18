// WriteProfile.cpp

#include <stdio.h>
#include <math.h>

extern FILE* OpenFile  (char *filename);

// ###################################
// Function to write CHEASE EXPEQ file
// ###################################
void WriteProfile (int NBPTS, double ASPCT, double ELONG, double TRIANG, 
		   int NPTS, double SC, double P0, double MUP, double NUP, 
		   double MU, double NU)
{
  FILE *file = OpenFile ("EXPEQ");

  int two = 2; double Zero = 0.;
  fprintf (file, "%18.8e\n%18.8e\n%18.8e\n", ASPCT, Zero, Zero);
  fprintf (file, "%5d\n", NBPTS);
  for (int i = 0; i < NBPTS; i++)
    {
      double th = double (i) * 2.*M_PI /double (NBPTS-1);
      double R  = 1. + ASPCT * cos(th + TRIANG * sin(th));
      double Z  = ASPCT * ELONG * sin(th);
      fprintf (file, "%18.8e %18.8e\n", R, Z);
    }
  fprintf (file, "%5d\n%5d\n", NPTS, two);
  for (int i = 0; i < NPTS; i++)
    {
      double s = double (i) /double (NPTS-1);
      fprintf (file, "%18.8e\n", s);
    }
  for (int i = 0; i < NPTS; i++)
    {
      double s = double (i) /double (NPTS-1);
      double x = s /SC;
      double val;
      if (x <= 1.)
	val = P0 * pow (1. - pow (x, 2.*MUP), NUP);
      else
	val = 0.;
      fprintf (file, "%18.8e\n", val);
    }
  for (int i = 0; i < NPTS; i++)
    {
      double s = double (i) /double (NPTS-1);
      double x = s /SC;
      double val;
      if (x <= 1.)
	val = pow (1. - pow (x, 2.*MU), NU);
      else
	val = 0.;
      fprintf (file, "%18.8e\n", val);
    }
  
  fclose (file);
}

// main.cpp

// Main function for program Flux.
// See Flux.h

#include "Flux.h"

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Pointers to right-hand side functions for passing to gsl adaptive integration routines

int pRhs1 (double r, const double y[], double dydr[], void *params)
{
  Flux flux = *(Flux *) params;

  int status = flux.Rhs1 (r, y, dydr, NULL);

  return status;
}

int pRhs2 (double r, const double y[], double dydr[], void *params)
{
  Flux flux = *(Flux *) params;

  int status = flux.Rhs2 (r, y, dydr, NULL);

  return status;
}

int pRhs3 (double r, const double y[], double dydr[], void *params)
{
  Flux flux = *(Flux *) params;

  int status = flux.Rhs3 (r, y, dydr, NULL);

  return status;
}
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main (int argc, char **argv)
{
  int control = 0;
  if (argc > 1)
    control = atoi(argv[1]);
  else
    {
      printf ("\nUsage:\n");
      printf ("flux   - print this message and exit\n");
      printf ("flux 0 - read EQDSK data and exit\n");
      printf ("flux 1 - read EQDSK data, construct flux coordinate system, and exit\n\n");
      return 0;
    }

  if (control == 0)
    {
      // Read EQDSK data and exit
      Flux flux;
      flux.Stage1 (0);
    }
  else
    {
      // Input Chease EQDSK metric data and output intermediate metric data.
      // Intermediate metric data in directory Stage1/
      Flux flux;
      flux.Stage1 (1);

      // Set global parameters
      flux.SetParameters ();

      // Input intermediate metric data and output Tomuhawc metric data.
      // Tomuhawc metric data in directory Stage2/
      flux.Stage2 ();
    }
}

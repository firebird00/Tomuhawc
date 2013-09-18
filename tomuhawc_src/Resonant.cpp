// Resonant.cpp

#include "Tomuhawc.h"

// ############################################
// Function to calculate resonant surface data.
// Assumes monotonically increasing q-profile.
//
// sw = 0: calculate resonant surface data 
// sw = 1: clean up 
//
// interactive .. set if in interactive mode
//
// ############################################
void Thawc::Resonant (int sw, int interactive)
{
  if (sw == 0)
    {
      FILE *logfile;

      if (interactive)
	{
	  logfile = OpenFilea ("Stage3/logfile");
	  fprintf (logfile, "Rational surfaces:\n");
	}

      // Find number of rational surfaces
      double q0 = gsl_spline_eval (Sq, r0, Aq);
      double qb = gsl_spline_eval (Sq, rb, Aq);
      double qa = gsl_spline_eval (Sq, ra, Aq);
      vac       = int (ntor * qa) - int (ntor * q0);
      num       = int (ntor * qb) - int (ntor * q0);
      if (num == 0) return;
      
      // Allocate memory
      Mres  = new int    [vac];
      Qres  = new double [vac];
      Rres  = new double [vac];
      Sres  = new double [vac];
      EFres = new double [vac];
      HHres = new double [vac];
      ETres = new double [vac];
      FTres = new double [vac];
      HTres = new double [vac];
      FAres = new double [vac];
      FRres = new double [vac];

      // Determine radii of rational surfaces
      for (int i = 0; i < vac; i++)
	{
	  Mres[i] = int (ntor * q0) + i + 1;
	  Qres[i] = double (Mres[i]) /double (ntor);
	  qval    = Qres[i];

	  Rres[i] = RootFind ();

	  // Determine properties of rational surfaces
	  double EF, HH, ET, FT, HT, FA, FR;
	  GGJCalc (Rres[i], EF, HH, ET, FT, HT, FA, FR);
	  double q  = gsl_spline_eval (Sq,  Rres[i], Aq);
	  double qp = gsl_spline_eval (Sqp, Rres[i], Aqp);
	  Sres [i] = Rres[i] * qp /Qres[i];
	  EFres[i] = EF;
	  HHres[i] = HH;
	  ETres[i] = ET;
	  FTres[i] = FT;
	  HTres[i] = HT;
	  FAres[i] = FA; 
	  FRres[i] = FR;

	  if (interactive)
	    fprintf (logfile, "m = %4d r/b = %11.4e q = %11.4e s = %11.4e E+F = %11.4e H = %11.4e KE-G = %11.4e KF+G = %11.4e KH = %11.4e res = %11.4e\n",
		     Mres[i], Rres[i]/rb, Qres[i], Sres[i], EFres[i], HHres[i], ETres[i], FTres[i], HTres[i], fabs(q - Qres[i]));
	  printf ("m = %4d r/a = %11.4e r/b = %11.4e q = %11.4e s = %11.4e nuS - nuL = %11.4e res = %11.4e\n",
		  Mres[i], Rres[i], Rres[i]/rb, Qres[i], Sres[i], sqrt(1.-4.*EFres[i]-4.*HHres[i]), fabs(q - Qres[i]));
	}
      if (interactive) fclose (logfile);

      // Output rational surface data
      logfile = OpenFile ("Stage3/rational.out");
      fprintf (logfile, "%d", vac);
      for (int i = 0; i < vac; i++)
	{
	  double nsl   = sqrt(1. - 4.*(EFres[i] + HHres[i]));
	  double fSres = sqrt(FAres[i]) * FRres[i] * abar*abar /double(ntor) /Sres[i] /Rres[i]/Rres[i] /a/a;
	  double fAres = double(ntor) * Sres[i] /sqrt(FAres[i]);
	  fprintf (logfile, " %d %e %e %e %e %e %e %e %e %e",
		   Mres[i], Rres[i], Sres[i], 
		   EFres[i], HHres[i], ETres[i], FTres[i], HTres[i], fAres, fSres);
	}
      fprintf (logfile, "\n");
      fclose (logfile);
      
      // Output more rational surface data
      logfile = OpenFile ("Stage3/ratsur.out");
      for (int i = 0; i  < vac; i++)
	fprintf (logfile, "%e %e ", Rres[i], Qres[i]);
      fprintf (logfile, "\n");
      fclose (logfile);
    }
  else
    {
      // Deallocate memory
      delete[] Mres;  delete[] Qres;  delete[] Rres;  delete[] Sres;  
      delete[] EFres; delete[] HHres; delete[] ETres; delete[] FTres;
      delete[] HTres; delete[] FAres; delete[] FRres; 
    }
}

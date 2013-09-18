// main.cpp

// Main function for program GGJ.

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <complex.h>
#include "GGJ.h"

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Prototypes of data input functions

void read_input (char[], char[]);
void read_data  (char[], char[], int, double[]);

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main (int argc, char *argv[])
{
  
  // Calling with argument triggers noniteractive mode
  if (argc > 1)
    {
      printf ("\n***********\nProgram GGJ\n***********\n\n");

      FILE *file  = fopen ("Stage3/rational.out", "r");
      FILE *file1 = fopen ("Stage3/Delta.out",    "r");
      FILE *file2 = fopen ("Stage4/ggj.out",      "w");
      
      int num;
      fscanf  (file,  "%d", &num);
      fscanf  (file1, "%d", &num);
      fprintf (file2, "%d ", num);
      
      int    *Mres  = new int[num];
      double *Rres  = new double[num];
      double *Sres  = new double[num];
      double *EFres = new double[num];
      double *HHres = new double[num];
      double *ETres = new double[num];
      double *FTres = new double[num];
      double *HTres = new double[num];
      double *fAres = new double[num];
      double *fSres = new double[num];
      double *Delta = new double[num];
      int dum;
      for (int i = 0; i < num; i++)
	{
	  fscanf (file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &Mres [i], &Rres [i], &Sres [i],
		  &EFres[i], &HHres[i], &ETres[i], &FTres[i],
		  &HTres[i], &fAres[i], &fSres[i]);
	  fscanf (file1, "%d %lf", &dum, &Delta[i]);
	  
	  GGJ ggj (EFres[i], HHres[i], ETres[i], FTres[i], HTres[i], fAres[i], fSres[i]);
	  double Deltar, Deltai, Qi, err1, err2, err3;
	  Qi = 0.1;
	  ggj.FindCrite (-1., Qi, Deltar, Deltai);
	  double nsl = sqrt(1.-4.*EFres[i]-4.*HHres[i]);
	  double Sc  = pow(fabs(Delta[i]/Deltar), 3./nsl);
	  Qi = (Deltar > 0.) ? Qi: 0.;
	  printf ("m = %3d  Ee = %11.4e  nsl = %11.4e  Deltaec = %11.4e  Sec = %11.4e  res = %11.4e\n",
		  Mres[i], Delta[i], nsl, Deltar, Sc, fabs(Deltai));
	  fprintf (file2, "%d %e %e %e %e %e %e ", 
		   Mres[i], Delta[i], nsl, Deltar, Qi, Deltai, Sc);
	}
      printf ("\n"); fprintf (file2, "\n");
      delete[] Mres;  delete[] Rres;  delete[] Sres; 
      delete[] EFres; delete[] HHres; delete[] ETres;
      delete[] FTres; delete[] HTres; delete[] fAres; delete[] fSres;
      
      fclose (file); fclose (file1); fclose (file2);
    }
  else
    {
      // Calling without argument triggers interactive mode
      GGJ ggj;
      char in_string[5];
      char *cursor = "GGJ> ";
      double data[5];
      
      // Control loop for interactive mode
      printf ("\n");
      while (1)
	{
	  read_input (cursor, in_string);
	  
	  if ((in_string[0] == 'q') && (in_string[1] == 'u'))
	    exit (0); 
	  else if ((in_string[0] == 'f') && (in_string[1] == 'a'))
	    {
	      read_data (cursor, "fA ? ", 1, data);
	      if (!ggj.SetfA (data[0])) 
		printf ("%-sFA = %11.4e\n", cursor, data[0]);
	    }  
	  else if ((in_string[0] == 'f') && (in_string[1] == 's'))
	    {
	      read_data (cursor, "fS ? ", 1, data);
	      if (!ggj.SetfS (data[0])) 
		printf ("%-sfS = %11.4e\n", cursor, data[0]);
	    }  
	  else if ((in_string[0] == 'e') && (in_string[1] == 'f'))
	    {
	      read_data (cursor, "EF ? ", 1, data);
	      if (!ggj.SetEF (data[0])) 
		printf ("%-sEF = %11.4e\n", cursor, data[0]);
	    }      
	  else if ((in_string[0] == 'h') && (in_string[1] == 's'))
	    {
	      read_data (cursor, "H ? ", 1, data);
	      if (!ggj.SetH (data[0])) 
		printf ("%-sH = %11.4e\n", cursor, data[0]);
	    }      
	  else if ((in_string[0] == 'k') && (in_string[1] == 'f'))
	    {
	      read_data (cursor, "KFG ? ", 1, data);
	      if (!ggj.SetKFG (data[0])) 
		printf ("%-sKFG = %11.4e\n", cursor, data[0]);
	    }      
	  else if ((in_string[0] == 'k') && (in_string[1] == 'e'))
	    {
	      read_data (cursor, "KEG ? ", 1, data);
	      if (!ggj.SetKEG (data[0])) 
		printf ("%-sKEG = %11.4e\n", cursor, data[0]);
	    }      
	  else if ((in_string[0] == 'k') && (in_string[1] == 'h'))
	    {
	      read_data (cursor, "KH ? ", 1, data);
	      if (!ggj.SetKH (data[0])) 
		printf ("%-sKH = %11.4e\n", cursor, data[0]);
	    }      
	  else if ((in_string[0] == 'q') && (in_string[1] == 'r'))
	    {
	      read_data (cursor, "Qr ? ", 1, data);
	      if (!ggj.SetQr (data[0])) 
		printf ("%-sQr = %11.4e\n", cursor, data[0]);
	    } 
	  else if ((in_string[0] == 'q') && (in_string[1] == 'i'))
	    {
	      read_data (cursor, "Qi ? ", 1, data);
	      if (!ggj.SetQi (data[0])) 
		printf ("%-sQi = %11.4e\n", cursor, data[0]);
	    } 
	  else if ((in_string[0] == 'n') && (in_string[1] == 's'))
	    {
	      read_data (cursor, "N ? ", 1, data);
	      if (!ggj.SetN (int(data[0]))) 
		printf ("%-sN = %4d\n", cursor, int(data[0]));
	    }
	  else if ((in_string[0] == 's') && (in_string[1] == 't'))
	    ggj.Status ();	
	  else if ((in_string[0] == 'h') && (in_string[1] == 'e'))
	    {
	      printf ("\n");  
	      printf ("ef   ... set layer parameter, EF\n"); 
	      printf ("hs   ... set layer parameter, H\n"); 
	      printf ("kfg  ... set layer parameter, KFG\n"); 
	      printf ("kef  ... set layer parameter, KEG\n"); 
	      printf ("kh   ... set layer parameter, KH\n"); 
	      printf ("fa   ... set layer parameter, fA\n"); 
	      printf ("fs   ... set layer parameter, fS\n"); 
	      printf ("qr   ... set layer parameter, Qr\n");
	      printf ("qr   ... set layer parameter, Qi\n\n");
	      printf ("ns   ... set number of grid-points, N\n\n"); 
	      printf ("solv ... solve Glasser-Green-Johnson layer equations\n");
	      printf ("es   ... scan Deltae\n");
	      printf ("os   ... scan Deltao\n"); 
	      printf ("rs   ... scan Qr\n");
	      printf ("is   ... scan Qi\n");
	      printf ("ce   ... find critical Deltae\n\n");
	      printf ("qe   ... solve tearing parity stability problem\n");
	      printf ("qo   ... solve twisting parity stability problem\n\n");
	      printf ("help ... print list of valid commands\n");
	      printf ("stat ... print current values of input parameters\n");
	      printf ("quit ... exit program\n\n");
	    }
	  else if ((in_string[0] == 's') && (in_string[1] == 'o') || 
		   (in_string[0] == 'r') && (in_string[1] == 'u'))
	    {
	      ggj.Solv();
	    }
	  else if ((in_string[0] == 'e') && (in_string[1] == 's'))
	    {
	      ggj.ScanDeltae ();
	    }  
	  else if ((in_string[0] == 'o') && (in_string[1] == 's'))
	    {
	      ggj.ScanDeltao ();
	    } 
	  else if ((in_string[0] == 'r') && (in_string[1] == 's'))
	    {
	      ggj.ScanQr ();
	    }  
	  else if ((in_string[0] == 'i') && (in_string[1] == 's'))
	    {
	      ggj.ScanQi ();
	    }
	  else if ((in_string[0] == 'c') && (in_string[1] == 'e'))
	    {
	      double x2 = -0.1, _Qi, Deltar, Deltai;
	      ggj.FindCrite (x2, _Qi, Deltar, Deltai);
	      
	      printf ("Qi = %11.4e  Deltar = %11.4e  Deltai = %11.4e\n",
		      _Qi, Deltar, Deltai);
	    }
	  else if ((in_string[0] == 'q') && (in_string[1] == 'e'))
	    {
	      double Delta, Qr, Qi, dQ, err; int iter;
	      printf ("Delta, Qr, Qi, dQ ?? ");
	      scanf ("%lf %lf %lf %lf", &Delta, &Qr, &Qi, &dQ);
	      ggj.FindQe (Delta, Qr, Qi, dQ, err, iter);
	      printf ("\nQr = %11.4e  Qi = %11.4e  err = %11.4e  iter = %4d\n\n",
		      Qr, Qi, err, iter);
	    }
	  else if ((in_string[0] == 'q') && (in_string[1] == 'o'))
	    {
	      double Delta, Qr, Qi, dQ, err; int iter;
	      printf ("Delta, Qr, Qi, dQ ?? ");
	      scanf ("%lf %lf %lf %lf", &Delta, &Qr, &Qi, &dQ);
	      ggj.FindQo (Delta, Qr, Qi, dQ, err, iter);
	      printf ("\nQr = %11.4e  Qi = %11.4e  err = %11.4e  iter = %4d\n\n",
		      Qr, Qi, err, iter);
	    }
	}
    }
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Function to read line of input from terminal and store it in array input[] of length 5

void read_input (char cursor[], char input[])
{
  int count = 0; char c;

  printf ("%-s", cursor);
  while ((c = getchar ()) && (c != '\n'))
    {
      input[count] = c;
      if (count < 4) ++count;
    }
  input[count] = '\0';
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Function to read n data items from terminal and store them in double array x[]

void read_data (char cursor[], char message[], int n, double x[])
{
  int count = 1; char c;

  printf ("%-s%-s", cursor, message);
  for (count = 1; count <= n; count++)
    scanf ("%lf", &x[count - 1]);
  while ((c = getchar ()) && (c != '\n')) {}
}


// Flux.cpp

#include "Flux.h"

// ###########
// Constructor
// ###########
Flux::Flux ()
{
  J     = 64;      // Number of points in r grid
  K     = 128;     // Number of points in theta grid
  sc    = 0.2;     // Controls central concentration of r grid points
  tc    = 1.;      // Controls inboard concentration of th grid points
  Gamma = 1.6667;  // Ratio of specific heats
  h0    = 1.e-6;   // Initial integration step-length
  acc   = 1.e-14;  // Integration accuracy

  ntor  = 1;       // Toroidal mode number
}

// #################################
// Function to set global parameters
// #################################
void Flux::SetParameters ()
{
  char   *c = new char[100];
  int    in;
  double db;
  FILE *file = OpenFiler ("flux.in");
  if (fscanf (file, "%s %d",  c, &in) == 2)
    if (in > 2) J = in;
  if (fscanf (file, "%s %d",  c, &in) == 2)
    if (in > 2) K = in;
  if (fscanf (file, "%s %lf", c, &db) == 2)
    if (db > 0.) sc = db;
  if (fscanf (file, "%s %lf", c, &db) == 2)
    if (db > 0.) tc = db;
  if (fscanf (file, "%s %lf", c, &db) == 2)
    if (db > 0.) Gamma = db;
  if (fscanf (file, "%s %lf", c, &db) == 2)
    if (db > 0.) h0 = db;
  if (fscanf (file, "%s %lf", c, &db) == 2)
    if (db > 0.) acc = db;
  fclose (file);
  delete[] c;
  
  printf ("Input Parameters (from flux.in):\n");
  printf ("J  = %4d         K   = %4d\n", J, K);
  printf ("sc = %11.4e  tc  = %11.4e\n", sc, tc);
  printf ("h0 = %11.4e  acc = %11.4e  Gamma = %11.4e\n", h0, acc, Gamma);

  file = fopen ("thawc.in", "r");
  if (file != NULL)
    {
      if (fscanf (file, "%s %d",  c, &in) != 2)
	{
	  printf ("Error reading thawc.in (line 1)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %d",  c, &ntor) != 2)
	{
	  printf ("Error reading thawc.in (line 2)\n");
	  exit (1);
	}
    }	
  
  printf ("\nInput Parameters (from thawc.in):\n");
  printf ("ntor  = %4d\n", ntor);

  file = OpenFile ("Stage2/grid.out");
  fprintf (file, "%d %d\n", J+1, K/2-1);
  fclose (file);

  file = OpenFile ("Stage2/jk.out");
  fprintf (file, "%d %d\n", J, K);
  fclose (file);
}

// #####################################
// Function to open new file for writing
// #####################################
FILE* Flux::OpenFile (char *filename)
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
FILE* Flux::OpenFiler (char *filename)
{
  FILE* file = fopen (filename, "r");
  if (file == NULL) 
    {
      printf ("Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

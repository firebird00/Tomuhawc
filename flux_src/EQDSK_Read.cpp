// EQDSK_Read.cpp

#include "Flux.h"

// #################################################
// Function to read EQDSK file. Returns 0 on success
// #################################################
int Flux::EQDSK_Read (char *FILENAME, 
		      int &INRBOX, int &INZBOX, 
		      double &RBOXLEN, double &ZBOXLEN, double &RBOXLFT, 
		      double &R0EXP, double &B0EXP,
		      double &RAXIS, double &ZAXIS, double &PSIAXIS, 
		      double xg[], double xP[], double xggp[], double xPp[], double xq[], 
		      double psi[][N],
		      int &NBPTS, double RBPTS[], double ZBPTS[],
		      double &IP, double &QC, double &QA, double &BETAP, double &LI)
{
  FILE *file = OpenFiler (FILENAME);
  if (file == NULL) 
    return 1;
  char  *c = new char[1000];
  int    dummy, NLPTS;
  double dum;
  double RLPTS[N], ZLPTS[N];

  if (fscanf (file, "%s %s %s %s %s %s %d %d %d", 
	      c, c, c, c, c, c, &dummy, &INRBOX, &INZBOX) != 9) 
    return 2;
  if (fscanf (file, "%lf %lf %lf %lf %lf", &RBOXLEN, &ZBOXLEN, &R0EXP, &RBOXLFT, &dum) != 5) 
    return 3;
  if (fscanf (file, "%lf %lf %lf %lf %lf", &RAXIS, &ZAXIS, &PSIAXIS, &dum, &B0EXP)     != 5) 
    return 4; 
  if (fscanf (file, "%lf %lf %lf %lf %lf", &dum, &dum, &dum, &dum, &dum)               != 5)
    return 5; 
  if (fscanf (file, "%lf %lf %lf %lf %lf", &dum, &dum, &dum, &dum, &dum)               != 5) 
    return 6; 
  if (INRBOX > N)
    return 7;
  if (INZBOX > N)
    return 8;

  for (int i = 0; i < INRBOX; i++)
    if (fscanf (file, "%lf", &xg[i])   != 1)
      return 9;
  for (int i = 0; i < INRBOX; i++)
    if (fscanf (file, "%lf", &xP[i])   != 1)
      return 10;
  for (int i = 0; i < INRBOX; i++)
    if (fscanf (file, "%lf", &xggp[i]) != 1)
      return 11;
  for (int i = 0; i < INRBOX; i++)
    if (fscanf (file, "%lf", &xPp[i])  != 1)
      return 12;

  for (int j = 0; j < INZBOX; j++)
    for (int i = 0; i < INRBOX; i++)
      if (fscanf (file, "%lf", &psi[i][j]) != 1)
	return 13;

  for (int i = 0; i < INRBOX; i++)
    if (fscanf (file, "%lf", &xq[i])       != 1)
      return 14;

  if (fscanf (file, "%d %d", &NBPTS, &NLPTS)           != 2)
    return 15;
  if (NBPTS > N)
    return 16;
  if (NLPTS > N)
   return 17;
  for (int i = 0; i < NBPTS; i++)
    if (fscanf (file, "%lf %lf", &RBPTS[i], &ZBPTS[i]) != 2)
      return 18;
  for (int i = 0; i < NLPTS; i++)
    if (fscanf (file, "%lf %lf", &RLPTS[i], &ZLPTS[i]) != 2)
      return 19;

  for (int n = 0; n < 9; n++)
    if (fgets (c, 1000, file) == NULL) 
      return 20;
  if (fscanf (file, "%lf", &IP) != 1)
    return 21;
  for (int n = 0; n < 3; n++)
    if (fgets (c, 1000, file) == NULL)
      return 22;
  if (fscanf (file, "%lf", &QC) != 1)
    return 23;
  for (int n = 0; n < 1; n++)
    if (fgets (c, 1000, file) == NULL)
      return 24;
  if (fscanf (file, "%lf", &QA) != 1)
    return 25;
  for (int n = 0; n < 1; n++)
    if (fgets (c, 1000, file) == NULL)
      return 26;
  if (fscanf (file, "%lf", &BETAP) != 1)
    return 27;
  for (int n = 0; n < 2; n++)
    if (fgets (c, 1000, file) == NULL)
      return 28;
  if (fscanf (file, "%lf", &LI) != 1)
    return 29;
  for (int n = 0; n < 2; n++)
    if (fgets (c, 1000, file) == NULL)
      return 30;
  if (fscanf (file, "%lf", &BETA) != 1)
    return 31;
  for (int n = 0; n < 6; n++)
    if (fgets (c, 1000, file) == NULL)
      return 32;
  if (fscanf (file, "%lf", &IASPCT) != 1)
    return 33;

  delete[] c;
  fclose (file);
  
  return 0;
}

// Tomuhawc.cpp

#include "Tomuhawc.h"

// ###########
// Constructor
// ###########
Thawc::Thawc ()
{ 
  // Set default values of class parameters
  twist   = 0;
  ntor    = 1;
  side    = 9;
  off     = 0;
 
  r0      = 1.e-4;
  r1      = 0.2;
  nfix    = 10;

  acc     = 1.e-9;
  del     = 1.e-8;
  eps     = 1.e-4; 

  h0      = 1.e-4;
  flg     = 2;
  meth    = 3;
  skip    = 2;
  nulc    = 5.e-9;

  Eta     = 1.e-16;
  Maxiter = 30;

  // Over-ride default values from input files
  FILE *file = fopen ("thawc.in", "r");
  char *c    = new char[100];
  if (file != NULL)
    {
      if (fscanf (file, "%s %d",  c, &twist) != 2)
	{
	  printf ("Error reading thawc.in (line 1)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %d",  c, &ntor) != 2)
	{
	  printf ("Error reading thawc.in (line 2)\n");
	  exit (1);
	}	
      if (fscanf (file, "%s %d",  c, &side) != 2)
	{
	  printf ("Error reading thawc.in (line 3)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %d",  c, &off) != 2)
	{
	  printf ("Error reading thawc.in (line 3)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &r0) != 2)
	{
	  printf ("Error reading thawc.in (line 4)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &r1) != 2)
	{
	  printf ("Error reading thawc.in (line 5)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %d",  c, &nfix) != 2)
	{
	  printf ("Error reading thawc.in (line 6)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &acc) != 2)
	{
	  printf ("Error reading thawc.in (line 7)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &del) != 2)
	{
	  printf ("Error reading thawc.in (line 8)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &eps) != 2)
	{
	  printf ("Error reading thawc.in (line 9)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %lf", c, &nulc) != 2)
	{
	  printf ("Error reading thawc.in (line 10)\n");
	  exit (1);
	}
      if (fscanf (file, "%s %d", c, &skip) != 2)
	{
	  printf ("Error reading thawc.in (line 11)\n");
	  exit (1);
	}
    }
  fclose (file);
  delete[] c;

  // Print welcome message
  printf ("\n****************\nProgram TOMUHAWC\n****************\n\n");
  printf ("Input Parameters (from thawc.in):\n");
  printf ("twist = %3d  ntor = %3d  side = %3d  off = %3d  nfix = %3d\n",
	  twist, ntor, side, off, nfix);
  printf ("r0  = %11.4e  r1  = %11.4e\n", 
	  r0, r1); 
  printf ("acc = %11.4e  del = %11.4e  eps = %11.4e  nulc = %11.4e  skip = %3d\n", 
	  acc, del, eps, nulc, skip); 

  // Initialize logfile
  FILE *logfile = OpenFile ("Stage3/logfile");
  fclose (logfile);
}

// ##########
// Destructor
// ##########
Thawc::~Thawc ()
{}

// ###########################
// Status function
//
//  sw = 0 - write to logfile
//  sw = 1 - write to terminal
//
// ###########################
void Thawc::Status (int sw)
{ 
  if (sw)
    {
      printf ("\ntwist = %3d  ntor = %3d  side = %3d  off = %3d\n",
	      twist, ntor, side, off);
      printf ("\nr0/a = %11.4e  r1/b = %11.4e  nfix = %3d\n",
	      r0, r1, nfix);
      printf ("\nh0 = %11.4e  acc = %11.4e  flg = %2d  meth = %2d  skip = %2d\n",
	      h0, acc, flg, meth, skip, del, nulc, eps);
      printf ("\ndel = %11.4e  nulc = %11.4e  eps = %11.4e\n",
	      del, nulc, eps);
      printf ("\nEta = %11.4e  Maxiter = %3d\n\n",
	      Eta, Maxiter); 
    }
  else
    {
      FILE *logfile = OpenFilea ("Stage3/logfile");
      fprintf (logfile, "\nCalculation Parameters:\n");
      fprintf (logfile, "twist = %3d  ntor = %3d  side = %3d  off = %3d\n",
	       twist, ntor, side, off);
      fprintf (logfile, "r0/a = %11.4e  r1/b = %11.4e  nfix = %3d\n",
	       r0, r1, nfix);
      fprintf (logfile, "h0 = %11.4e  acc = %11.4e  flg = %2d  meth = %2d  skip = %2d\n",
	       h0, acc, flg, meth, skip);
      fprintf (logfile, "del = %11.4e  nulc = %11.4e  eps = %11.4e\n",
	       del, nulc, eps);
      fprintf (logfile, "Eta = %11.4e  Maxiter = %3d\n",
	       Eta, Maxiter); 
      fclose (logfile);
    }	
}

// #################################
// Functions to set class parameters
// #################################
int Thawc::Settwist (int _twist)
{
  if (_twist < 0 || _twist > 1)
    {
      printf ("Thawc::Settwist: Error - invalid data value\n");
      return 1;
    }
  twist = _twist;
  return 0;
}
int Thawc::Setntor (int _ntor)
{
  if (_ntor < 1)
    {
      printf ("Thawc::Setntor: Error - invalid data value\n");
      return 1;
    }
  ntor = _ntor;
  return 0;
}
int Thawc::Setr0 (double _r0)
{
  if (_r0 < 0. || _r0 > 1.)
    {
      printf ("Thawc::Setr0: Error - invalid data value\n");
      return 1;
    }
  r0 = _r0;
  return 0;
}
int Thawc::Setside (int _side)
{
  if (_side < 0)
    {
      printf ("Thawc::Setside: Error - invalid data value\n");
      return 1;
    }
  side = _side;
  return 0;
}
int Thawc::Setoff (int _off)
{
  if (_off < 0)
    {
      printf ("Thawc::Setoff: Error - invalid data value\n");
      return 1;
    }
  off = _off;
  return 0;
}
int Thawc::Setnfix (int _nfix)
{
  if (_nfix < 0)
    {
      printf ("Thawc::Setnfix: Error - invalid data value\n");
      return 1;
    }
  nfix = _nfix;
  return 0;
}
int Thawc::Seth0 (double _h0)
{
  if (_h0 < 0.)
    {
      printf ("Thawc::Seth0: Error - invalid data value\n");
      return 1;
    }
  h0 = _h0;
  return 0;
}
int Thawc::SetEta (double _Eta)
{
  if (_Eta < 0.)
    {
      printf ("Thawc::SetEta: Error - invalid data value\n");
      return 1;
    }
  Eta = _Eta;
  return 0;
}
int Thawc::Setr1 (double _r1)
{
  if (_r1 > 1.)
    {
      printf ("Thawc::Setr1: Error - invalid data value\n");
      return 1;
    }
  r1 = _r1;
  return 0;
}
int Thawc::Setacc (double _acc)
{
  if (_acc < 0.)
    {
      printf ("Thawc::Setacc: Error - invalid data value\n");
      return 1;
    }
  acc = _acc;
  return 0;
}
int Thawc::Setflg (int _flg)
{
  if (_flg < 0 || _flg > 2)
    {
      printf ("Thawc::Setflg: Error - invalid data value\n");
      return 1;
    }
  flg = _flg;
  if (flg == 0)
    printf ("Absolute truncation error\n");
  else if (flg == 1)
    printf ("Relative truncation error\n");
  else if (flg == 2)
    printf ("Mixed truncation error\n");
  return 0;
}
int Thawc::Setmeth (int _meth)
{
  if (_meth < 0 || _meth > 6)
    {
      printf ("Thawc::Setmeth: Error - invalid data value\n");
      return 1;
    }
  meth = _meth;
  if (meth == 0)
    printf ("4th order (classical) Runge-Kutta method\n");
  else if (meth == 1)
    printf ("Embedded Runge-Kutta-Fehlberg (4,5) method\n");
  else if (meth == 2)
    printf ("Embedded Runge-Kutta Cash-Karp (4,5) method\n");
  else if (meth == 3)
    printf ("Embedded Runge-Kutta Prince-Dormand (8,9) method\n");
  else if (meth == 4)
    printf ("Implicit 4th order Runge-Kutta at Gaussian points\n");
  else if (meth == 5)
    printf ("M=2 implicit Gear method\n");
  else if (meth == 6)
    printf ("Variable order Adams/BDF method (lsoda)\n");
  return 0;
}
int Thawc::Setskip (int _skip)
{
  if (_skip < 1)
    {
      printf ("Thawc::Setskip: Error - invalid data value\n");
      return 1;
    }
  skip = _skip;
  return 0;
}
int Thawc::Setdel (double _del)
{
  if (_del < 0.)
    {
      printf ("Thawc::Setdel: Error - invalid data value\n");
      return 1;
    }
  del = _del;
  return 0;
}
int Thawc::Setnulc (double _nulc)
{
  if (_nulc < 0.)
    {
      printf ("Thawc::Setnulc: Error - invalid data value\n");
      return 1;
    }
  nulc = _nulc;
  return 0;
}
int Thawc::Seteps (double _eps)
{
  if (_eps < 0.)
    {
      printf ("Thawc::Seteps: Error - invalid data value\n");
      return 1;
    }
  eps = _eps;
  return 0;
}




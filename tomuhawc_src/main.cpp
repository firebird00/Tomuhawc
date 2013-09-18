// main.cpp

// Main function for program Tomuhawc.
// See Tomuhawc.h for program description.

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "Tomuhawc.h"

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Pointers to right-hand side functions for passing to gsl adaptive integration routines

int pRhs (double r, const double y[], double dydr[], void *params)
{
  Thawc thawc = *(Thawc *) params;

  int status = thawc.Rhs (r, y, dydr, NULL);
  
  func++;

  return status;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Prototypes of data input functions

void read_input (char[], char[]);
void read_data  (char[], char[], int, double[]);

// %%%%%%%%%%%%%%%%%
// Global parameters
int step, func;

// %%%%%%%%%%%%%
// Main function
int main (int argc, char *argv[])
{
  // Calling with argument triggers noniteractive mode
  if (argc > 1)
    {
      Thawc thawc;
      int err = thawc.Solv (0);
      exit (0);
    }

  // Calling without argument triggers interactive mode
  Thawc thawc;
  char in_string[5];
  char *cursor = "Tomuhawc> ";
  double data[5];
  
  // Main control loop
  printf ("\n");
  while (1)
    {
      read_input (cursor, in_string);

      if ((in_string[0] == 'q') && (in_string[1] == 'u'))
	exit (0); 
      else if ((in_string[0] == 'h') && (in_string[1] == '0'))
	{
	  read_data (cursor, "h0 ? ", 1, data);
	  if (!thawc.Seth0 (data[0])) 
	    printf ("%-sh0 = %11.4e\n", cursor, data[0]);
	} 
      else if ((in_string[0] == 't') && (in_string[1] == 'w'))
	{
	  read_data (cursor, "twist ? ", 1, data);
	  if (!thawc.Settwist (int(data[0]))) 
	    printf ("%-stwist = %4d\n", cursor, int(data[0]));
	} 
      else if ((in_string[0] == 'n') && (in_string[1] == 't'))
	{
	  read_data (cursor, "ntor ? ", 1, data);
	  if (!thawc.Setntor (int(data[0]))) 
	    printf ("%-sntor = %4d\n", cursor, int(data[0]));
	}    
      else if ((in_string[0] == 's') && (in_string[1] == 'i'))
	{
	  read_data (cursor, "side ? ", 1, data);
	  if (!thawc.Setside (int(data[0]))) 
	    printf ("%-sside = %4d\n", cursor, int(data[0]));
	}
      else if ((in_string[0] == 'o') && (in_string[1] == 'f'))
	{
	  read_data (cursor, "off ? ", 1, data);
	  if (!thawc.Setoff (int(data[0]))) 
	    printf ("%-soff = %4d\n", cursor, int(data[0]));
	}
      else if ((in_string[0] == 'n') && (in_string[1] == 'f'))
	{
	  read_data (cursor, "nfix ? ", 1, data);
	  if (!thawc.Setnfix (int(data[0]))) 
	    printf ("%-snfix = %4d\n", cursor, int(data[0]));
	}   
      else if ((in_string[0] == 'e') && (in_string[1] == 't'))
	{
	  read_data (cursor, "Eta ? ", 1, data);
	  if (!thawc.SetEta (data[0])) 
	    printf ("%-sEta = %11.4e\n", cursor, data[0]);
	} 
      else if ((in_string[0] == 'h') && (in_string[1] == '0'))
	{
	  read_data (cursor, "h0 ? ", 1, data);
	  if (!thawc.Seth0 (data[0])) 
	    printf ("%-sh0 = %11.4e\n", cursor, data[0]);
	} 
      else if ((in_string[0] == 'r') && (in_string[1] == '1'))
	{
	  read_data (cursor, "r1 ? ", 1, data);
	  if (!thawc.Setr1 (data[0])) 
	    printf ("%-sr1 = %11.4e\n", cursor, data[0]);
	}  
      else if ((in_string[0] == 'r') && (in_string[1] == '0'))
	{
	  read_data (cursor, "r0 ? ", 1, data);
	  if (!thawc.Setr0 (data[0])) 
	    printf ("%-sr0 = %11.4e\n", cursor, data[0]);
	} 
      else if ((in_string[0] == 'a') && (in_string[1] == 'c'))
	{
	  read_data (cursor, "acc ? ", 1, data);
	  if (!thawc.Setacc (data[0])) 
	    printf ("%-sacc = %11.4e\n", cursor, data[0]);
	}    
      else if ((in_string[0] == 'f') && (in_string[1] == 'l'))
	{
	  read_data (cursor, "flg ? ", 1, data);
	  if (!thawc.Setflg (int(data[0]))) 
	    printf ("%-sflg = %4d\n", cursor, int(data[0]));
	}    
      else if ((in_string[0] == 'm') && (in_string[1] == 'e'))
	{
	  read_data (cursor, "meth ? ", 1, data);
	  if (!thawc.Setmeth (int(data[0]))) 
	    printf ("%-smeth = %4d\n", cursor, int(data[0]));
	}   
      else if ((in_string[0] == 's') && (in_string[1] == 'k'))
	{
	  read_data (cursor, "skip ? ", 1, data);
	  if (!thawc.Setskip (int(data[0]))) 
	    printf ("%-sskip = %4d\n", cursor, int(data[0]));
	}   
      else if ((in_string[0] == 'd') && (in_string[1] == 'e'))
	{
	  read_data (cursor, "del ? ", 1, data);
	  if (!thawc.Setdel (data[0])) 
	    printf ("%-sdel = %11.4e\n", cursor, data[0]);
	}  
      else if ((in_string[0] == 'e') && (in_string[1] == 'p'))
	{
	  read_data (cursor, "eps ? ", 1, data);
	  if (!thawc.Seteps (data[0])) 
	    printf ("%-seps = %11.4e\n", cursor, data[0]);
	}     
      else if ((in_string[0] == 'n') && (in_string[1] == 'u'))
	{
	  read_data (cursor, "nulc ? ", 1, data);
	  if (!thawc.Setnulc (data[0])) 
	    printf ("%-snulc = %11.4e\n", cursor, data[0]);
	} 
      else if ((in_string[0] == 's') && (in_string[1] == 't'))
	thawc.Status (1);	
      else if ((in_string[0] == 'h') && (in_string[1] == 'e'))
	{
	  printf ("\n");  
	  printf ("twis ... set mode parameter, twist\n");
	  printf ("ntor ... set mode parameter, ntor\n");
	  printf ("side ... set mode parameter, side\n");
	  printf ("off  ... set mode parameter, off\n\n");
	  printf ("r1   ... set fixup parameter, r1\n");
	  printf ("nfix ... set fixup parameter, nfix\n\n");
	  printf ("r0   ... set intergration parameter, r0\n");
	  printf ("h0   ... set integration parameter, h0\n");
	  printf ("acc  ... set integration parameter, acc\n");
	  printf ("flg  ... set integration parameter, flg\n");
	  printf ("meth ... set integration parameter, meth\n");
	  printf ("skip ... set integration parameter, skip\n\n");
	  printf ("del  ... set jump parameter, del\n");
	  printf ("nucl ... set jump parameter, nulc\n");
	  printf ("eps  ... set finite difference parameter, eps\n\n");
	  printf ("eta  ... set zero finding parameter, Eta\n\n");
	  printf ("solv ... perform stability calculation\n\n");
	  printf ("help ... print list of valid commands\n");
	  printf ("stat ... print current values of input parameters\n");
	  printf ("quit ... exit program\n\n");
	}
      else if ((in_string[0] == 's') && (in_string[1] == 'o'))
	{
	  int err = thawc.Solv (1);
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

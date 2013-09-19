// WriteNamelist.cpp

// NSURF    = 6, NEQDSK = 0,             :- plasma boundary data read from EXPEQ
// NPPFUN   = 4, NEQDSK = 0              :- pressure profile data read from EXPEQ
// NFUNC    = 4, NEQDSK = 0              :- current profile data read from EXPEQ
// RC       = 1.                         :- centroid of plasma boundary at R = 1.
// NSTTP    = 1                          :- current profile is TT'
// NTFM0    = 0                          :- T_edge = 1.
// NRSCAL   = 1                          :- equilibrium rescaled such that magnetic axis at R=1
// NCSCAL   = 1, CSSPEC = 0., QSPEC = QC :- equilibrium scaled such that Q=QC at magnetic axis
// NBAL     = 0, NOPT = 0                :- Disable balloning optimization      
// NVERBOSE = 0                          :- Minimal output
// NS                                    :- number of radial grid points
// NT                                    :- number of poloidal grid points
// NIDEAL   = 6                          :- EQDSK output (for Flux/Tomuhawc)
// NRBOX, NZBOX                          :- EQDSK resolution
// R0EXP    = 10. (m), B0EXP = 2.5 (T)   :- EQDSK scaling parameters
// NIDEAL   = 3                          :- INP1 output (for PEST)
// NEGP     = 0, NER = 2                 :- Choose correct PEST Jacobian
// NPSI                                  :- number of radial grid points in PEST output
// NTNOVA                                :- number of poloidal grid points in PEST output

#include <stdio.h>
#include <math.h>

extern FILE* OpenFile  (char *filename);

// ######################################
// Function to write CHEASE namelist file
// ######################################
void WriteNamelist (int NIDEAL, double QC, int NS, int NT, int NRBOX, int NZBOX)
{
  FILE *file = OpenFile ("chease_namelist");

  fprintf (file, "***\n*** Output from Profile\n***\n***\n");
  fprintf (file, "&EQDATA\nNEQDSK=0,\nNSURF = 6,\nNPPFUN = 4,\nNFUNC = 4,\nNSTTP = 1,\nRC = 1.,\n");
  fprintf (file, "NTMF0 = 0,\nNRSCAL = 1,\nNEGP = 0,\nNER = 2,\nNBAL = 0,\nNOPT = 0,\n");
  fprintf (file, "NVERBOSE = 0,\nNS = %d,\nNT = %d,\nNPSI = %d,\nNTNOVA = %d,\n", 
	   NS, NT, NS, NT);	
  fprintf (file, "NIDEAL = %d,\nNRBOX = %d,\nNZBOX = %d,\nR0EXP = 10.,\nB0EXP = 2.5,\n", 
	   NIDEAL, NRBOX, NZBOX);
  fprintf (file, "NCSCAL = 1,\nCSSPEC = 0.,\nQSPEC = %f,\n&END\n",
	   QC);

  fclose (file);
}

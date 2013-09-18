// Segment_Fixup.cpp

#include "Tomuhawc.h"

// ####################################################################################
// Function to integrate multiple solution vectors to r = rx
// while performing nf fixups
//
// r                               .. flux surface label
// rx                              .. integrate until r = rx
// nf                              .. number of fixups
// YY  (dim1, dim+vac)             .. solution vector at radius r
//  YY(i=0,  dim -1;k=0,dim-1)        ... psi: i     .. poloidal harmonic index
//  YY(i=dim,dim1-1;k=0,dim-1)        ... Z  : i-dim .. poloidal harmonic index
//                                             k     .. index of independent solutions
//                                                       launched from axis
//  YY(i=0,  dim -1;j=dim,dim+vac-1)  ... psi: i     .. poloidal harmonic index
//  YY(i=dim,dim1-1;j=dim,dim+vac-1)  ... Z  : i-dim .. poloidal harmonic index
//                                             j     .. index of small solutions
//                                                       launched from rational surfaces 
// Psi(vac, dim+vac)               .. coefficient of large solution
//  Psi(i, j)                         ... coefficient at ith rational surface due to jth solution
// dPsi(vac, dim+vac)              .. coefficient of small solution
//  Psi(i, j)                         ... coefficient at ith rational surface due to jth solution
// logf                            .. fixup spacing flag
// flag                            .. truncation error flag
// edge                            .. set if rx is wall radius
// interactive                     .. set if in iteractive mode
// _adap                           .. name of adpative integration monitor file
// _soln                           .. name of solution file
// _edge                           .. name of edge data file
//
// side                            .. number of sideband harmonics
// vac                             .. total number of rational surfaces (including vacuum surfaces)
// num                             .. number of rational surfaces in plasma
// dim  = vac + 2*side             .. dimension of coupling matrices
// dim1 = 2*dim                    .. dimension of single solution vector 
//
// ####################################################################################
void Thawc::Segment_Fixup (double &r, double rx, int nf, 
			   gsl_matrix *YY, gsl_matrix *Psi, gsl_matrix *dPsi, 
			   int logf, int flag, int edge, int interactive, 
			   char *_adap, char *_soln, char *_edge)
{
  double rw = r, RE;
  int lmax = (nf == 0) ? 1 : nf;
  for (int l = 1; l <= lmax; l++)
    {
      // Determine next fixup radius
      if (nf < 2)
	RE = rx;
      else
	{
	  if (logf)
	    RE = exp(log(rw) + (log(rx) - log(rw)) *double(l) /double(nf));
	  else
	    RE = rw + (rx - rw) *double(l) /double(nf);
	}
      
      // Integrate solution vectors to next fixup radius
      printf ("Integrate solution vectors to fixup radius: r/b = %11.4e\n", RE/rb);
      Segment (r, RE, YY, flag, edge, interactive, _adap, _soln, _edge);

      // Fixup solution vector
      if (interactive) LogVector (r, YY);
      if (nf > 0)
	{
	  Fixup (YY, Psi, dPsi, edge, interactive, _soln, _edge);
	  if (interactive) LogVector (r, YY);
	}
    }
}

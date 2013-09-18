// Jump1.cpp

#include "Tomuhawc.h"

// #####################################################################################################
// Function to evolve single solution vector across jth rational surface
// while calculating reconnected flux
// 
// delta               .. distance of closest approach to rational surface 
//                         (rational surface at r+delta)
// i                   .. index of independent solution vector
//                         i = 0,   dim-1      .. independent solution vectors launched from axis
//                         i = dim, dim+vac-1  .. small solution vectors launched from rational surfaces
// j                   .. index of rational surface to jump across (j = 0, vac-1)
// Y[dim1]             .. solution vector at radius r
//  Y(k=0,  dim -1)         psi: k     .. poloidal harmonic index
//  Y(k=dim,dim1-1)         Z  : k-dim .. poloidal harmonic index
// Psi                 .. coefficient of large solution at rational surface
// dPsi                .. coefficient of small solution at rational surface
// _tear               .. 1/0 for tearing/twisting parity
// interactive         .. set if in interactive mode
//
// side                .. number of sideband harmonics
// vac                 .. total number of rational surfaces (including vacuum surfaces)
// num                 .. number of rational surfaces in plasma
// dim  = vac + 2*side .. dimension of coupling matrices
// dim1 = 2*dim        .. dimension of single solution vector 
//
// #####################################################################################################
void Thawc::Jump1 (double delta, int i, int j, double Y[], double& Psi, double& dPsi, int _tear, int interactive)
{
  // Iteration switch for finite pressure case
  int iter = 1;

  FILE *logfile;

  if (interactive)
    logfile = OpenFilea ("Stage3/logfile");

  // ...............................
  // Find index of resonant harmonic
  // ...............................
  int km;
  for (int k = 0; k < dim; k++)
    if (mpol[k] == Mres[j]) 
      km = k;

  // .....................................
  // Calculate m, r, q at rational surface
  // .....................................
  double mm = double (Mres[j]);
  double rm = Rres[j];
  double qm = Qres[j];

  // ............................
  // Calculate L1, P1, T1, and sm
  // ............................
  double L1, P1, T1, sm;
  L1P1T1 (rm, mm, km, L1, P1, T1, sm);

  // ...........................
  // Calculate Mercier indicies
  // ..........................
  Couple (rm);
 
  double L0   = - gsl_matrix_get (Lmat, km, km) /mm/sm;
  double P0   = - gsl_matrix_get (Pmat, km, km) /mm/sm;
  double RT   = sqrt (0.25 + L0 *P0);

  double nuL  = 0.5 - RT;
  double nuS  = 0.5 + RT;

  double dnuL = pow (delta, nuL);
  double dnuS = pow (delta, nuS);
 
  // ............
  // Calculate CC
  // ............
  double CC = 0.;
  for (int k = 0; k < dim; k++)
    {
      if (k != km)
	{
	  double kk = double (mpol[k] - mpol[km]);
	  CC += (gsl_matrix_get (Pmat, km, k) *gsl_matrix_get (Lmat, km, k) 
		 - gsl_matrix_get (Mmat, km, k) *gsl_matrix_get (Nmat, km, k)) /kk;
	}
    }
  CC /= - mm *sm *rm;
  
  double  AL, AS, bS;
  double *akt = new double[dim];
  double *bkt = new double[dim];

  // ++++++++++++++++++++
  // Finite pressure case
  // ++++++++++++++++++++
  if (fabs(nuL) > nulc)
    {
      // Calculate bL, bS
      double bL = nuL/L0;
             bS = nuS/L0;

      // Calculate lamL, gamL
      double lamL = (L0 *P1 /nuL + T1 + nuL *(L1 /L0 - 2.)) /2./rm + CC /nuL;
      double gamL = ((1. + nuL) *(P1 /nuL + T1 /L0 - nuL /L0) + P0 *(L1 /L0 - 1.)) /2./rm + CC /L0;

      // Find the psi and Z
      double *psi = new double[dim];
      double *Z   = new double[dim];
      for (int k = 0; k < dim; k++)
	{
	  psi[k] = Y[k];
	  Z  [k] = Y[dim+k]; 
	}
      
      // Calculate the ak, bk, akt, bkt
      double *ak  = new double[dim];
      double *bk  = new double[dim];
      for (int k = 0; k < dim; k++)
	{
	  if (k != km)
	    {
	      ak [k] = - (gsl_matrix_get (Mmat, k, km) /nuL
			  + gsl_matrix_get (Lmat, k, km) /L0) /mm/sm;
	      bk [k] = - (gsl_matrix_get (Nmat, k, km) /L0
			  + gsl_matrix_get (Pmat, k, km) /nuL) /mm/sm;
	      akt[k] = - (gsl_matrix_get (Mmat, k, km) /nuS
			  + gsl_matrix_get (Lmat, k, km) /L0) /mm/sm;
	      bkt[k] = - (gsl_matrix_get (Nmat, k, km) /L0
			  + gsl_matrix_get (Pmat, k, km) /nuS) /mm/sm;
	    }
	  else
	    {
	      ak [k] = 0.;
	      bk [k] = 0.;
	      akt[k] = 0.;
	      bkt[k] = 0.;
	    }
	}
 
      // Initialize AL
      double AC, BC;
      double *Pbar = new double[dim];
      double *Zbar = new double[dim];
     
      if (iter)
	{
	  // Iterate calculation of AL and AS
	  AL = psi[km] /dnuL;
	  for (int iii = 0; iii < 10; iii++)
	    {
	      // Calculate the Pbar and Zbar
	      for (int k = 0; k < dim; k++)
		{
		  if (k != km)
		    {
		      Pbar[k] = psi[k] - ak[k] *AL *dnuL;
		      Zbar[k] = Z  [k] - bk[k] *AL *dnuL;
		    }
		  else
		    {
		      Pbar[k] = 0.;
		      Zbar[k] = 0.;
		    }
		}
	      
	      // Calculate sums
	      AC = 0.; BC = 0.;
	      for (int k = 0; k < dim; k++)
		{
		  if (k != km)
		    {
		      double kk = double (mpol[k] - mpol[km]);
		      AC += (gsl_matrix_get (Nmat, km, k) *Zbar[k] 
			     + gsl_matrix_get (Pmat, km, k) *Pbar[k]) /kk;
		      BC += (gsl_matrix_get (Lmat, km, k) *Zbar[k] 
			     + gsl_matrix_get (Mmat, km, k) *Pbar[k]) /kk;
		    }
		}
	      AC /= - rm*P0;
	      BC /= - rm*L0;
	      BC +=   AC/L0;
	      
	      AS = - ((Z[km] - bL*psi[km]) /dnuS + dnuL *(BC - bL *AC) 
		      + dnuL *(gamL - bL *lamL) *AL *dnuL) /(bS - bL);
	      AL = (psi[km] + AS *dnuS + AC *delta) 
		/dnuL/(1. - delta *lamL);
	    }
	}
      else
	{
	  // Calculate the Pbar and Zbar
	  for (int k = 0; k < dim; k++)
	    {
	      Pbar[k] = psi[k] - ak[k] *psi[km];
	      Zbar[k] = Z  [k] - bk[k] *psi[km];
	    }
	  
	  // Calculate sums
	  AC = 0.; BC = 0.;
	  for (int k = 0; k < dim; k++)
	    {
	      if (k != km)
		{
		  double kk = double (mpol[k] - mpol[km]);
		  AC += (gsl_matrix_get (Nmat, km, k) *Zbar[k] 
			 + gsl_matrix_get (Pmat, km, k) *Pbar[k]) /kk;
		  BC += (gsl_matrix_get (Lmat, km, k) *Zbar[k] 
			 + gsl_matrix_get (Mmat, km, k) *Pbar[k]) /kk;
		}
	    }
	  AC /= - rm*P0;
	  BC /= - rm*L0;
	  BC +=   AC/L0;
	  
	  // Calculate AL and AS
	  AS = - ((Z[km] - bL*psi[km]) + delta *(BC - bL *AC) 
		  + delta *(gamL - bL *lamL) *psi[km]) /dnuS /(bS - bL);
	  AL = (psi[km] + AS *dnuS + AC *delta) 
	    /dnuL /(1. - delta *lamL);
	}

      // Calculate L1k, M1k, N1k, P1k
      double *L1k = new double[dim]; double *M1k = new double[dim];
      double *N1k = new double[dim]; double *P1k = new double[dim];
      L1kP1k (rm, mm, km, L1k, M1k, N1k, P1k);

      // Calculate ck, dk, Pbar1, Zbar1
      double *ck    = new double[dim]; double *dk    = new double[dim];
      double *Pbar1 = new double[dim]; double *Zbar1 = new double[dim];
      for (int k = 0; k < dim; k++)
	{
	  ck[k] = - nuL *ak[k]
	    + L1k[k] *bL + M1k[k]
	    - (rm /mm/sm) *(gsl_matrix_get (Lmat, k, km) *gamL + gsl_matrix_get (Mmat, k, km) *lamL);
	  dk[k] = - ((mm*sm /(double (mpol[k]) - mm)) + nuL) *bk[k]
	    + N1k[k] *bL + P1k[k]
	    - (rm /mm/sm) *(gsl_matrix_get (Nmat, k, km) *gamL + gsl_matrix_get (Pmat, k, km) *lamL);
	  Pbar1[k] = - (rm /mm/sm) *(gsl_matrix_get (Lmat, k, km) *BC + gsl_matrix_get (Mmat, k, km) *AC);
	  Zbar1[k] = - (mm *sm /(double (mpol[k]) - mm)) *Zbar[k]
	    - (rm/mm/sm) *(gsl_matrix_get (Nmat, k, km) *BC + gsl_matrix_get (Pmat, k, km) *AC);
	  for (int kk = 0; kk < dim; kk++)
	    {
	      if (kk != km)
		{
		  double kx = double (mpol[kk] - mpol[km]);
		  ck[k]    += (gsl_matrix_get (Lmat, k, kk)   *bk[kk]   
			       + gsl_matrix_get (Mmat, k, kk) *ak[kk])   /kx;
		  dk[k]    += (gsl_matrix_get (Nmat, k, kk) *bk[kk]   
			       + gsl_matrix_get (Pmat, k, kk) *ak[kk])   /kx;
		  Pbar1[k] += (gsl_matrix_get (Lmat, k, kk) *Zbar[kk] 
			       + gsl_matrix_get (Mmat, k, kk) *Pbar[kk]) /kx;
		  Zbar1[k] += (gsl_matrix_get (Nmat, k, kk) *Zbar[kk] 
			       + gsl_matrix_get (Pmat, k, kk) *Pbar[kk]) /kx;
		}
	    }
	  ck[k]    /= (rm*(1.+nuL));
	  dk[k]    /= (rm*(1.+nuL));
	  Pbar1[k] /= rm;
	  Zbar1[k] /= rm;
	}

      // Perform jump
      for (int k = 0; k < dim; k++)
	{
	  if (_tear == 1 || j > num-1)
	    {
	      // Tearing parity: vacuum tearing parity
	      if (k == km)
		{
		  if (iter)
		    {
		      Y[k]     += 2. *delta *(AL *dnuL *lamL + AC) + 2.     *AS *dnuS;
		      Y[dim+k] += 2. *delta *(AL *dnuL *gamL + BC) + 2. *bS *AS *dnuS;
		    }
		  else
		    {
		      Y[k]     += 2. *delta *(psi[km] *lamL + AC) + 2.     *AS *dnuS;
		      Y[dim+k] += 2. *delta *(psi[km] *gamL + BC) + 2. *bS *AS *dnuS;
		    }
		}
	      else
		{
		  if (iter)
		    {
		      Y[k]     += 2. *delta *(AL *dnuL *ck[k] + Pbar1[k]) + 2. *akt[k] *AS *dnuS;
		      Y[dim+k] += 2. *delta *(AL *dnuL *dk[k] + Zbar1[k]) + 2. *bkt[k] *AS *dnuS;
		    }
		  else
		    {
		      Y[k]     += 0.;
		      Y[dim+k] += 0.;
		    }
		} 
	    }
	  else
	    {
	      // Twisting parity (excluding vacuum)
	      if (k == km)
		{
		  if (iter)
		    {
		      Y[k]     += - 2.     *AL *dnuL + 2. *delta *AC;
		      Y[dim+k] += - 2. *bL *AL *dnuL + 2. *delta *BC;
		    }
		  else
		    {
		      Y[k]     += - 2.     *psi[km] + 2. *delta *AC;
		      Y[dim+k] += - 2. *bL *psi[km] + 2. *delta *BC;
		    }	
		}
	      else
		{
		  if (iter)
		    {
		      Y[k]     += - 2. *AL *dnuL *ak[k] + 2. *delta *Pbar1[k];
		      Y[dim+k] += - 2. *AL *dnuL *bk[k] + 2. *delta *Zbar1[k];
		    }
		  else
		    {
		      Y[k]     += - 2. *psi[km] *ak[k];
		      Y[dim+k] += - 2. *psi[km] *bk[k];
		    }	
		}
	    }
	}
       
      // Log diagnostic data
      AC = (fabs(AC) < Eta) ? 0. : AC;
      BC = (fabs(BC) < Eta) ? 0. : BC;
      AL = (fabs(AL) < Eta) ? 0. : AL;
      AS = (fabs(AS) < Eta) ? 0. : AS;
      int mcen = (i < dim) ? mpol[i] : Mres[i-dim];
      if (interactive)
	{
	  if (i == 0)
	    {
	      fprintf (logfile, "nuL = %+10.3e nuS = %+10.3e L1 = %+10.3e P1 = %+10.3e T1 = %+10.3e\n",
		       nuL, nuS, L1, P1, T1);
	      fprintf (logfile, "bS = %+10.3e bL = %+10.3e gamL = %+10.3e lamL = %+10.3e dnuL = %+10.3e dnuS = %10.3e delta = %+10.3e\n",
		       bS, bL, gamL, lamL, dnuL, dnuS, delta);
	    }	
	  double rat = (fabs(AL) > 1.e-6) ? - AS /AL /(P1*L0/2./rm/nuL) : 0.;
	  fprintf (logfile, "m_cen = %+3d AC = %+10.3e BC = %+10.3e AL = %+11.4e AS = %+11.4e rat = %+11.4e\n",
		   mcen, AC, BC, AL, AS, rat);
	  for (int k = 0; k < dim; k++)
	    {
	      int mm = (k < dim) ? mpol[k] : Mres[k-dim];
	      fprintf (logfile, "m = %+3d psi[k] = %+10.3e Z[k] = %+10.3e Pbar[k] = %+10.3e Zbar[k] = %+10.3e ak[k] = %+10.3e bk[k] = %+10.3e\n",
		       mm, psi[k], Z[k], Pbar[k], Zbar[k], ak[k], bk[k]);
	    }
	}
      
      // Clean up
      delete[] psi;  delete[] Z;     
      delete[] ak;   delete[] bk;   delete[] ck;    delete[] dk;   
      delete[] L1k;  delete[] M1k;  delete[] N1k;   delete[] P1k;
      delete[] Pbar; delete[] Zbar; delete[] Pbar1; delete[] Zbar1;
    }

  // ++++++++++++++++++++
  // Zero pressure case
  // ++++++++++++++++++++
  else
    {     
      // Flag use of zero pressure jump condition in plasma
      if (j < num) 
	{
	  nflag = 1;
	  if (i == 0) 
	    printf ("WARNING: Zero pressure jump condition\n");
	}

      // Calculate bL, bS
      double bL = nuL /L0;
             bS = 1.  /L0;

      // Calculate lamL and gamL
      double lamL = L0 *P1 /rm + CC;
      double gamL =     P1 /rm;

      // Find the psi and Z
      double *psi = new double[dim];
      double *Z   = new double[dim];
      for (int k = 0; k < dim; k++)
	{
	  psi[k] = Y[k];
	  Z  [k] = Y[dim+k]; 
	}
      
      // Calculate the akh, bkh, akt, bkt
      double *akh = new double[dim];
      double *bkh = new double[dim];
      for (int k = 0; k < dim; k++)
	{
	  if (k != km)
	    {
	      akh[k] = - gsl_matrix_get (Mmat, k, km) /mm/sm;
	      bkh[k] = - gsl_matrix_get (Pmat, k, km) /mm/sm;
	      akt[k] = - (gsl_matrix_get (Mmat, k, km) /nuS
			  + gsl_matrix_get (Lmat, k, km) /L0) /mm/sm;
	      bkt[k] = - (gsl_matrix_get (Nmat, k, km) /L0
			  + gsl_matrix_get (Pmat, k, km) /nuS) /mm/sm;
	    }
	  else
	    {
	      akh[k] = 0.;
	      bkh[k] = 0.;
	      akt[k] = 0.;
	      bkt[k] = 0.;
	    }
	}
        
      // Calculate the Pbar and Zbar
      double *Pbar = new double[dim];
      double *Zbar = new double[dim];
      for (int k = 0; k < dim; k++)
	{
	  Pbar[k] = psi[k] - akh[k] *log(delta) *psi[km];
	  Zbar[k] = Z  [k] - bkh[k] *log(delta) *psi[km];
	}
      
      // Calculate sums
      double AC = 0., AD = 0., BD;
      for (int k = 0; k < dim; k++)
	{
	  if (k != km)
	    {
	      double kk = double (mpol[k] - mpol[km]);
	      AC += (gsl_matrix_get (Lmat, km, k) *Zbar[k] 
		     + gsl_matrix_get (Mmat, km, k) *Pbar[k]) /kk;
	      AD += (gsl_matrix_get (Nmat, km, k) *Zbar[k] 
		     + gsl_matrix_get (Pmat, km, k) *Pbar[k]) /kk;
	    }
	}
      AC /= rm;
      AD *= L0 /rm;
      BD  = AD /L0;
      
      // Calculate AS and AL
      double logd = log(delta);
      AS  = - (Z[km]/delta + (BD + gamL *psi[km]) *logd) /bS;
      AL  = (psi[km] + delta *(AS + AC + AD *(logd - 1.))) 
	/(1. - lamL *delta *(logd - 1.));
      
      // Perform jump
      for (int k = 0; k < dim; k++)
	{
	  if (_tear == 1 || j > num-1)
	    {
	      // Tearing parity: vacuum always tearing parity
	      if (k == km)
		{
		  Y[k]     += 2. *delta *(psi[k] *lamL *(logd - 1.) + AC + AD *(logd - 1.)) 
		    + 2. *AS *delta;
		  Y[dim+k] += 2. *delta *logd *(psi[k] *gamL             + BD)             
		    + 2. *AS *bS *delta;
		}
	      else
		{
		  Y[k]     += 0.;
		  Y[dim+k] += 0.;
		}
	    }
	  else
	    {
	      // Twisting parity (excluding vacuum)
	      if (k == km)
		{
		  Y[k]     += - 2. *AL + 2. *delta *(AC + AD *(logd - 1.));
		  Y[dim+k] +=            2. *delta *logd *BD;
		}
	      else
		{
		  Y[k]     += - 2. *logd *AL *akh[k];
		  Y[dim+k] += - 2. *logd *AL *bkh[k];
		}
	    }
	}
 
      // Log diagnostic data
      AC = (fabs(AC) < Eta) ? 0. : AC;
      AD = (fabs(AD) < Eta) ? 0. : AD;
      BD = (fabs(BD) < Eta) ? 0. : BD;
      AL = (fabs(AL) < Eta) ? 0. : AL;
      AS = (fabs(AS) < Eta) ? 0. : AS;
      int mcen = (i < dim) ? mpol[i] : Mres[i-dim];
      if (interactive) 
	{
	  if (i == 0)
	    {
	      fprintf (logfile, "nuL = %+10.3e nuS = %+10.3e L1 = %+10.3e P1 = %+10.3e T1 = %+10.3e\n",
		       nuL, nuS, L1, P1, T1);
	      fprintf (logfile, "bS = %+10.3e bL = %+10.3e gamL = %+10.3e lamL = %+10.3e dnuL = %+10.3e dnuS = %10.3e delta = %+10.3e\n",
		       bS, bL, gamL, lamL, dnuL, dnuS, delta);
	    }	
	  fprintf (logfile, "mcen = %+3d AC = %+10.3e AD = %+10.3e BD = %+10.3e AL = %+11.4e AS = %+11.4e\n",
		   mcen, AC, AD, BD, AL, AS);
	  for (int k = 0; k < dim; k++)
	    {
	      int mm = (k < dim) ? mpol[k] : Mres[k-dim];
	      fprintf (logfile, "m = %+3d psi[k] = %+10.3e Z[k] = %+10.3e Pbar[k] = %+10.3e Zbar[k] = %+10.3e akh[k] = %+10.3e bkh[k] = %+10.3e\n",
		       mm, psi[k], Z[k], Pbar[k], Zbar[k], akh[k], bkh[k]);
	    }
	}
      
      // Clean up
      delete[] psi;  delete[] Z;   delete[] akh;  delete[] bkh; 
      delete[] Pbar; delete[] Zbar;
    }

  // ......................
  // Calculate Psi and dPsi
  // ......................
  double Lmm = gsl_matrix_get (Lmat, km, km);
  if (_tear == 1 || j > vac-1)
    {	
      Psi  =   pow(rm, nuL) *sqrt((nuS - nuL) /Lmm) *AL;
      dPsi =   pow(rm, nuS) *sqrt((nuS - nuL) /Lmm) *2. *AS;
    }	
  else
    {
      Psi  = - pow(rm, nuL) *sqrt((nuS - nuL)/Lmm) *AL;
      dPsi = - pow(rm, nuS) *sqrt((nuS - nuL)/Lmm) *2. *AS;
    }

  // ......................................
  // Launch small solution when appropriate
  // ......................................
  if (i >= dim && i-dim == j)
    {
      double Lmm = gsl_matrix_get (Lmat, km, km);
      AS = 1. /pow(rm, nuS) /sqrt((nuS - nuL) /Lmm);
      
      for (int k = 0; k < dim; k++)
	{
	  if (k == km)
	    {
	      Y[k]     = AS     *dnuS;
	      Y[dim+k] = AS *bS *dnuS;
	    }
	  else
	    {
	      Y[k]     = AS *akt[k] *dnuS;
	      Y[dim+k] = AS *bkt[k] *dnuS;
	    }
	}
      dPsi = 1.;
    }

  if (interactive) 
    {
      fprintf (logfile, "Psi = %+10.3e dPsi = %+10.3e\n",
	       Psi, dPsi);
      fclose (logfile);
    }

  // Clean up
  delete[] akt, delete[] bkt; 
}


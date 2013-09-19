// Stage2.cpp

#include "Flux.h"

// ######################################################
// Function to input stage1 data and output Tomuhawc data
// ######################################################
void Flux::Stage2 ()
{
  // ..............
  // Read R, Z grid
  // ..............
  FILE *file = OpenFiler ("Stage1/Points.out");
  if (fscanf (file, "%d %d %d", &II, &JJ, &KK) != 3)
    {
      printf ("Error reading Stage1/Points.out\n");
      exit (1);
    }
  fclose (file);
  printf ("\nConstructing flux coordinate system:\n");

  R = new double [II];  // R array
  Z = new double [JJ];  // Z array

  file = OpenFiler ("Stage1/R.out");
  for (int i = 0; i < II; i++)
    if (fscanf (file, "%lf", &R[i]) != 1)
      {
	printf ("Error reading Stage1/R.out\n");
	exit (1);
      }
  fclose (file);
  file = OpenFiler ("Stage1/Z.out");
  for (int j = 0; j < JJ; j++)
    if (fscanf (file, "%lf", &Z[j]) != 1)
      {
	printf ("Error reading Stage1/Z.out\n");
	exit (1);
      }
  fclose (file);

  // ........
  // Set psib
  // ........
  char *cc = new char[100]; int dum; double sb;
  file = fopen ("profile.in", "r");
  if (fscanf (file, "%s %d",  cc, &dum) != 2)
    {
      printf ("Error reading profile.in (line 1)\n");
      exit (1);
    }
  if (fscanf (file, "%s %d",  cc, &dum) != 2)
    {
      printf ("Error reading profile.in (line 2)\n");
      exit (1);
    }
  if (fscanf (file, "%s %d",  cc, &dum) != 2)
    {
      printf ("Error reading profile.in (line 3)\n");
      exit (1);
    }
  if (fscanf (file, "%s %d",  cc, &dum) != 2)
    {
      printf ("Error reading profile.in (line 4)\n");
      exit (1);
    }
  if (fscanf (file, "%s %d",  cc, &dum) != 2)
    {
      printf ("Error reading profile.in (line 5)\n");
      exit (1);
    }
  if (fscanf (file, "%s %d",  cc, &dum) != 2)
    {
      printf ("Error reading profile.in (line 6)\n");
      exit (1);
    }
  if (fscanf (file, "%s %lf", cc, &sb) != 2)
   {
     printf ("Error reading profile.in (line 7)\n");
      exit (1);
   }
  fclose (file);
  delete[] cc;

  psib = 1. - sb*sb;

  // ........
  // Read Psi
  // ........
  Psi = gsl_matrix_alloc (II, JJ);  // Psi(R, Z)

  double val;
  file = OpenFiler ("Stage1/Psi.out");
  for (int i = 0; i < II; i++)
    for (int j = 0; j < JJ; j++)
      {	
	if (fscanf (file, "%lf", &val) != 1)
	  {
	    printf ("Error reading Stage1/Psi.out\n");
	    exit (1);
	  }
	gsl_matrix_set (Psi, i, j, val);
      }
  fclose (file);
  file = OpenFile ("Stage2/Psi.out");
  for (int i = 0; i < II; i++)
    {
      for (int j = 0; j < JJ; j++)
	{
	  double val = gsl_matrix_get (Psi, i, j);
	  if (val < 0.) 
	    val = 0.;
	  fprintf (file, "%17.10e ", val);
	}
      fprintf (file, "\n");
    }
  fclose (file);
  file = OpenFile ("Stage2/PsiR.out");
  for (int i = 0; i < II; i++)
    {
      for (int j = 0; j < JJ; j++)
	fprintf (file, "%17.10e ", GetPsiR (R[i], Z[j]));
      fprintf (file, "\n");
    }
  fclose (file);
  file = OpenFile ("Stage2/PsiZ.out");
  for (int i = 0; i < II; i++)
    {
      for (int j = 0; j < JJ; j++)
	fprintf (file, "%17.10e ", GetPsiZ (R[i], Z[j]));
      fprintf (file, "\n");
    }
  fclose (file);  

  file = OpenFiler ("Stage1/Axis.out");
  if (fscanf (file, "%lf %lf %lf", &Raxis, &Zaxis, &psic) != 3)
    {
      printf ("Error reading Stage1/Axis.out\n");
      exit (1);
    }	
  fclose (file);

  file = OpenFile ("Stage2/J.out");
  for (int i = 0; i < II; i++)
    {
      for (int j = 0; j < JJ; j++)
	{
	  double val = psic * (GetPsiRR (R[i], Z[j]) + GetPsiZZ (R[i], Z[j])
			       - GetPsiR (R[i], Z[j]) /R[i]) /R[i];
	  if (gsl_matrix_get (Psi, i, j) < 0.5*psib)
	    val = 0.;
	  fprintf (file, "%17.10e ", val);
	}
      fprintf (file, "\n");
    }
  fclose (file);  

  // ..................
  // Find magnetic axis
  // ..................
  double rmin = 1.e6;
  for (int i = 0; i < II; i++)
    if (fabs (R[i] - Raxis) < rmin)
      {
	rmin = fabs (R[i] - Raxis);
	ic   = i;
      }
  rmin = 1.e6;
  for (int j = 0; j < JJ; j++)
    if (fabs (Z[j] - Zaxis) < rmin)
      {
	rmin = fabs (Z[j] - Zaxis);
	jc   = j;
      }

  L = ic+1; // Number of points in Psi(R,0) array

  // .............
  // Read boundary
  // .............
  Rb = new double [KK];   // R (boundary)
  Zb = new double [KK];   // Z (boundary)

  file = OpenFiler ("Stage1/Boundary.out");
  for (int k = 0; k < KK; k++)
    if (fscanf (file, "%lf %lf", &Rb[k], &Zb[k]) != 2)
      {
	printf ("Error reading Stage1/Boundary.out\n");
	exit (1);
      }	
  fclose (file);
  double Rmin = 1.e6, Rmax = -1.e6, Zmin = 1.e6, Zmax = -1.e6;
  for (int k = 0; k < KK; k++)
    {
      if (Rb[k] < Rmin) Rmin = Rb[k];
      if (Rb[k] > Rmax) Rmax = Rb[k];
      if (Zb[k] < Zmin) Zmin = Zb[k];
      if (Zb[k] > Zmax) Zmax = Zb[k];
    }

  // .................
  // Read profile data
  // .................
  double *PSI = new double [II];   // PSI array (EQDSK)
  double *g   = new double [II];   // g(PSI)
  double *p   = new double [II];   // p(PSI)
  double *ggp = new double [II];   // g dg/dPSI
  double *pp  = new double [II];   // dp/dPSI
  double *q   = new double [II];   // q(PSI)

  file = OpenFiler ("Stage1/Profiles.out");
  for (int i = 0; i < II; i++)
    if (fscanf (file, "%lf %lf %lf %lf %lf %lf", &PSI[i], &g[i], &p[i], &ggp[i], &pp[i], &q[i]) != 6)
      {
	printf ("Error reading Stage1/Profiles.out\n");
	exit (1);
      }
  fclose (file);

  // ...................
  // Calculate q profile
  // ...................
 
  // Setup Psi(R,0) profile
  s = new double [L]; // Array of s = sqrt[1 - Psi(R,0)] values

  for (int l = 0; l < L-1; l++)
    s[L-1-l] = sqrt (1. - gsl_matrix_get (Psi, l, jc));
  s[0] = 0.;

  // Setup R(Z=0) profile
  Rs = new double [L]; // Array of R(s) values

  for (int l = 0; l < L-1; l++)
    Rs[L-1-l] = R[l];
  Rs[0] = Raxis;

  file = OpenFile ("Stage2/rs.out");
  for (int l = 0; l  < L; l++)
    fprintf (file, "%17.10e %17.10e\n", Rs[l], s[l]);
  fclose (file);

  // Set up Psi grid
  P   = new double [J];  // Psi array
  RP  = new double [J];  // R(Psi)
  rP  = new double [J];  // r(Psi)
  GP  = new double [J];  // g(Psi)
  QGP = new double [J];  // q(psi)/g(psi) 
  QP  = new double [J];  // q(Psi)
  PP  = new double [J];  // P(Psi)
  GPP = new double [J];  // dg/dPsi
  PPP = new double [J];  // dP/dPsi
  S   = new double [J];  // sqrt(1 - Psi)
  QX  = new double [J];  // q(Psi) from EQDSK

  for (int j = 0; j < J; j++)
    {
      double s = double (j) /double (J);
      P[j] = 1. - s * tanh(s/sc) /tanh(1./sc);
      S[j] = sqrt(1. - P[j]);
    }

  // Calculate g(Psi), P(Psi) profile
  for (int j = 0; j < J; j++)
    {
      double pval = psic * P[j];

      GP [j] = Interpolate (II, PSI, g,   pval, 0);
      PP [j] = Interpolate (II, PSI, p,   pval, 0);
      GPP[j] = Interpolate (II, PSI, ggp, pval, 0) / GP[j];
      PPP[j] = Interpolate (II, PSI, pp,  pval, 0);
      QX [j] = Interpolate (II, PSI, q,   pval, 0);
    }

  // Calculate R(Psi) profile
  for (int j = 0; j < J; j++)
    RP[j] = Interpolate (L, s, Rs, S[j], 0);

  // Calculate q(Psi)/g(Psi) profile
  printf ("Calculating q(Psi)/g(Psi) profile\n");
  CalcQGP ();
  QGP[0] = q[0] /g[0];
  QP [0] = q[0];
  
  // Calculate r(P) profile
  printf ("Calculating r(Psi) profile\n");
  CalcrP ();
  rP[0] = 0.;
  
  file = OpenFile ("Stage2/rpsi1.out");
  for (int j = 0; j  < J; j++)
    fprintf (file, "%17.10e %17.10e\n", RP[j], S[j]);
  fclose (file);

  file = OpenFile ("Stage2/pqg.out");
  for (int j = 0; j < J; j++)
    fprintf (file, "%17.10e %17.10e\n", S[j], QGP[j]);
  fclose (file);

  file = OpenFile ("Stage2/rp.out");
  for (int j = 0; j < J; j++)
    fprintf (file, "%17.10e %17.10e\n", rP[j], S[j]);
  fclose (file);
 
  // Output q(r) profile
  file = OpenFile ("Stage2/qr.out");
  for (int j = 0; j < J; j++)
    fprintf (file, "%17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e\n", 
	     rP[j], QP[j], QGP[j], P[j], GP[j], PP[j], GPP[j], PPP[j], QX[j]);
  fclose (file);

  // Calculate qb
  double rb = epsb/epsa;
  if (fabs(rb-1.) < 1.e-4)
    qb = QA;
  else if (rb > rP[J-1])
    qb = Extrapolate1 (J, rP, QP, rb);
  else
    qb = Interpolate (J, rP, QP, rb, 0);
  file = OpenFile ("Stage2/Edge.out");
  fprintf (file, "%17.10e %17.10e %17.10e\n", psib, rb, qb);
  fclose (file);

  // ........................
  // Calculate straight angle
  // ........................
  printf ("Calculating straight angle\n");
  th  = new double [K];           // Theta array
  Rst = gsl_matrix_alloc (J, K);  // [R(r,theta) - Rc] /r
  Zst = gsl_matrix_alloc (J, K);  // Z(r,theta) /r
  Rt  = gsl_matrix_alloc (J, K);  // (dR/dtheta) /r
  Zt  = gsl_matrix_alloc (J, K);  // (dR/dtheta) /r
  Rr  = gsl_matrix_alloc (J, K);  // (dR/dr) 
  Zr  = gsl_matrix_alloc (J, K);  // (dZ/dr) 

  for (int k = 0; k < K; k++)
    {
      double t = double (k) /double (K-1);
      th[k]    = M_PI * sin (0.5*M_PI * t/tc) /sin (0.5*M_PI * 1./tc);
    }

  Calctst ();
  double zero = 0.;
  for (int k = 0; k < K; k++)
    {
      val = Extrapolate3 (J, rP, Rst, zero, k);
      gsl_matrix_set (Rst, 0, k, val);
      val = Extrapolate3 (J, rP, Zst, zero, k);
      gsl_matrix_set (Zst, 0, k, val);
      val = Extrapolate3 (J, rP, Rt, zero, k);
      gsl_matrix_set (Rt, 0, k, val);
      val = Extrapolate3 (J, rP, Zt, zero, k);
      gsl_matrix_set (Zt, 0, k, val);
    }
  
  // Calculate boundary data
  Boundary ();

  // Output flux coordinates
  file = OpenFile ("Stage2/rx.out");
  for (int j = 0; j < J; j++)
    fprintf (file, "%17.10e ", rP[j]);     
  fprintf (file, "\n");
  fclose (file);
  file = OpenFile ("Stage2/tt.out");
  for (int k = 0; k < K; k++)
    fprintf (file, "%17.10e ", th[k]);  
  fprintf (file, "\n");
  fclose (file);
  file = OpenFile ("Stage2/Rr.out");
  for (int k = 0; k < K; k++)
    {
      for (int j = 0; j < J; j++)
	fprintf (file, "%17.10e ", Raxis + rP[j] * gsl_matrix_get (Rst, j, k));
      fprintf (file, "\n");
    }
  fprintf (file, "\n");
  fclose (file);
  file = OpenFile ("Stage2/Zr.out");
  for (int k = 0; k < K; k++)
    {
      for (int j = 0; j < J; j++)
	fprintf (file, "%17.10e ", rP[j] * gsl_matrix_get (Zst, j, k));
      fprintf (file, "\n");
    }
  fprintf (file, "\n");
  fclose (file);
  file = OpenFile ("Stage2/Rt.out");
  for (int j = 0; j < J; j++)
    {
      for (int k = 0; k < K; k++)
	fprintf (file, "%17.10e ", Raxis + rP[j] * gsl_matrix_get (Rst, j, k));
    }
  fprintf (file, "\n");
  fclose (file);
  file = OpenFile ("Stage2/Rtasy.out");
  for (int j = 0; j < J; j++)
    {
      for (int k = 0; k < K; k++)
	fprintf (file, "%17.10e ", Raxis + rP[j] * gsl_matrix_get (Rst, j, k));
      fprintf (file, "\n");
    }
  fprintf (file, "\n");
  fclose (file);
  file = OpenFile ("Stage2/Zt.out");
  for (int j = 0; j < J; j++)
    {
      for (int k = 0; k < K; k++)
	fprintf (file, "%17.10e ", rP[j] * gsl_matrix_get (Zst, j, k));
    }
  fprintf (file, "\n");
  fclose (file);
  file = OpenFile ("Stage2/Ztasy.out");
  for (int j = 0; j < J; j++)
    {
      for (int k = 0; k < K; k++)
	fprintf (file, "%17.10e ", rP[j] * gsl_matrix_get (Zst, j, k));
      fprintf (file, "\n");
    }
  fprintf (file, "\n");
  fclose (file);

  // ...................
  // Calculate Rr and Zr
  // ...................
  double *Rval = new double [J];
  double *Zval = new double [J];

  for (int k = 0; k < K; k++)
    {
      for (int j = 0; j < J; j++)
	{
	  Rval[j] = Raxis + rP[j] * gsl_matrix_get (Rst, j, k);
	  Zval[j] =         rP[j] * gsl_matrix_get (Zst, j, k);
	}

      for (int j = 0; j < J; j++)
      	{
	  val = Interpolate (J, rP, Rval, rP[j], 1);
      	  gsl_matrix_set (Rr, j, k, val);   
	  val = Interpolate (J, rP, Zval, rP[j], 1);
      	  gsl_matrix_set (Zr, j, k, val); 
      	}
    }
  delete[] Rval; delete[] Zval;

  // ......................
  // Output metric elements
  // ......................
  igradr2 = gsl_matrix_alloc (J, K); // |grad r|^{-2}
  gradrt  = gsl_matrix_alloc (J, K); // r grad r . grad th / |grad r|^{-2}
  R2      = gsl_matrix_alloc (J, K); // (R^2 - R_c^2) /r

  file = OpenFile ("Stage2/r.out");
  for (int k = 0; k < K; k++)
    {
      for (int j = 0; j < J; j++)
	{
	  double r = gsl_matrix_get (Rst, j, k);

	  fprintf (file, "%17.10e ", r);
	}
      fprintf (file, "\n");
    }
  fclose (file);
  file = OpenFile ("Stage2/z.out");
  for (int k = 0; k < K; k++)
    {
      for (int j = 0; j < J; j++)
	{
	  double z = gsl_matrix_get (Zst, j, k);

	  fprintf (file, "%17.10e ", z);
	}
      fprintf (file, "\n");
    }
  fclose (file);
  file = OpenFile ("Stage2/rrt.out");
  for (int k = 0; k < K; k++)
    {
      for (int j = 0; j < J; j++)
	{
	  double rt = gsl_matrix_get (Rt, j, k);

	  fprintf (file, "%17.10e ", rt);
	}
      fprintf (file, "\n");
    }
  fclose (file);
  file = OpenFile ("Stage2/zzt.out");
  for (int k = 0; k < K; k++)
    {
      for (int j = 0; j < J; j++)
	{
	  double zt = gsl_matrix_get (Zt, j, k);

	  fprintf (file, "%17.10e ", zt);
	}
      fprintf (file, "\n");
    }
  fclose (file);
  file = OpenFile ("Stage2/rrr.out");
  for (int k = 0; k < K; k++)
    {
      for (int j = 0; j < J; j++)
	{
	  double rr = gsl_matrix_get (Rr, j, k);

	  fprintf (file, "%17.10e ", rr);
	}
      fprintf (file, "\n");
    }
  fclose (file);
  file = OpenFile ("Stage2/zzr.out");
  for (int k = 0; k < K; k++)
    {
      for (int j = 0; j < J; j++)
	{
	  double zr = gsl_matrix_get (Zr, j, k);

	  fprintf (file, "%17.10e ", zr);
	}
      fprintf (file, "\n");
    }
  fclose (file);
  file = OpenFile ("Stage2/igradr2.out");
  for (int k = 0; k < K; k++)
    {
      for (int j = 0; j < J; j++)
	{
	  double Rv = Raxis + rP[j] * gsl_matrix_get (Rst, j, k);
	  double rt = gsl_matrix_get (Rt, j, k);
	  double zt = gsl_matrix_get (Zt, j, k);
	  double rr = gsl_matrix_get (Rr, j, k);
	  double zr = gsl_matrix_get (Zr, j, k);

	  val = Rv*Rv *epsa*epsa*epsa*epsa /(rt*rt + zt*zt);
	  gsl_matrix_set (igradr2, j, k, val);

	  fprintf (file, "%17.10e ", val);
	}
      fprintf (file, "\n");
    }
  fclose (file);
  file = fopen ("Stage2/gradrt.out", "w");
  for (int k = 0; k < K; k++)
    {
      for (int j = 0; j < J; j++)
	{
	  double rt = gsl_matrix_get (Rt, j, k);
	  double zt = gsl_matrix_get (Zt, j, k);
	  double rr = gsl_matrix_get (Rr, j, k);
	  double zr = gsl_matrix_get (Zr, j, k);

	  val = - (zr*zt + rr*rt) /(zt*zt + rt*rt);
	  gsl_matrix_set (gradrt, j, k, val);

	  fprintf (file, "%17.10e ", val);
	}
      fprintf (file, "\n");
    }
  fclose (file);
  file = fopen ("Stage2/R2.out", "w");
  for (int k = 0; k < K; k++)
    {
      for (int j = 0; j < J; j++)
	{
	  double R1 = Raxis + rP[j] * gsl_matrix_get (Rst, j, k);

	  if (j == 0)
	    val = 2. * Raxis * gsl_matrix_get (Rst, j, k);
	  else
	    val = (R1*R1 - Raxis*Raxis) /rP[j];
	  gsl_matrix_set (R2, j, k, val);

	  fprintf (file, "%17.10e ", val);
	}
      fprintf (file, "\n");
    }
  fclose (file);
  double Jac = 0.;
  file = fopen ("Stage2/Jac.out", "w");
  for (int k = 0; k < K; k++)
    {
      for (int j = 0; j < J; j++)
	{
	  double R1 = Raxis + rP[j] * gsl_matrix_get (Rst, j, k);

	  double rt = gsl_matrix_get (Rt, j, k);
	  double zt = gsl_matrix_get (Zt, j, k);
	  double rr = gsl_matrix_get (Rr, j, k);
	  double zr = gsl_matrix_get (Zr, j, k);

	  val = (rt*zr - rr*zt) /R1 /epsa/epsa;
	  Jac += val;

	  fprintf (file, "%17.10e ", val);
	}
      fprintf (file, "\n");
    }
  fclose (file);
  Jac /= double (K*J);

  // ...........................
  // Calculate profile functions
  // ...........................
  printf ("Calculating profile functions\n");
  Profile ();

  // ...........................
  // Calculate metric functions
  // ...........................
  printf ("Calculating metric functions\n");
  Metric ();

  // ...........
  // Output data
  // ...........
  printf ("II = %3d  JJ = %3d  KK = %3d  ic = %3d  jc = %3d\n", 
	  II, JJ, KK, ic, jc);
  printf ("psib = %11.4e  a = %11.4e  b = %11.4e  b/a = %11.4e  qb = %11.4e\n", 
	  psib, epsa, epsb, rb, qb);
  printf ("Jacobian residual = %11.4e\n\n", 
	  fabs(Jac-1.));

  file = OpenFile ("Stage2/info.out");
  fprintf (file, "%17.10e %17.10e %17.10e %17.10e\n", Gamma, 1./IASPCT, epsa, epsb/epsa);
  fclose (file);

  // ........
  // Clean up
  // ........
  delete[] Rb;  delete[] Zb;  delete[] R;  delete[] Z;
  delete[] PSI; delete[] g;   delete[] p;  delete[] ggp; 
  delete[] s;   delete[] Rs;  delete[] pp; delete[] q;
  delete[] P;   delete[] RP;  delete[] rP; delete[] th;  
  delete[] S;   delete[] GP;  delete[] PP; delete[] GPP;
  delete[] PPP; delete[] QGP; delete[] QP; delete[] QX; 
  gsl_matrix_free (Rt);      gsl_matrix_free (Zt);     gsl_matrix_free (Rst);
  gsl_matrix_free (Zst);     gsl_matrix_free (Rr);     gsl_matrix_free (Zr);
  gsl_matrix_free (igradr2); gsl_matrix_free (gradrt); gsl_matrix_free (Psi);     
  gsl_matrix_free (R2);
}

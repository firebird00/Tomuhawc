import graph;
import palette;
import contour;

file    in4 = input("../Stage1/Boundary.out").line();
real[][] A4 = in4.dimension (0,0);
A4          = transpose(A4);
real[] R    = A4[0];
real[] Z    = A4[1];
     
file    in2 = input("../Stage1/Box.out").line();
real[][] A2 = in2.dimension (0,0);
A2          = transpose(A2);
real[] aa   = A2[0];
real rmin   = (real) aa[0];
real[] bb   = A2[1];
real zmin   = (real) bb[0]; 
real[] cc   = A2[2];
real rmax   = (real) cc[0]; 
real[] dd   = A2[3];
real zmax   = (real) dd[0]; 
real zsize  = 400;
real rsize  = 800*(rmax-rmin)/(zmax-zmin);

size(rsize,zsize,(rmin,0.),(rmax,zmax));

file    in3 = input("../Stage1/Axis.out").line();
real[][] A3 = in3.dimension (0,0);
A3          = transpose(A3);
real[] xx   = A3[0];
real rw     = (real) xx[0];
real[] yy   = A3[1];
real zw     = (real) yy[0]; 
     
file    inr   = input("Rr.out").line();
real[][]  Ar  = inr.dimension (0,0);
int K         = (Ar[0]).length;
          Ar  = transpose(Ar); 

file    inz   = input("Zr.out").line();
real[][] Az   = inz.dimension (0,0);
         Az   = transpose(Az); 

file     inrt = input("Rtasy.out").line();
real[][]  Art = inrt.dimension (0,0);
int Kt        = (Art[0]).length;
          Art = transpose(Art); 

file    inzt  = input("Ztasy.out").line();
real[][] Azt  = inzt.dimension (0,0);
         Azt  = transpose(Azt); 

pen s;

s = solid+red+0.25;
for (int k = 0; k < Kt; ++k)
{
draw(graph(Art[k],Azt[k]),s);		
}
s = solid+blue+0.25;
for (int k = 0; k < K; ++k)
{
draw(graph(Ar[k],Az[k]),s);		
}

dot ((rw, zw));

s = black + 2;
draw ((rmin,0.)--(rmin,zmax)--(rmax,zmax)--(rmax,0.),s);

s = black + 0.5;
draw(graph(R,Z),s);

limits ((0.99*rmin,-0.01*zmax),(1.01*rmax,1.01*zmax),Crop);

pen q = fontsize(20.);
defaultpen (q);
xaxis("$R$",  Bottom, RightTicks, above=true);
xaxis(Top, NoTicks, above=true);
yaxis("$Z$", Right, RightTicks, above=true);
yaxis(Left, NoTicks, above=true);

import graph;
import palette;
import contour;
   
file    in = input("J.out").line();
real[][] a = in.dimension (0,0);
a=-a;
  
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
real zsize  = 500;
real rsize  = 500*(rmax-rmin)/(zmax-zmin);

file    in4 = input("../Stage1/Boundary.out").line();
real[][] A4 = in4.dimension (0,0);
A4          = transpose(A4);
real[] rb   = A4[0];
real[] zb   = A4[1];

file    in3 = input("../Stage1/Axis.out").line();
real[][] A3 = in3.dimension (0,0);
A3          = transpose(A3);
real[] xx   = A3[0];
real rw     = (real) xx[0];
real[] yy   = A3[1];
real zw     = (real) yy[0]; 

size(rsize,zsize,(rmin,zmin),(rmax,zmax));

pen[] Palette = BWRainbow();
bounds range  = image (a, Full, (rmin,zmin), (rmax,zmax), Palette);

pen s = black+3.;
draw ((rmin,zmin)--(rmax,zmin)--(rmax,zmax)--(rmin,zmax)--(rmin,zmin),s);

s = red+1;
draw(graph(rb,zb),s);

real[] cvals = uniform(range.min,range.max,10);
draw(contour(a, (rmin,zmin), (rmax,zmax), cvals));

file    in4 = input("../Stage2/Edge.out").line();
real[][] A4 = in4.dimension (0,0);
A4          = transpose(A4);
real[] xxx  = A4[0];
real psib   = (real) xxx[0];

in = input("J.out").line();
a = in.dimension (0,0);

file    inx = input("../Stage1/Psi.out").line();
real[][] ax = inx.dimension (0,0);

s = blue + 1;
real[] cvals = {psib};
draw(contour(ax, (rmin,zmin), (rmax,zmax), cvals),s);

s=white;
dot ((rw, zw),s);

pen q = fontsize(20.);
defaultpen (q);
xaxis("$R$",  Bottom, RightTicks, above=true);
xaxis(Top, NoTicks, above=true);
yaxis("$Z$", Right, RightTicks, above=true);
yaxis(Left, NoTicks, above=true);

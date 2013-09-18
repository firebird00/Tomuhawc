import graph;
import palette;
import contour;
     
file    in2 = input("Box.out").line();
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

size(rsize,zsize,(rmin,zmin),(rmax,zmax));

file    in1 = input("Boundary.out").line();
real[][] A1 = in1.dimension (0,0);
A1          = transpose(A1);
     
real[] r  = A1[0];
real[] z  = A1[1];

file    in3 = input("Axis.out").line();
real[][] A3 = in3.dimension (0,0);
A3          = transpose(A3);
     
real[] Raxis = A3[0];
real[] Zaxis = A3[1];
real raxis   = (real) Raxis[0];
real zaxis   = (real) Zaxis[0];

file    in = input("psi.out").line();
real[][] a = in.dimension (0,0);
     
pen[] Palette = Rainbow();
bounds range  = image (a, Full, (rmin,zmin), (rmax,zmax), Palette);

real[] cvals = uniform(0.,range.max,10);
draw(contour(a, (rmin,zmin), (rmax,zmax), cvals));

dot ((raxis,zaxis));

pen s = solid + red + 1.;
draw(graph(r,z),s);

pen q = fontsize(20.);
defaultpen (q);
xaxis("$R$", Bottom, RightTicks, above=true);
xaxis(Top,  NoTicks, above=true);
yaxis("$Z$", Right, RightTicks, above=true);
yaxis(Left, NoTicks, above=true);

s = black+3.;
draw ((rmin,zmin)--(rmax,zmin)--(rmax,zmax)--(rmin,zmax)--(rmin,zmin),s);
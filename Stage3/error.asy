import graph;
     
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

size(rsize,zsize,(rmin,zmin),(rmax,zmax));

file    in1 = input("error.out").line();
real[][] A1 = in1.dimension (0,0);
A1          = transpose(A1);
     
real[] rb  = A1[0];
real[] zb  = A1[1];
real[] rx  = A1[2];
real[] zx  = A1[3];
real[] rp  = A1[4];
real[] zp  = A1[5];
real[] rm  = A1[6];
real[] zm  = A1[7];

pen s = dashed + black + 1.;
draw(graph(rb,zb),s);
s=dotted+black+2;
draw(graph(rx,zx),s);
s = solid + black + 1.;
draw(graph(rp,zp),s);
s = solid + black + 1.;
draw(graph(rm,zm),s);

dot((1.,0.));

limits ((0.95*rmin,1.1*zmin), (1.05*rmax,1.1*zmax), Crop);

pen q = fontsize(20.);
defaultpen (q);
xaxis("$R/R_0$", Bottom, RightTicks, above=true);
xaxis(Top,  NoTicks, above=true);
yaxis("$Z/R_0$", Right, RightTicks, above=true);
yaxis(Left, NoTicks, above=true);


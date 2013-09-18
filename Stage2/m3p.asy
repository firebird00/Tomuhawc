import graph;
     
size(500,500,IgnoreAspect);

file    in1 = input("profile.out").line();
real[][] A1 = in1.dimension (0,0);
A1          = transpose(A1);
     
real[] r = A1[0];

file    in = input("M3p.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] m0 = A[0];

pen s = solid + 1.5;
draw(graph(r,m0),s,marker(scale(1.5mm)*polygon(4)));

s = dotted + 1.;
xequals(1.,s);

s = solid + 1.;
file    in4 = input("../Stage2/Edge.out").line();
real[][] A4 = in4.dimension (0,0);
A4          = transpose(A4);
real[] xxx  = A4[1];
real b   = (real) xxx[0];
xequals (b, s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r$",BottomTop,LeftTicks);
yaxis("$d(M_3)/dr$",LeftRight,RightTicks);

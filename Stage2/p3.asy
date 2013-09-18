import graph;
     
size(500,500,IgnoreAspect);

file    in = input("profile.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] r = A[0];
real[] q = A[5];

pen s = solid + 1.;
draw(graph(r,q),s,marker(scale(1.5mm)*polygon(4)));

file    in4 = input("../Stage2/Edge.out").line();
real[][] A4 = in4.dimension (0,0);
A4          = transpose(A4);
real[] xxx  = A4[1];
real b   = (real) xxx[0];
xequals (b, s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r$",BottomTop,LeftTicks);
yaxis("$\alpha_p$",LeftRight,RightTicks);

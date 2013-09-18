import graph;
     
size(500,500,IgnoreAspect);

file    in = input("Profiles.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] p  = A[0];
real[] g  = A[3];

pen s = solid + 1.5;
draw(graph(p,g),s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$\psi$",BottomTop,LeftTicks);
yaxis("$g\,dg/d\psi$",LeftRight,RightTicks);

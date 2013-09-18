import graph;
     
size(500,500,IgnoreAspect);

file    in = input("Profiles.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] p  = A[0];
real[] g  = A[5];

pen s = solid + 1.5;
draw(graph(p,g),s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$\Psi$",BottomTop,LeftTicks);
yaxis("$q$",LeftRight,RightTicks);

import graph;
     
size(500,500,IgnoreAspect);

file    in = input("pqg.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);

     
real[] r  = A[0];
real[] p  = A[1];

pen s = solid + 1.5;
draw(graph(r,p),s,marker(scale(1.5mm)*polygon(4)));

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$(1-\Psi)^{1/2}$",BottomTop,LeftTicks);
yaxis("$q/g$",LeftRight,RightTicks);

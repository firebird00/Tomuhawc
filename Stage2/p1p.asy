import graph;
     
size(500,500,IgnoreAspect);

file    in = input("profile.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] r = A[0];
real[] q = A[9];

pen s = solid + 1.;
draw(graph(r,q),s,marker(scale(1.5mm)*polygon(4)));

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r$",BottomTop,LeftTicks);
yaxis("$d\alpha_\epsilon/dr$",LeftRight,RightTicks);

import graph;
     
size(500,500,IgnoreAspect);

file    in = input("ldete.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] r = A[0];
real[] h = A[1];

pen s = solid + 1.;
draw(graph(r,h),s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$\log(\Delta)$",LeftRight,RightTicks);

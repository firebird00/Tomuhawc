import graph;
     
size(500,500,IgnoreAspect);

file    in = input("rs.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);

file    in1 = input("rpsi1.out").line();
real[][] A1 = in1.dimension (0,0);
A1          = transpose(A1);
     
real[] r  = A[0];
real[] p  = A[1];

real[] r1  = A1[0];
real[] p1  = A1[1];

pen s = solid + 1.5;
draw(graph(p,r),s);

s = solid + 0 + red;
draw(graph(p1,r1),s,marker(scale(1.5mm)*polygon(4)));	

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$(1-\Psi)^{1/2}$",BottomTop,LeftTicks);
yaxis("$R$",LeftRight,RightTicks);

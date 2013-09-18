import graph;
     
size(500,500,IgnoreAspect);

file    in = input("GGJ.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] r = A[0];
real[] d = A[4];

pen s;

file    in0  = input("info.out").line();
real[][] A0  = in0.dimension (0,0);
A0           = transpose(A0);

real[] x0 = A0[0];
real[] x1 = A0[1];
real[] x2 = A0[2];
real[] x3 = A0[3];
real[] x4 = A0[4];
int dim = (int)  x0[0];
int num = (int)  x1[0];
int vac = (int)  x2[0];
int xff = (int)  x3[0];
real b  = (real) x4[0];

s = solid + 1.;
draw(graph(r,d),s);

xlimits(0.,1.,Crop);

s = dotted + 0.5;
yequals(0.,s);

// Read in rational surface data
file    in1 = input("ratsur.out").line();
real[][] A1 = in1.dimension (0,0);
A1          = transpose(A1);

s = black + dashed + 1.;
for (int i = 0; i < vac; ++i)
xequals(A1[2*i][0],s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r$",BottomTop,LeftTicks);
yaxis("$K\,F+G$",LeftRight,RightTicks);

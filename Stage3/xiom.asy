import graph;
     
size(500,500,IgnoreAspect);

// Read in dim and num
file    in  = input("info.out").line();
real[][] A0 = in.dimension (0,0);
A0          = transpose(A0);

real[] x0 = A0[0];
real[] x1 = A0[1];
real[] x2 = A0[2];
int dim = (int) x0[0];
int num = (int) x1[0];
int vac = (int) x2[0];

// Read in data
file    in = input("xio.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);

file    in1 = input("modes.out").line();
real[][] A1 = in1.dimension (0,0);
A1          = transpose(A1);

// Determine rational surface
int res;
if (num == 1)
{
res = 1;
}
else
{
res = getint(1,"rational surface number (1,2,3, etc) ? ");
}
int off = (res-1)*dim;

real[] x = new real[dim];
real[] y = new real[dim];

real sum1 = 0.;
for (int i = 0; i < dim; ++i)
{
real fac;
sum1 += A[i+off][0]* A[i+off][0];
}

for (int i = 0; i < dim; ++i)
{
x[i] = A1[i+off][0];
y[i] = A[i+off][0] /sqrt(sum1);
}

pen s = solid + 0.5;
draw (graph(x,y),s,marker(scale(1.5mm)*polygon(4)));

s = black + dotted + 1;
yequals(0.,s);
xequals(0.,s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$m$",BottomTop,LeftTicks);
yaxis("$\xi_m^o$",LeftRight,RightTicks);

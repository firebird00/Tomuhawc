import graph;
     
size(500,500,IgnoreAspect);

// Read in dim and num
file    in  = input("info.out").line();
real[][] A0 = in.dimension (0,0);
A0          = transpose(A0);

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

// Read in data
file    in = input("psiFe.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
    
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

int skip = (int)((dim - vac)/2) - xff;

pen s = solid + black + 1.;
for (int i = 0; i < num; ++i)
{
draw(graph(A[0],A[1+i+skip+off]),s);
}

s = dashed + green + 1.;
draw(graph(A[0],A[1+res-1+skip+off]),s);

s = solid + red + 0.5;
for (int i = 0; i < skip; ++i)
{
draw(graph(A[0],A[1+i+off]),s);
}

s = solid + blue + 0.5;
for (int i = 0; i < skip; ++i)
{
draw(graph(A[0],A[1+i+num+skip+off]),s);
}

// Read in rational surface data
file    in1 = input("ratsur.out").line();
real[][] A1 = in1.dimension (0,0);
A1          = transpose(A1);

s = black + dashed + 0.5;
for (int i = 0; i < vac; ++i)
xequals(A1[2*i][0],s);
s = solid + 1;
xequals (b, s);

s = dotted + 1;
yequals(0.,s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r/a$",BottomTop,LeftTicks);
yaxis("$\psi$",LeftRight,RightTicks);

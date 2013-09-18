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

write("dim = ", dim); write("num = ", num);

// Read in data
file    in = input("solne.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);

int res = getint(1,"solution no. ? ");

pen s = solid + black + 1.;
for (int i = 0; i < dim; ++i)
{
draw(graph(A[0],A[1+2*dim*res+dim+i]),s);
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

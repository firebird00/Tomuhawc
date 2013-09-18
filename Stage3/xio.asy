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

real xicos(real t)
{
real tt = t*pi/180.;

real sum1 = 0.;
for (int i = 0; i < dim; ++i)
{
sum1 += A[i+off][0]* A[i+off][0];
}

real sum = 0.;
for (int i = 0; i < dim; ++i)
{
sum += A[i+off][0]*cos(A1[i+off][0]*tt);
}
return sum/sum1;
}

real xisin(real t)
{
real tt = t*pi/180.;

real sum1 = 0.;
for (int i = 0; i < dim; ++i)
{
sum1 += A[i+off][0]* A[i+off][0];
}

real sum = 0.;
for (int i = 0; i < dim; ++i)
{
sum += A[i+off][0]*sin(A1[i+off][0]*tt);
}
return sum/sum1;
}

real xi(real t)
{
real xc = xicos(t);
real xs = xisin(t);
return sqrt(xc*xc+xs*xs);
}

real xi1(real t)
{
return -xi(t);
}

pen s = solid + dashed + 1.;
draw (graph(xi,0.,360.,1000),s);
draw (graph(xi1,0.,360.,1000),s);
pen s = solid + red + 1.;
draw (graph(xicos,0.,360.,1000),s);
s = solid + blue + 1.;
draw (graph(xisin,0.,360.,1000),s);

xlimits(0.,360.,Crop);
s = black + dotted + 1;
yequals(0.,s);
xequals(180.,s);	
xequals(90.,s);	
xequals(270.,s);	

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$\theta(^\circ)$",BottomTop,LeftTicks);
yaxis("$\xi^o$",LeftRight,RightTicks);

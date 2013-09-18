import graph;
     
size(500,500,IgnoreAspect);

file    in1 = input("profile.out").line();
real[][] A1 = in1.dimension (0,0);
A1          = transpose(A1);
     
real[] r = A1[0];

file    in = input("M5.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);

int k = getint (1, "k ? ");
write(0,1*k,2*k,3*k,4*k,5*k,6*k,7*k,8*k,9*k);
     
real[] m0 = A[0];
real[] m1 = A[1*k];
real[] m2 = A[2*k];
real[] m3 = A[3*k];
real[] m4 = A[4*k];
real[] m5 = A[5*k];
real[] m6 = A[6*k];
real[] m7 = A[7*k];
real[] m8 = A[8*k];
real[] m9 = A[9*k];
          
pen s = solid + 1.5;
draw(graph(r,m0),s);
s = red + 0.5;
draw(graph(r,m1),s);
s = green + 0.5;
draw(graph(r,m2),s);
s = blue+ 0.5;
draw(graph(r,m3),s);
s = dashed+red + 0.5;
draw(graph(r,m4),s);
s = dashed+green + 0.5;
draw(graph(r,m5),s);
s = dashed+blue+ 0.5;
draw(graph(r,m6),s);
s = dotted+red + 0.5;
draw(graph(r,m7),s);
s = dotted+green + 0.5;
draw(graph(r,m8),s);
s = dotted+blue+ 0.5;
draw(graph(r,m9),s);

s = dotted + 1.;
yequals(0.,s);
xequals(1.,s);

s = solid + 1.;
file    in4 = input("../Stage2/Edge.out").line();
real[][] A4 = in4.dimension (0,0);
A4          = transpose(A4);
real[] xxx  = A4[1];
real b   = (real) xxx[0];
xequals (b, s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$r$",BottomTop,LeftTicks);
yaxis("$M_5$",LeftRight,RightTicks);

import graph;
     
size(500,500,IgnoreAspect);

file    in = input("scanqi.out").line();
real[][] A = in.dimension (0,0);
A          = transpose(A);
     
real[] qr = A[0];
real[] qi = -A[1];
real[] dr = A[2];
real[] di = A[3];
real[] dx = sqrt(dr*dr+di*di);

pen s = solid + 1.;
draw(graph(qi,dr, qi>0),s);
s = dashed + 1.;
draw(graph(qi,di, qi>0),s);
s = dotted + 2.;
draw(graph(qi,dx, qi<0),s);
s = dotted + 1;
xequals (0.,s);
yequals (0.,s);

pen qq = fontsize(25.);
defaultpen (qq);
xaxis("$\tilde{\mit\Omega}_1$",BottomTop,LeftTicks);
yaxis("$\tilde{\Delta}_1^e$",LeftRight,RightTicks);

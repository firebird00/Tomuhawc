import graph;
import palette;
import contour;

size(500,500,IgnoreAspect);
          
file    in1 = input("qr.out").line();
real[][] A1 = in1.dimension (0,0);
A1          = transpose(A1);
     
real[] r  = A1[0];

file    in = input("zr.out").line();
real[][] A = in.dimension (0,0);

int    j = getint(0, "j ? ");
real[] q = A[j];
     
pen s = solid + 1. + green;
draw(graph(r,q),s,MarkFill[0]);

pen q = fontsize(20.);
defaultpen (q);
xaxis("$r$",BottomTop,LeftTicks);
yaxis("$Z/r$",LeftRight,RightTicks);

import graph;
import palette;
import contour;

size(500,500,IgnoreAspect);
          
file    in1 = input("../Stage1/R.out").line();
real[][] A1 = in1.dimension (0,0);
A1          = transpose(A1);
     
real[] r  = A1[0];

file    in = input("J.out").line();
real[][] A = in.dimension (0,0);

int    j = getint(0, "j ? ");
real[] q = A[j];
    
pen s = solid + 1.;
draw(graph(r,q),s,marker(scale(1.5mm)*polygon(4)));

s = dotted +1.;
yequals (0.);

pen q = fontsize(20.);
defaultpen (q);
xaxis("$R$",BottomTop,LeftTicks);
yaxis("$J$",LeftRight,RightTicks);

import graph;
import palette;
     
size(500,500,Aspect);
     
file    in = input("Vmatres.out").line();
real[][] a = in.dimension (0,0);
     
int m = (int)(a[0].length);

pen[] Palette = BWRainbow();

image (a, Automatic, (0,0), (m,m), Palette);

limits ((0,0), (m,m), Crop);

pen q = fontsize(18.);
defaultpen (q);
xaxis("$m$",  BottomTop, LeftTicks, above=true);
yaxis("$m'$", LeftRight, RightTicks, above=true);
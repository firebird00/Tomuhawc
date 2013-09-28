import graph;
import palette;
     
size(500,500,Aspect);
     
file    in = input("Ps.out").line();
real[][] a = in.dimension (0,0);
     
int m = (int)(a[0].length);

real[][] b = new real[m][m];
for (int i = 0; i < m; ++i)
for (int j = 0; j < m; ++j)
b[i][j] = -fabs(a[i][j]);

pen[] Palette = BWRainbow();

image (b, Automatic, (0,0), (m,m), Palette);

limits ((0,0), (m,m), Crop);

pen q = fontsize(18.);
defaultpen (q);
xaxis("$m$",  BottomTop, LeftTicks, above=true);
yaxis("$m'$", LeftRight, RightTicks, above=true);
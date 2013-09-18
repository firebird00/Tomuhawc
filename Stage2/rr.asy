import graph;
import palette;
import contour;

size(500,500,IgnoreAspect);

file    in = input("rrr.out").line();
real[][] a = in.dimension (0,0);
     
pen[] Palette = Rainbow();
bounds range  = image (a, Full, (0,0), (1,1), Palette);

real[] cvals = uniform(range.min,range.max,20);
draw(contour(a, (0,0), (1,1), cvals));

limits ((0.,0.), (1.,1.), Crop);

pen q = fontsize(20.);
defaultpen (q);
xaxis("$\theta/\pi$",  Bottom, RightTicks, above=true);
xaxis(Top, NoTicks, above=true);
yaxis("$r$", Right, RightTicks, above=true);
yaxis(Left, NoTicks, above=true);


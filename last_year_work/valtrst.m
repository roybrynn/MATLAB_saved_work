h = 0.01;
[Xvals,Yvals] = meshgrid(0:h:1, 0:h:1);
g = meshgrid(0:0.01:1);

fun = @(x, y) x + y;

hold = (Xvals.*Xvals) + (Yvals.*Yvals);

int = (trapz(trapz(hold)))*(h*h)


function [gb] = Int2DPt2(xmax, xmin, ymax, ymin)
    [Xvals,Yvals] = meshgrid(0:0.01:1);

    fun = @(x, y) norm((Xvals-Yvals));
    all_values = integral2(fun, xmin, xmax, ymin, ymax);
    disp(all_values)
    gb = all_values(0:100);
end
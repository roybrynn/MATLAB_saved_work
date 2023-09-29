function [g] = NumericalTest(xmin,xmax,ymin,ymax)
%the step size
h = 0.01;

%x and y lists for calculating the integral
xmm = xmin:(xmax-xmin)/100:xmax;
ymm = ymin:(ymax-ymin)/100:ymax;

%the 2D matrices for the corresponding x and y values for the boundary
[Xvals,Yvals] = meshgrid(xmm, ymm);

%the plain grid for vals that will be altered
vals = meshgrid(xmm, ymm);

%for loop to that goes over the 2D matrix
for k = 1:size(Xvals,1)
    for l = 1:size(Yvals,1)
        %stores the function value in the corresponding spot
        vals(k,l) = ((Xvals(k,l))*(Yvals(k,l))); %function values
    end
end
%returns the calculated integral
g = trapz(ymm, trapz(xmm, vals, 2));

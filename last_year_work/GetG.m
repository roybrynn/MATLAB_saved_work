function [g] = GetG(xmin,xmax,ymin,ymax)
%num steps
h = 9;

%x and y lists for calculating the integral
xmm = xmin:(xmax-xmin)/h:xmax;
ymm = ymin:(ymax-ymin)/h:ymax;

%the 2D matrices for the corresponding x and y values for the boundary
[Xvals,Yvals] = meshgrid(xmm, ymm);

%the plain grids for g & vals that will be altered
g = meshgrid(xmm, ymm);
vals = meshgrid(xmm, ymm);

%for loop to that goes over the 2D matrix
for i = 1:size(Xvals,1)
    for j = 1:size(Yvals,1)
        for k = 1:size(Xvals,1)
            for l = 1:size(Yvals,1)
                %stores the function value in the corresponding spot
                %(((Xvals(k,l) - Yvals(k,l))^2) + ((Xvals(i,j) - Yvals(i,j))^2))
                vals(k,l) = (((Xvals(k,l) - Yvals(k,l))^2) + ((Xvals(i,j) - Yvals(i,j))^2));
            end
        end
        %saving the integral on the boundary to the corresponding spot in g
        %vals = vals.*vals;

        g(i, j) = trapz(ymm, trapz(xmm, vals, 2));
        %reseting vals to save the next values
        vals = meshgrid(xmm, ymm);
    end
end
h = 0.01;

[Xvals,Yvals] = meshgrid(0:h:1, 0:h:1);
g = meshgrid(0:0.01:1);

xmin = 0.25;
xmax = 0.75;
ymin  = 0.25;
ymax = 0.75;
% the original file

vals =meshgrid(0:0.01:1);

%calculate the integral for each matrix of one of the fixed IJs
%store that value in G

for i = 1:size(Xvals,1)
    for j = 1:size(Yvals,1)
        for k = 1:size(Xvals,1)
            for l = 1:size(Yvals,1)
                %store in matrix, not vals
                vals(k,l) = ((Xvals(i,j) - Xvals(k,l))^2) + ((Yvals(i,j) - Yvals(k,l))^2); %function value
                %hold = hold + ((val, xmin, xmax, ymin, ymax)); %integral value
                %integral is over the boundary? maybe not the xmin/max
                %am I calculating this correctly?
            end
        end
        g(i, j) = trapz(trapz(vals))*h*h;
        vals = meshgrid(0:0.01:1);
    end
end
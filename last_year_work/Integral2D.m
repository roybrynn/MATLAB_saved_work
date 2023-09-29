function [g] = Integral2D(bxmin, bxmax, bymin, bymax)
    x = 0:0.01:1;
    y = 0:0.01:1;
    g = 0:0.01:4;
    count = 1;
    vals = 1:1:4;
    vals(1:4) = 0;
    %change based off of x
    %mesh, how is that used for integration?
    while (count < 401)
       fun = @(x,y) abs(x-y);
       g(count) = integral2(fun, bxmin, bxmax, bymin, bymax);
       count = count + 1;
    end
end

%My problem is that I do not know what the function is supposed to do, I
%need clarity on what input should be used and how the integral should be
%calculated. When I have calculated the integral none of the values change.
%Am I just supposed to be calculating g for different values of B and
%storing them in the vector from 1-400?

%I was able to loop through the values along the given 1 by 1 square. 

%I currently have an error using 'integral2' because the function will not
%take the integral value of the distance between P1 and P2. 
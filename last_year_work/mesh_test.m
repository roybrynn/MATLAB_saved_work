    [Xvals,Yvals] = meshgrid(0:0.01:1);
    
    diff = Xvals - Yvals;
    diff(2:100, 2:100) = 0; %only values along the boundary

    disp(norm(diff))

    n_diff = norm(diff);

    xmax = 0.75;
    xmin = 0.25;
    ymax = 0.75;
    ymin = 0.25;

    g = ((n_diff*(xmax) - n_diff*(xmin))*ymax) - ((n_diff*(xmax) - n_diff*(xmin))*ymin);

    %integral(@(x)integral(@(y)(fun, 0.25, 0.75, 'ArrayValued', 1)), 0.25, 0.75, 'ArrayValued', 1);
    
    %all_values = integral2(fun, 0.25, 0.75, 0.25, 0.75);
    disp(g)
    %gb = all_values(0:100);
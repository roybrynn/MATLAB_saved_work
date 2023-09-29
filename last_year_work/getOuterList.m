function [final] = getOuterList(g1)
    %meshgrid used to take the x and y values for the grid
    [X, Y] = meshgrid(0:0.01:1);

    %reshape used to transform the values to lists
    %I chose to keep the length as 101*101 becuase my function to create g
    %will always return a matrix of that size
    gvec = reshape(g1, 101*101, 1);
    xvec = reshape(X, 101*101, 1);
    yvec = reshape(Y, 101*101, 1);

    %index is used to grab the appropriate indeces based on the boundary
    %conditions
    index = find((xvec == 0) | (xvec == 1) | (yvec == 1) | (yvec == 0));
    final = gvec(index);
end
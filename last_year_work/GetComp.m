function [comb] = GetComp(xmin,xmax,ymin,ymax)

%constructs the 4 seperate outer rectangles
%uses GetG for the function values
pt1 = GetG(0,xmin,0,1);
pt2 = GetG(xmin,xmax,0,ymin);
pt3 = GetG(xmin,xmax,ymax,1);
pt4 = GetG(xmax,1,0,1);

%combines the matrices to one
comb = pt1 + pt2 + pt3 + pt4;
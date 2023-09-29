%Test 1, grabs the values for the boundary rectangle
% this only seems to work when the lists are covering the same size
% rectangle

%Grabs the matrices for our given b0 and g

%tests for b = (0.72,0.78,0.72,0.78), (0.65,0.8,0.65,0.8)
%work as suspected, they produce disp1, and disp2
%if they are the same size, no matter where on the grid,
%they produce the same result (maybe need to alter equation?)
b0 = GetG(5,8,6,8);
b = GetG(6,8,5,8);

%structure from https://www.mathworks.com/matlabcentral/answers/275988-how-could-find-multiples-of-a-given-number-n
%gives the multiples in order of smallest to largest
siz = 1;
range = 1:siz;
vals = range(mod(siz,range)==0);
%using these values & center of mass we have drastically narrowed down the
%% 
%available options of what b0 could be

%center of mass
com = [7.5,7.5];

%length of vals
len = length(vals);
%list to hold the b0s
%https://www.mathworks.com/matlabcentral/answers/599413-how-do-you-create-a-matrix-of-matrices
b0s = cell(round(len/2), 2);
%function that will list all the options for b0 given the multiples and
%center of mass
[b0s, bool] = b0Options(b0,vals,com,b0s);

%if statements with corresponding print statements
%none should print
if(b0 >= b)
    disp("All of b0 is greater than or equal to b, b0 is within or equal to b")
else
    disp("other")
    %disp(b0 ~= b)
end
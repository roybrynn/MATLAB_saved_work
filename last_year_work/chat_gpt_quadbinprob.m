% Define the quadratic objective function
H = [2 -1; -1 4];
f = [-4; -8];
% Define the linear inequality constraints
A = [-3 -4; -1 -2; 1 0; 0 1];
b = [-24; -6; 0; 0];
% Define the variable bounds and type
lb = [0; 0];
ub = [Inf; Inf];
ctype = optimoptions('BB');
% Solve the integer quadratic programming problem
[x, fval] = quadprog(H, f, A, b, [], [], lb, ub, [], ctype);
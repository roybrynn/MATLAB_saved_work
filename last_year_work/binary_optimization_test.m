%build our b0, our total domain is going to be 0-10
b0 = GetG(2,5,7,10);

%area & center of mass
area = 9;
com = [3.5,8.5];
upperbound = [1,6,6,11];

%found from: https://www.mathworks.com/help/optim/ug/office-assignments-binary-integer-programming-solver-based.html
%also used chatgpt :)
A = b0;
b = ones(length(b0),1);

%what is f and how do i use it to work with the minimization
%c = trapz(, trapz(,, 2));
% sum(sijaijx - data)^2
c = [5, 4, 3,2 4, 53,24, 3,4,53];
f = -c;

lb = zeros(size(f)); %lower bound of 0
ub = lb + 1; %upper bound of 1
intvars = 1:length(f);

[x,fval] = intlinprog(f,intvars,A,b,[],[],lb,ub);

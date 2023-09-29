function [y] = f2(b)
%Double integral version of F
z = 0:1/100:1; %chi range
x = 0.25:(b - 0.25)/100:b; %range of integral
b = cast((b*100), "int64");

z(1:25) = 0;%initializing chi
z(25:(b)) = 1;
z(((b) + 1):end) = 0;

y = trapz(x, z); %given function value for b
end
%b0 = Integral2D(0.65, 0.67, 0.12, 0.21);

%int over B abs(x-y) chiB(y)dy
%Computation Domain = omega = [0,1] x [0,1]
%Ground Truth B = a shape within 1 by 1 square
%x = 0:0.01:1;
%y = 0:0.01:1;


%go on the boundRY
%store g(x1 - x400) as a vector
g = Integral2D(0.25, 0.75, 0.25, 0.75); %basic g value
disp(g);
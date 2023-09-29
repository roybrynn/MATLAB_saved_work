%Represents the case where b0 is a square/rectangle

%build our b0, our total domain is going to be 0-10
b0 = GetG(3,5,2,4);

area = 3;

%center of mass
com = [4,3];

%structure from https://www.mathworks.com/matlabcentral
%/answers/275988-how-could-find-multiples-of-a-given-number-n
%list of area multiples
range = 1:area;
vals = range(mod(area,range)==0);
%length of vals
len = length(vals);

%https://www.mathworks.com/matlabcentral/answers/5994
%13-how-do-you-create-a-matrix-of-matrices
%list to hold the b0s
b0s = cell(round(len/2), 2);

%will run until a match/upper bound is found for b0
match = false;

while match == false
    [b0s, bool] = b0Options(b0,vals,com,b0s);
    if bool == true
        match = true;
    else
        area = area + 1;
        range = 1:area;
        vals = range(mod(area,range)==0);
        len = length(vals);
        b0s = cell(round(len/2), 2);
    end
end

%https://www.mathworks.com/help/matlab/ref/cell2mat.html
%xdata = cell2mat([coms(1,:)]); ydata = cell2mat([coms(2,:)]);
%https://www.mathworks.com/matlabcentral/answers/336002-
%how-to-plot-a-single-collection-of-points
%plot(xdata,ydata,'-*')
%disp(count)
%disp(test_com)



function [b0s, bool] = b0Options(b0, vals, com, b0s)

%length of vals
len = size(vals);
bool = false;

%for loop through the list (divided by 2 because they are pairs)
for x = 1:(round(len(2)/2))
    %the first b0 is the smaller multiple as the x length, and the larger
    %as the y length
    %disp(x)
    a = com(1) - (0.5*vals(x));
    b = com(1) + (0.5*vals(x));
    c = com(2) - (0.5*(vals(len(2)-x+1)));
    d = com(2) + (0.5*(vals(len(2)-x+1)));

    b0s{x, 1} = GetG(a, b, c(1), d(1));

    if(b0 >= b0s{x,1})
        disp(x + "x1" + " " + a + " " + b + " " + c(1) + " " + d(1))
        bool = true;
    end
    %the second is the reverse
    b0s{x, 2} = GetG(c(1), d(1), a, b);
    
    if(b0 >= b0s{x,2})
        disp(x + "x2" + " " + c(1) + " " + d(1) + " " + a + " " + b)
        bool = true;
    end
end



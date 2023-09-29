%com plotting
hold on
plot(1.75,6,'b*')

%plotting for the lines of a plain setup
%need to setup lines different from plotting 
%EX: (xstart, xstop, ystart, ystop)
x = [1,3];
y = 3;
a = line(x,[y,y]);
ax = linspace(1,3);
ay = linspace(3,3);

x = [1,2];
y = 7;
b = line(x,[y,y]);
bx = linspace(1,2);
by = linspace(7,7);

x = [2,3];
y = 5;
c = line(x,[y,y]);
cx = linspace(2,3);
cy = linspace(5,5);

x = 1;
y = [3,7];
d = line([x,x],y);
dx = linspace(1,1);
dy = linspace(3,7);

x = 3;
y = [3,5];
e = line([x,x],y);
ex = linspace(3,3);
ey = linspace(3,5);

x = 2;
y = [5,7];
f = line([x,x],y);
fx = linspace(2,2);
fy = linspace(5,7);

lines = cell(6,1);
lines{1,1} = ax;
lines{1,2} = ay;
lines{2,1} = bx;
lines{2,2} = by;
lines{3,1} = cx;
lines{3,2} = cy;
lines{4,1} = dx;
lines{4,2} = dy;
lines{5,1} = ex;
lines{5,2} = ey;
lines{6,1} = fx;
lines{6,2} = fy;

%center of mass
com = [1.75,6];

%loop through each line's individual points, connect to com and see if
%there is overlap
%for x = 1:6
    x_vals = lines{1, 1};
    y_vals = lines{1, 2};
    %disp(line_testh.XData(end)-line_test.XData(start))
    for y = 1:(100)
        xh = [com(1),x_vals(y)];
        yh = [com(2),y_vals(y)];
        line_x = linspace(com(1), x_vals(y)); 
        line_y = linspace(com(2), y_vals(y));
        count = 0;
        for l = 2:(6)
            hold1 = lines{l, 1};
            hold2 = lines{l, 2};
            [x1,y1] = polyxpoly(hold1, hold2, line_x, line_y);
            if(isempty(x1))
                count = count + 1;
            end
            if(count >= 2)
                plot(x1, y1, 'r*');
                disp("Line: " + l + " vals: " + x1 + " " + y1)
            end
        end
    end
%end

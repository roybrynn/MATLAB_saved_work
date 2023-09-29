    x = 0:0.01:1;
    y = 0:0.01:1;
    g = 0:0.01:4;
    count = 1;
    vals = 1:1:4;
    vals(1:4) = 0;
while (count < 400)
        count2 = 1;
        while(count2 < 3)
            if ((count <= 100) || ((count >= 200) && (count <= 300)))
               if (count == 300)
                   x1 = 1;
                   y1 = y(100-(count - 201));
               elseif(count >= 200)
                   x1 = 1;
                   y1 = y(100-(count - 200));
               else
                   x1 = 0;
                   y1 = y(count);
               end
            else
               if ((count > 100) && (count < 200))
                   x1 = x(count-100);
                   y1 = 1;
               else
                   x1 = x(100-(count - 301));
                   y1 = 0;
               end
            end
       vals(count2) = x1;
       vals(count2 + 2) = y1;
       if (count2 == 1)
            count = count + count2;
       end
       count2 = count2 + 1;
        end
        disp(vals);
end

%Program that initializes G to be the 0.01 differences along the 1x1 square
%boundary and their corresponding x and y values
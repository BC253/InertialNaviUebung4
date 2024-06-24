function [C] = C(number,w)

if number == 1
    C = [1 0 0
        0 cos(w) sin(w)
        0 -sin(w) cos(w)];
elseif number == 2
    C = [cos(w) 0 -sin(w)
        0 1 0
        sin(w) 0 cos(w)];
elseif number == 3
    C = [cos(w) sin(w) 0 
        -sin(w) cos(w) 0
        0 0 1];
else 
    disp('wrong number,only 1 2 3')
end

   


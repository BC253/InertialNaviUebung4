function y = Heun(func,y0,x0,t)
y_P = y0 + t*func(y0,x0);
y = y0 + 1/2*t*(func(y0,x0)+func(y_P,x0));
y(7:10) = compact(normalize(quaternion(y(7:10)')));
end

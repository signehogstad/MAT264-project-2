function f_int = midpoint_integral_f(x0,x1,y0,y1,f)
nx = 15;
ny = 15;
hx = (x1-x0)/nx;
hy = (y1-y0)/ny;
f_int = 0;
for i = 0:(nx-1)
    for j = 0:(ny-1)
        xi = x0 + hx/2 + i*hx;
        yj = y0 + hy/2 + j*hy;
        f_int = f_int + hx*hy*f(xi,yj);
    end
end
end
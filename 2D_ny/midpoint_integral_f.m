function int_val = midpoint_integral_f(x0,x1,y0,y1,f)
A = (x1-x0)*(y1-y0);
x_avg = (x1+x0)/2;
y_avg = (y1+y0)/2;
int_val = f(x_avg,y_avg)*A;
end
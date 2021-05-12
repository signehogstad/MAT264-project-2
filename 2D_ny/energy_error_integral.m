function e_e_i_val = energy_error_integral(x0,x1,y0,y1,rx0,rx1,ry0,ry1,a)

k = 1;
dx = x1-x0; dy = y1-y0;
a1 = a(2); a2 = a(3); a3 = a(4);

% Integral 1:
t1 = 1/3*dx*dy*rx0^2;
t2 = 2*rx1*rx0*dy*(0.5*x1^3-x1^2*x0-1/3*x1^3+0.5*x1^2*x0)/dx^2;
t3 = -2*dy*rx0*rx1*(0.5*x0^2*(x0+x1)-x1*x0^2-1/3*x0^3)/dx^2;
t4 = 1/3*dx*dy*rx1^2;

t5 = 1/3*dy*ry0^2*dx;
t6 = 2*ry1*ry0*dx*(0.5*y1^2*(y1+y0)-y1^2*y0-1/3*y1^3)/dy^2;
t7 = -2*ry0*ry1*dx*(0.5*y0^2*(y1+y0)-y0^2*y1-1/3*y0^3)/dy^2;
t8 = 1/3*ry1^2*dx*dy;

integral_1 = t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8;

% Integral 2:
d1 = dx*(rx1+rx0)*(a1*dy+0.5*a3*(y1^2-y0^2));
d2 = dy*(ry1+ry0)*(a2*dx+0.5*a3*(x1^2-x0^2));

integral_2 = d1 + d2;

% Integral 3

s1 = a1^2*dx*dy+dx*a3*a1*(y1^2-y0^2);
s2 = 1/3*a3^2*dx*(y1^3-y0^3);
s3 = a2^2*dx*dy + a3*a2*dy*(x1^2-x0^2);
s4 = 1/3*a3^2*dy*(x1^3-x0^3);

integral_3 = s1 + s2 + s3 + s4;

e_e_i_val = integral_1 + integral_2 + integral_3;
end
function [answer] = stjernenorm_flux(q,x,N)
% Star-norm, flux error part
answer = 0;
k = 1;
for i = 1:N
    dx = x(i+1)-x(i);
    b1 = (q(i+1)-q(i))/dx;
    b0 = q(i);
    xi = (x(i+1) + x(i))/2;

    term1 = 0.25*k*pi*(2*pi*x(i+1) + sin(2*pi*x(i+1)));
    term2 = -0.25*k*pi*(2*pi*x(i) + sin(2*pi*x(i)));
    term3 = 2*sin(pi*x(i+1))*(b0+b1*(x(i+1)-xi));
    term4 = -2*sin(pi*x(i))*(b0+b1*(x(i)-xi));
    term5 = 2*b1/pi*(cos(pi*x(i+1))-cos(pi*x(i)));
    term6 = 1/k * b0^2*dx;
    term7 = 1/k*b0*b1*((x(i+1)-xi)^2-(x(i)-xi)^2);
    term8 = 1/3*1/k*b1^2*((x(i+1)-xi)^3 - (x(i)-xi)^3);

   
    answer = answer + (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8);

end
end
function cons_int = energy_conservation_integral(q,xh,N)
% Conservation integral

cons_int = 0;

for i = 1:N
    
    dx = xh(i+1)-xh(i);
    b_i1 = (q(i+1)-q(i))/dx;
   
    t1 = dx*((pi^4)/2 + b_i1^2);
    t2 = ((pi^3)/4)*(sin(2*pi*xh(i)) - sin(2*pi*xh(i+1)));
    t3 = 2*pi*b_i1*cos(pi*xh(i+1));
    t4 = 2*pi*b_i1*cos(pi*xh(i));
    
    cons_int = cons_int + (t1 + t2 + t3 - t4);
end
end
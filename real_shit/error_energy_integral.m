function val = error_energy_integral(xh,uh,q,N,k)

x = linspace(0,1,N);

for i = 1:N
   b0 = q(i);
   b1 = (q(i+1)-q(i))/(xh(i+1)-xh(i));
   b1 = 6;
   a1 = (uh(i+1)-uh(i))/(xh(i+1)-xh(i));
   
   x_ip = xh(i+1);
   x_i = xh(i);
   

%    x_ig = 0.5*(x_ip+x_i);
%    x_ig = x_i;
   x_ig = x(i);
   
   t1 = b0^2*(x_ip - x_i);
   t2 = b0*b1*(x_ip - x_ig)^2 - b0*b1*(x_i-x_ig)^2;
   t3 = 2*b0*a1*k*(x_ip - x_i);
   t4 = 1/3*b1^2*(x_ip - x_ig)^3 - 1/3*b1^2*(x_i - x_ig)^3;
   t5 = k*b1*a1*((x_ip-x_ig)^2-(x_i-x_ig)^2);
   t6 = k^2*a1^2*(x_ip - x_i);

   
   error_energy_vector(i) = t1 + t2 + t3 + t4 + t5 + t6;
   
end
val = sum(error_energy_vector);

end
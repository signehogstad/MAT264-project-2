function val = error_energy_integral(xh,uh,q,N,k)

for i = 1:N
   b0 = q(i); 
   b1 = (q(i+1)-q(i))/(xh(i+1)-xh(i));
   a1 = (uh(i+1)-uh(i))/(xh(i+1)-xh(i));
   
   x_ip = xh(i+1);
   x_i = xh(i);
   
   x_ig = 0.5*(x_ip+x_i);
   
   term1 = b0^2*(x_ip - x_i);
   term2 = b0*b1*(x_ip - x_ig)^2 - b0*b1*(x_i-x_ig)^2;
   term3 = 2*b0*a1*k*(x_ip-x_i);
   term4 = 1/3*b1^2*(x_ip-x_ig)^3 - 1/3*b1^2*(x_i-x_ig)^3;
   term5 = k^2*a1^2*(x_ip-x_i);
   
   error_energy_vector(i) = term1 + term2 + term3 + term4 + term5; 
end
val = sum(error_energy_vector);
end
function error = real_error(uh,xh,N)

error = 0;

for i = 1:N
    
    a_i1 = (uh(i+1)-uh(i))/(xh(i+1)-xh(i));
    
    dx = xh(i+1) - xh(i); 
    
    error = error + 0.25*(2*(2*a_i1^2 + pi^2)*dx + 8*a_i1*(sin(pi*xh(i)) ...
        - sin(pi*xh(i+1))) + pi*(sin(2*pi*xh(i+1)) - sin(2*pi*xh(i))));
    
end

end

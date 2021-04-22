function q = flux(u,dx,k)
% Calculating and plotting q
% k = 1 is a kind constant. In heat transfer situations, higher k will mean
% higher transport of heat or something.

for i = 2:N
   q(i) = -k * (u(i) - u(i-1))/dx;
end

q(1) = q(2) + pi*cos(pi*dx*0.5)-pi;

%NBNBNB JUKS
q(N+1) = q(1);

end

function q = flux(u,dx,k,X,N)
% Calculating and plotting q
% k = 1 is a kind constant. In heat transfer situations, higher k will mean
% higher transport of heat or something.

for i = 2:N
   q(i) = -k * (u(i) - u(i-1))/dx;
end

q(1) = q(2) + pi*cos(pi*dx*0.5)-pi;

%NBNBNB IKKJE JUKS LENGRE ME E PRO YAY
q(N+1) = q(N) + pi*(cos(pi*X)-cos(pi*(X + 0.5*dx)));

end

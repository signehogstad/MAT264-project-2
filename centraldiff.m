function [x,u_est,dx] = centraldiff(X,N)

f = @(x) pi*pi*sin(pi*x);
u_eksakt = @(x) sin(pi*x);
    
x_val = linspace(0,X,N);
dx = x_val(2)-x_val(1);

% x_ny for mean value of edge-potential
count = 0;
for i = 2:N
    x_ny(i) = dx*(0.5 + count);
    count = count + 1;
end
x_ny(1) = 0;
x_ny(N+1) = x_ny(N)+0.5*dx;

x = x_ny;

A = -1;
B = 2;
A = diag(A*ones(1,N-3),1) + diag(B*ones(1,N-2)) + diag(A*ones(1,N-3),-1);

b = zeros(N,1);
X = zeros(N,1);

for i = 1:N
    b(i) = dx^2 * f(x_val(i));
end

b(1) = 0;
b(N) = 0;

u(2:N-1) = A\b(2:N-1); % Matrisen er kun gyldig for ikke-randpunktene
u(1) = 0;
u(N) = 0;

% Creates a v-vector with average values of u between cells
for i = 2:N
    v(i) = (u(i-1)+u(i))/2;
end
v(1) = 0;
v(N+1) = 0;

u_est = v;

plot(x_ny,v,'r')
hold on
% Calculating and plotting q
% k = 1 is a kind constant. In heat transfer situations, higher k will mean
% higher transport of heat or something.
k = 1; 
for i = 2:N
   q(i) = -k * (u(i) - u(i-1))/dx;
end

q(1) = q(2) + pi*cos(pi*dx*0.5)-pi;

%NBNBNB JUKS KANSKJE
q(N+1) = q(1);

plot(x_ny,q,'b')
hold on

plot(x_val,u_eksakt(x_val))

end
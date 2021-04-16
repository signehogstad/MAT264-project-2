% clear;
% clc;
hold off;

f = @(x) pi*pi*sin(pi*x);
u_eksakt = @(x) sin(pi*x);

X = 2;
N = 10;
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
% Calculating the "energy error interval"
error_grad_greier = 0;

for i = 1:N
    
   b_i0 = q(i+1); 
   b_i1 = (q(i+1)-q(i))/dx;
   a_i1 = (v(i+1)-v(i))/dx;
   if i == 1
       dX = 0.5*dx;
   elseif i == N 
       dX = 0.5*dx;    
   else 
       dX = dx;
   end  
   error_grad_greier(i) = dX*(1/k)*((b_i0 + k*a_i1)^2 + ... 
   + b_i1*dX*(b_i0 + k*a_i1) * (dX^2)/3 * b_i1^2);
end



% Conservation integral

cons_int = 0;

for i = 2:N
    
    b_i1 = (q(i+1)-q(i))/dx;
    
    cons_int = cons_int + dx*((pi^4)/2 + b_i1^2) + ((pi^3)/4)*(sin(2*pi*x_ny(i)) - ...
        sin(2*pi*x_ny(i+1))) + 2*pi*b_i1*(cos(pi*x_ny(i+1)) - cos(pi*x_ny(i)));
    
end

cons_int

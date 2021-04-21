clear;
clc;
hold off;
format long;

f = @(x) pi*pi*sin(pi*x);
u_eksakt = @(x) sin(pi*x);

X = 2;
N = 800;
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
    
   b_i0 = q(i); 
   b_i1 = (q(i+1)-q(i))/(x_ny(i+1)-x_ny(i));
   a_i1 = (v(i+1)-v(i))/(x_ny(i+1)-x_ny(i));
   
   b_0v(i) = b_i0;
   b_1v(i) = b_i1;
   a_1v(i) = a_i1;
   
   x_ip = x_ny(i+1);
   x_i = x_ny(i);
   
   x_ig = 0.5*(x_ip+x_i);
   
   term1 = b_i0^2*(x_ip - x_i);
   term2 = b_i0*b_i1*(x_ip - x_ig)^2 - b_i0*b_i1*(x_i-x_ig)^2;
   term3 = 2*b_i0*a_i1*k*(x_ip-x_i);
   term4 = 1/3*b_i1^2*(x_ip-x_ig)^3 - 1/3*b_i1^2*(x_i-x_ig)^3;
   term5 = k^2*a_i1^2*(x_ip-x_i);
   
   error_grad_greier(i) = term1 + term2 + term3 + term4 + term5; 
   
end

energy_error = sum(error_grad_greier)


% Conservation integral

cons_int = 0;

for i = 1:N
    
    b_i1 = (q(i+1)-q(i))/(x_ny(i+1) - x_ny(i));
    
    dX = x_ny(i+1)-x_ny(i);
   
    cons_int = cons_int + (dX/2)*(dX*((pi^4)/2 + b_i1^2) + ((pi^3)/4)*(sin(2*pi*x_ny(i)) - ...
        sin(2*pi*x_ny(i+1))) + 2*pi*b_i1*(cos(pi*x_ny(i+1)) - cos(pi*x_ny(i))));
    
end

cons_int

majorant = sqrt(energy_error) + sqrt(cons_int)

% Testing testing

test = 0;

for i = 1:N
    
    a_i1 = (v(i+1)-v(i))/(x_ny(i+1)-x_ny(i));
    
    dX = x_ny(i+1) - x_ny(i); 
    
    test = test + 0.25*(2*(2*a_i1^2 + pi^2)*dX + 8*a_i1*(sin(pi*x_ny(i)) ...
        - sin(pi*x_ny(i+1))) + pi*(sin(2*pi*x_ny(i+1)) - sin(2*pi*x_ny(i))));
    
end

error = sqrt(test)





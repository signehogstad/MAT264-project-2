clear;
clc;
hold off;

f = @(x) pi*pi*sin(pi*x);
u_eksakt = @(x) sin(pi*x);

X = 2;
N = 10;
x_val = linspace(0,X,N);
dx = x_val(2)-x_val(1);

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

plot(x_val,u,'m')
hold on

% plot(x_val,u_eksakt(x_val),'k')
% hold on

% Creates a v-vector with average values og u between cells
for i = 2:N-1
    v(i) = (u(i-1)+u(i))/2;
end
v(1) = 0;
v(N) = 0;


% Idk if this is overthinking but here is an x-value that looks like the
% one in Jan's slides ([0 0.5 1.5 2.5 ......])
count = 0;
for i = 2:N
    x_ny(i) = dx*(0.5 + count);
    count = count + 1;
end
x_ny(1) = 0;
% plot(x_ny,v,'r')

% Calculating and plotting q
k = 1; % #Yolo-value of k bc we don't know what it is
for i = 2:N
   q(i) = -k * (u(i) - u(i-1))/dx;
end



plot(x_ny,q,'b')



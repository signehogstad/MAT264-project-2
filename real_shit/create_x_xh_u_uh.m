function [x,xh,u,uh,dx] = create_x_xh_u_uh(X,N,f)
% Matrix diagonal values
a = -1;
bb = 2;

% x-vector for "whole" cells
x = linspace(0,X,N);
dx = x(2) - x(1);
count = 0;

% xh-vector for half cells
for i = 2:N
    xh(i) = dx*(0.5 + count);
    count = count + 1;
end
xh(1) = 0;
xh(N+1) = xh(N)+0.5*dx;

% Matrix A
A = diag(a*ones(1,N-3),1) + diag(bb*ones(1,N-2)) + diag(a*ones(1,N-3),-1);

b = zeros(N,1);

% Create b-vector. Used for calculating u later
for i = 1:N
    b(i) = dx^2 * f(x(i));
end

b(1) = 0;
b(N) = 0;

% Potential for "whole" cells, u:
u(2:N-1) = A\b(2:N-1);
u(1) = 0;
u(N) = 0;

% Potential evaluated in the half points, uh:
for i = 2:N
    uh(i) = (u(i-1)+u(i))/2;
end
uh(1) = 0;
uh(N+1) = 0; 
end

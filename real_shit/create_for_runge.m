function [x,u] = create_for_runge(X,N,f)
% Matrix diagonal values
a = -1;
bb = 2;

% x-vector
x = linspace(0,X,N);
dx = x(2) - x(1);

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

end
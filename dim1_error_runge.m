
X = 2;
N1 = 200;
N2 = 20;
[x1,u1,dx1] = centraldiff(X,N1);
[x2,u2,dx2] = centraldiff(X,N2);

k = dx2/dx1;

alpha = 2;

error_bound = 0;


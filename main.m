clear;
clc;
hold off

k = 1;

f = @(x) pi*pi*sin(pi*x);
u_eksakt = @(x) sin(pi*x);

N = 40;
X = 2;

% Getting stuff from create_stuff
[x,xh,u,uh,dx] = create_x_xh_u_uh(X,N,f);

% Calc. flux. Function looks like this q = flux(u,dx,k,X)
q = flux(u,dx,k,X,N);

% Calc. error_energy_interval. Function looks like this val = error_energy_integral(xh,q,N)
energy_int = sqrt(error_energy_integral(xh,uh,q,N,k));

% Calc. error_conservation_integral. Function looks like this
cons_int = sqrt(error_conservation_integral(q,xh,N));

% Calculating majorant
majorant = energy_int + cons_int

% Calculating real error for comparison
error = sqrt(real_error(uh,xh,N))


plot(xh,uh,'m','LineWidth',1.2)
hold on

plot(xh,q,'g','LineWidth',1.2)
hold on

plot(x,u_eksakt(x),'k','LineWidth',1.2,'LineStyle','--')





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

% Runge error, possibly terribly wrong rip.
% NB: Only works for Nh1 = n * Nh2, where n is an integer.
Nh1 = 200; Nh2 = 50;

[xh1,uh1] = create_for_runge(X,Nh1,f);
[xh2,uh2] = create_for_runge(X,Nh2,f);

[runge_error_bound, uh2_new] = Runge_error(uh1,uh2);
runge_error_bound
real_error_L2 = norm(u_eksakt(xh1) - uh1,2)

% Make a longer xh2 (xh2_new) to match uh2_new
xh2_new = make_longer(xh2,xh1);

% Plotting xh2/uh2 and xh2_new/uh2_new for comparison just for lolz
plot(xh2_new,uh2_new,'b','LineWidth',1.5)
hold on

plot(xh2,uh2,'k','LineStyle','--','LineWidth',1.5)
hold on

plot(xh1,uh2_new,'m','LineWidth',1.5)



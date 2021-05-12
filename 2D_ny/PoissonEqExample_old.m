% Poisson equation: -Div(k(x,y)Grad(u)) = f
hold off;
clear; clc;
format short

%Diffusivity
% k = @(x,y) 1+(x-1).*(y-1);
k = @(x,y) ones(length(x),length(y));

% Exact solution
u_fabricated = @(x,y) sin(2*pi*x).*sin(2*pi*y);
% u_fabricated = @(x,y) x.^2.*y.^2-x.^2.*y-x.*y.^2+x.*y;

%RHS function
f = @(x,y) 8*pi*pi*sin(2*pi*x).*sin(2*pi*y)-(2*pi*cos(2*x*pi).*sin(2*pi*y).*(y-1)+2*pi*sin(2*pi*x).*cos(2*pi*y).*(x-1));
% f = @(x,y) -2*y^2+2*y-2*x^2+2*x;
% fi =  @(x,y) (8*pi*pi*sin(2*pi*x).*sin(2*pi*y).*k(x,y)-(2*pi*cos(2*x*pi).*sin(2*pi*y).*(y-1)+2*pi*sin(2*pi*x).*cos(2*pi*y).*(x-1)))^2;

%Determine number of cells in each direction
nx = 40;
ny = 40;
dx = 1/(nx-1);
dy = 1/(ny-1);

% Build matrices
[A, b, G, D, K, cells, edges] = assembleMatrices(nx,ny,f,k);

% Solve linear system
u = A\b;

% Reshape and plot
[X,Y] = meshgrid(0:dy:1,0:dx:1);
U = reshape(u,nx,ny);

% Calculate fabricated solution vector and compute error in "L2-norm"
u_fabricated_vect = u_fabricated(cells(:,1),cells(:,2));
error7 = norm(u-u_fabricated_vect,2)*sqrt(dx*dy)

% Team Ingrid and Signe sitt verk herfra og nedover vvvvv

q = -K.*G*u;
u = reshape(u,nx,ny);

[avg_pot,U_int,a_mat] = avg_pot_int_pot(nx,ny,U,X,Y);


error_L2 = 0;
for i = 1:nx
    for j = 1:ny
        error_L2 = error_L2 + (U_int(i,j) - U(i,j))^2;
    end
end
error = sqrt(dx*dy*error_L2)



figure(2)
surf(X,Y,U_int)
axis([0 1 0 1])
xlabel('x')
ylabel('y')
title('Interpolated potential yay')

colorbar()

conservation_integral = 0;
real_e_e = 0;
energy_error_int = 0;
stjernenorm_midt = 0;
j = 1;

for i = 1:nx*ny

    [coord_vals,flux_vals,j,a_vec,maj_int] = get_flux_and_coord(cells,edges,q,a_mat,i,f,j);
        
    a1 = a_vec(2); a2 = a_vec(3); a3 = a_vec(4);
 
    energy_error_int = energy_error_int + maj_int;
    
    xv = coord_vals(1,1); xe = coord_vals(1,2);
    ys = coord_vals(2,1); yn = coord_vals(2,2);
    
    dy = yn-ys; dx = xe-xv;

    rv = flux_vals(1,1); 
    re = flux_vals(1,2);
    rs = flux_vals(2,1); 
    rn = flux_vals(2,2);
    
    fi = @(x,y) (8*pi*pi*sin(2*pi*x).*sin(2*pi*y)-(2*pi*cos(2*x*pi).*sin(2*pi*y).*(y-1)+2*pi*sin(2*pi*x).*cos(2*pi*y).*(x-1)) - ((re-rv)/dx + (rn-rs)/dy)).^2;

    conservation_integral = conservation_integral + midpoint_integral_f(xv,xe,ys,yn,fi);
      
    f_jevlig = @(x,y) (2*pi*cos(2*pi*x)*sin(2*pi*y))^2-4*pi*cos(2*pi*x)*sin(2*pi*y)*(a1+a3*y)+(a1+a3*y)^2 + ...
        + (2*pi*cos(2*pi*y)*sin(2*pi*x))^2-4*pi*cos(2*pi*y)*sin(2*pi*x)*(a2 + a3*x) + (a2+a3*x)^2;
    
    real_e_e = real_e_e + midpoint_integral_f(xv,xe,ys,yn,f_jevlig);
    
    f_merjevlig = @(x,y) (2*pi*cos(2*pi*x)*sin(2*pi*y))^2+4*pi*cos(2*pi*x)*sin(2*pi*y)*((xe-x)*rv/dx+(re*(x-xv))/dx)+(((xe-x)*rv/dx+(re*(x-xv))/dx))^2 + ...
        + (2*pi*cos(2*pi*y)*sin(2*pi*x))^2+4*pi*cos(2*pi*y)*sin(2*pi*x)*((yn-y)*rs/dy+(rn*(y-ys))/dy)+(((yn-y)*rs/dy+(rn*(y-ys))/dy))^2;
    
    stjernenorm_midt = stjernenorm_midt + midpoint_integral_f(xv,xe,ys,yn,f_merjevlig);
    
end
% energy_error_int = sqrt(energy_error_int);
% 
% conservation_integral = 1/(pi*sqrt(2))*sqrt(conservation_integral);
% 
% real_e_e = sqrt(real_e_e)
% 
% majorant = energy_error_int + conservation_integral
% 
% stjernenorm_midt = sqrt(stjernenorm_midt);
% 
% stjernenorm = real_e_e + conservation_integral + stjernenorm_midt
% 
% efficiency_index = 3*majorant/stjernenorm




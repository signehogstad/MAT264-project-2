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
f = @(x,y) 8*pi*pi*sin(2*pi*x).*sin(2*pi*y).*k(x,y)-(2*pi*cos(2*x*pi).*sin(2*pi*y).*(y-1)+2*pi*sin(2*pi*x).*cos(2*pi*y).*(x-1));
% f = @(x,y) -2*y^2+2*y-2*x^2+2*x;
fi =  @(x,y) (8*pi*pi*sin(2*pi*x).*sin(2*pi*y).*k(x,y)-(2*pi*cos(2*x*pi).*sin(2*pi*y).*(y-1)+2*pi*sin(2*pi*x).*cos(2*pi*y).*(x-1)))^2;

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
figure(1)
set(gcf,'Position',[100 100 1200 500])
subplot(1,2,1);
surf(X,Y,U)
axis([0 1 0 1])


% Calculate fabricated solution vector and compute error in "L2-norm"
u_fabricated_vect = u_fabricated(cells(:,1),cells(:,2));
error = norm(u-u_fabricated_vect,2)*sqrt(dx*dy)
xlabel('x')
ylabel('y')
title('u')
colorbar()


% Compute the flux vector [F_1,F_2] at the cells and plot with quiver
[F_1, F_2] = flux(G,K,u,cells,nx,ny);
subplot(1,2,2)
quiver(cells(:,1),cells(:,2),F_1,F_2)
axis([0 1 0 1])
xlabel('x')
ylabel('y')
title('Flux')

% Plot the fabricated solution
U_fabricated = reshape(u_fabricated_vect,nx,ny);
figure(2)
subplot(1,1,1)
surf(X,Y,U_fabricated)
axis([0 1 0 1])
xlabel('x')
ylabel('y')
title('u_{fabricated}')
colorbar()

% Team Ingrid and Signe sitt verk herfra og nedover vvvvv

q = -K.*G*u;

[avg_pot,int_pot_mat,a] = avg_pot_int_pot(nx,ny,U,X,Y);

cells;
edges;

surf(X,Y,int_pot_mat)
axis([0 1 0 1])
xlabel('x')
ylabel('y')
title('Interpolated potential yay')

colorbar()

maj_int = 0;
int_val = 0;
for i = 1:nx*ny
    [coord_vals,flux_vals] = get_flux_and_coord(cells,edges,q,a,i,f);
    
    x0 = coord_vals(1,1); x1 = coord_vals(1,2);
    y0 = coord_vals(2,1); y1 = coord_vals(2,2);
    
    %IDK? F?r for syk konvergens av dette. Men tror det blir feil ? ikke ha
    %det med, s???? idk....
    dy = y1-y0; dx = x1-x0;
    
    rx0 = flux_vals(1,1); rx1 = flux_vals(1,2);
    ry0 = flux_vals(2,1); ry1 = flux_vals(2,2);
    
    maj_int = maj_int + energy_error_integral(x0,x1,y0,y1,rx0,rx1,ry0,ry1,a);
    
    fi = @(x,y) (8*pi*pi*sin(2*pi*x).*sin(2*pi*y).*k(x,y)-(2*pi*cos(2*x*pi).*sin(2*pi*y).*(y-1)+2*pi*sin(2*pi*x).*cos(2*pi*y).*(x-1)) - ((rx1-rx0)/dx + (ry1-ry0)/dy))^2;
    
    int_val = int_val + midpoint_integral_f(x0,x1,y0,y1,fi);
end
maj_int = sqrt(maj_int);
int_val = 1/(pi*sqrt(2))*sqrt(int_val)

majorant = maj_int + int_val;
colormap winter

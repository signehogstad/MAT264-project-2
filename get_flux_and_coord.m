function [coord_vals,flux_vals,maj_int,j,a_vec] = get_flux_and_coord(cells,edges,q,a_mat,i,f,j)


a_vec = zeros(4,1);
w = 6; n = 5; e = 4; s = 3;

% Cell_info looks like: [x y south east north west]
cell_info = cells(i,:);
W = cell_info(w); N = cell_info(n); E = cell_info(e); S = cell_info(s);
% Getting x0, y0, x1, y1 by fetching edges(s,:), edges(e,:) and so on.
% Erlend har vært supernice and given us all these matriser. edges(s,:)
% will be information abouth the southern end of the cell. It will
% contain a "center"-value as x-coordinate, and will be y0 as
% y-coordinate.

% Corners:
if W == 0 && S == 0
    x0 = 0; x1v = edges(E,:); x1 = x1v(1); 
    y0 = 0; y1v = edges(N,:); y1 = y1v(2); 
    rx0 = 0; rx1 = 0;
    ry0 = 0; ry1 = 0;
elseif S == 0 && E == 0
    x0v = edges(W,:); x0 = x0v(1); x1 = 1;
    y0 = 0; y1v = edges(N,:); y1 = y1v(2);
    rx0 = 0; rx1 = 0;
    ry0 = 0; ry1 = 0;
elseif E == 0 && N == 0
    x0v = edges(W,:); x0 = x0v(1); x1 = 1;
    y0v = edges(S,:); y0 = y0v(2); y1 = 1;
    rx0 = 0; rx1 = 0;
    ry0 = 0; ry1 = 0;
elseif W == 0 && N == 0
    x0 = 0; x1v = edges(E,:); x1 = x1v(1);
    y0v = edges(S,:); y0 = y0v(2); y1 = 1;
    rx0 = 0; rx1 = 0;
    ry0 = 0; ry1 = 0;

% All edges:

% WEST
elseif W == 0
    x0 = 0; x1v = edges(E,:); x1 = x1v(1); 
    y0v = edges(S,:); y0 = y0v(2); y1v = edges(N,:); y1 = y1v(2); 
    
    dy = y1-y0;
    
    rx1 = q(E); rx0 = rx1 - 1/dy*midpoint_integral_f(x0,x1,y0,y1,f); 
    ry0 = 0; ry1 = 0;
    
% SOUTH
elseif S == 0
    x0v = edges(W,:); x0 = x0v(1); x1v = edges(E,:); x1 = x1v(1);
    y0 = 0; y1v = edges(N,:); y1 = y1v(2);
    
    dx = x1-x0;
    
    rx0 = 0; rx1 = 0;
    ry1 = q(N); ry0 = ry1 - 1/dx*midpoint_integral_f(x0,x1,y0,y1,f); 
    
% EAST
elseif E == 0
    x0v = edges(W,:); x0 = x0v(1); x1 = 1;
    y0v = edges(S,:); y0 = y0v(2); y1v = edges(N,:); y1 = y1v(2);
    
    dy = y1-y0;
    
    rx0 = q(W); rx1 = rx0 + 1/dy*midpoint_integral_f(x0,x1,y0,y1,f);
    ry0 = 0; ry1 = 0;
    
% NORTH
elseif N == 0
    x0v = edges(W,:); x0 = x0v(1); x1v = edges(E,:); x1 = x1v(1);
    y0v = edges(S,:); y0 = y0v(2); y1 = 1;
    
    dx = x1-x0;
    
    rx0 = 0; rx1 = 0;
    ry0 = q(S); ry1 = ry0 + 1/dx*midpoint_integral_f(x0,x1,y0,y1,f);

% All other cells:
else
    x0v = edges(W,:); x0 = x0v(1); x1v = edges(E,:); x1 = x1v(1);
    y0v = edges(S,:); y0 = y0v(2); y1v = edges(N,:); y1 = y1v(2);
    rx0 = q(W); rx1 = q(E);
    ry0 = q(S); ry1 = q(N);
    
    a_vec = a_mat(j,:);
    j = j + 1;
    
end

% rx0
% rx1

coord_vals = [x0 x1; y0 y1;];
flux_vals = [rx0 rx1; ry0 ry1;];

maj_int = energy_error_integral(x0,x1,y0,y1,rx0,rx1,ry0,ry1,a_vec);

a1 = a_vec(2); a2 = a_vec(3); a3 = a_vec(4);


end

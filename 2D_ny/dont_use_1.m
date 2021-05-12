function [avg_pot,a_mat,interpolated_potential] = interpolate_potential(nx,ny,U,cells,edges)

for i = 1:nx-1
    for j = 1:ny-1
        avg_pot(i,j) = (U(i,j)+U(i+1,j)+U(i+1,j+1)+U(i,j+1))/4;
    end
end
count = 0;
for i = 2:nx-1
    for j = 2:nx-1
        
        count = count + 1;
        
        edge_x1 = edges((nx-1)+j,:);
        edge_x2 = edges(nx+j,:);
        edge_y1 = edges(nx*ny-1+i,:);
        edge_y2 = edges(nx*ny+i,:);
        
        x1 = edge_x1(1);
        x2 = edge_x2(1);
        y1 = edge_y1(2);
        y2 = edge_y2(2);
        
        Q11 = avg_pot(i-1,j-1); 
        Q12 = avg_pot(i-1,j);
        Q21 = avg_pot(i,j-1);
        Q22 = avg_pot(i,j);
        
        b = [Q11; Q21; Q12; Q22];

    
          
        A = [1 x1 y1 x1*y1; 1 x1 y2 x1*y2; 1 x2 y1 x2*y1; 1 x2 y2 x2*y2];
        
        a = A\b;
        a_mat(count,:) = a;
        a0 = a(1); a1 = a(2); a2 = a(3); a3 = a(4);
        
        xa = (x1+x2)/2;
        ya = (y1+y2)/2;
        
        interpolated_potential(i-1,j-1) = a0 + a1*xa + a2*ya + a3*xa*ya;
    
    end
end
interpolated_potential = conv2(interpolated_potential,[0,0,0;0,1,0;0,0,0]);

end
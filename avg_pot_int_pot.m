function [avg_pot,int_pot_mat,a_mat] = avg_pot_int_pot(nx,ny,U,X,Y)

dx = 1/(nx-1); dy = 1/(ny-1); k = 1;
for i = 1:nx-1
    for j = 1:ny-1
        avg_pot(i,j) = (U(i,j)+U(i+1,j)+U(i+1,j+1)+U(i,j+1))/4;
    end
end
cell_number = 1;
for i = 2:nx-1
    for j = 2:ny-1       
        [int_pot,a] = v_interpolation_value(i,j,avg_pot,X(1,:),Y(:,1),dx,dy);        
        est_pot = U(i,j);
        int_pot_mat(i-1,j-1) = int_pot;
        int_pot;
        a_mat(cell_number,:) = a;
        cell_number = cell_number + 1;
        a_mat;
    end
end
avg_pot;
int_pot_mat = conv2(int_pot_mat,[0,0,0;0,1,0;0,0,0]);
end



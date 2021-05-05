function int_pot = v_interpolation_value(i,j,avg_pot,x,y,dx,dy)

x1 = x(i)-dx/2; x2 = x(i+1)-dx/2; y1 = y(j)-dy/2; y2 = y(j+1)-dy/2;
A = [1 x1 y1 x1*y1; 1 x1 y2 x1*y2; 1 x2 y1 x2*y1; 1 x2 y2 x2*y2];
        
Q11 = avg_pot(i-1,j-1);
Q12 = avg_pot(i-1,j);
Q21 = avg_pot(i,j-1);
Q22 = avg_pot(i,j);
b = [Q11 Q12 Q21 Q22]';

a = A\b;

term1 = a(1);
term2 = a(2)*x(i)+a(3)*y(j);
term3 = a(4)*x(i)*y(j);
int_pot = term1 + term2 + term3;


end

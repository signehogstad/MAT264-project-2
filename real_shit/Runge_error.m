function [error,uh2_new] = Runge_error(uh1,uh2,N,n,dx)
% Takes in uh1 and uh2, which are two different approximations of u, with
% different stepsize h (h1 and h2). Calculates the runge error bound.
k = N/n;
k;
error = 0;
if mod(length(uh1),length(uh2)) == 0
    uh2_new = make_longer(uh2,uh1);
    for i = 1:length(uh2_new)
       error = error + (uh1(i)-uh2_new(i))^2;  
    end
    error = sqrt(dx*error);
    error = 1/(k^2-1) * error;
else
    disp('NBNBNB: uh1 must be n times the size of uh2, where n is an integer #Tryagain')
end

end




function [error_bound,uh2_new] = Runge_error(uh1,uh2)
% Takes in uh1 and uh2, which are two different approximations of u, with
% different stepsize h (h1 and h2). Calculates the runge error bound.
if mod(length(uh1),length(uh2)) == 0
    uh2_new = make_longer(uh2,uh1);
    error_bound = 1/3*norm(uh1 - uh2_new,2);
else
        disp('NBNBNB: uh1 must be n times the size of uh2, where n is an integer #Tryagain')
end

end




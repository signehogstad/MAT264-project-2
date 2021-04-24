function vectn = make_longer(vect,targvect)
% This function takes in two vectors. One vect, which is short, and another
% one, targvect, which is an integer (n) number of times longer than vect 
% #goodinglish.
% It returns the vector vectn, which is a "longification" of vect, such as
% the length of vectn and targvect is the same. vectn is kind of equal to
% vect, but consists of n times the same value instead of only one. Ex:
% vect = [1 2 3 4], n = 2, vectn = [1 1 2 2 3 3 4 4 5 5].
vectn = []; 
j = 0;
n = length(targvect)/length(vect);
while length(vectn) < length(targvect)
    for i = (1 + j*n):(1 + n*j + n)
        vectn(i) = vect(j + 1);
    end
    j = j + 1;
    vectn = vectn(1:end-1);
end



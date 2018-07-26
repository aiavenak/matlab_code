function [x fx] = f_convolve_gaussKernel(lamda,sigma)

% convolve the eigenvalues of the normalized Laplacian of a Graph
% with a gaussian distribution kernel
%
% input     lamda - eigevalues
%           sigma - std of gaussian den. finc. 
%                   smaller sigma emphasize the finer details 
%                   whereas larger values bring out the global patterns
%
% Reference  Benerjee PhD thesis 2008

[d1 d2] = size(lamda);
if d1 > d2
    lamda = lamda';
end

x = [0:.001:2]';
nx = length(x);
lamda = lamda(ones(nx,1),:);

var = sigma^2;
C =(2*pi*var)^(-1/2);

fx = bsxfun(@(x,lamda)(abs(x-lamda)).^2 ,x,lamda);
fx = exp(fx./(-2*var));
fx = C.*sum(fx,2);

fx = fx./sum(fx);
function [pr bins] = f_eigvals_distr(aij,nbins)

[d1 d2] = size(aij);
if (d1 == d2)
    % aij is adyacency matrix
    vaij = f_get_spectrum_lap(aij);
else  % aij is vector of eigenvalues
    vaij = aij;
end

if nargin < 2
    nbins = sqrt(length(vaij));
end
bins = linspace(0,2,nbins);

pr = hist(vaij,bins);
pr = pr./length(vaij);


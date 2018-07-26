function [pr hh bins] = f_get_distr(X,rmv_zeros,bins)

if rmv_zeros
    ff = X ~= 0;
    X = X(ff);
end

if nargin < 3
    [d1 d2] = size(X);
    nbins = ceil(sqrt(max([d1 d2])));    
    maxX = max(X(:));
    minX = min(X(:));
    bins = linspace(minX,maxX,nbins);
end

hh = hist(X(:),bins);
pr = hh./numel(X);


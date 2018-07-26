function [vn v] = fcn_vi(x,y)
%FCN_VI     variation of information
%
%   V = FCN_VI(X,Y) returns the variation of information between two
%       partitions, X and Y. Variation of information is comptued as
%       follows:
%
%           V(X,Y) = 2*H(X,Y) - H(X) - H(Y)
%
%       Here, H(X,Y) is the joint entropy of X and Y. H(X) and H(Y) are the
%       marginal entropies of X and Y, respectively. To normalize V(X,Y),
%       we divide by log(N) where N is the length of the vectors X and Y,
%       respectively.
%
%       Inputs:     X,      partition membership vector
%                   Y,      partition membership vector
%
%       Outputs:   VN,      normalized variation of information
%                   V,      variation of information
%
%       Note: The output V has units bits.
%
%   Richard Betzel, Indiana University, 2012
%

%modification history
%04.14.2012 - original version

nx = length(x);
ny = length(y);

if nx ~= ny
    error('input vector lengths do not match');
end

ux = unique(x);
uy = unique(y);

px = hist(x,ux)./nx;
py = hist(y,uy)./ny;

hx = -sum(px.*log(px));
hy = -sum(py.*log(py));

pxy = zeros(length(ux),length(uy));
for i = 1:length(ux);
    ind = x == ux(i);
    yyy = y(ind);
    pxy(i,:) = hist(yyy,uy);
end
pxy = nonzeros(pxy./(nx));
hxy = -sum(pxy.*log(pxy));

v = 2*hxy - hx - hy;
vn = v./log(nx);
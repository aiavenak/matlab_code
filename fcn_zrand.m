function zr = fcn_zrand(x,y)
% FCN_ZRAND         z-score rand index
%
%   Reference: Traud et al (2010). Comparing community structure to
%   characteristics in online collegiate social networks. arxiv:0809.0690v3
%



% clear all
% close all
% clc

% two partitions
% x = [1 1 1 2 2 2 3 3 3 4 4 5]';
% y = [1 1 2 2 3 3 4 4 5 5 6 6]';

% build agreement matrices
indx = dummyvar(x);
indy = dummyvar(y);

% get size of partitions
n = length(x);

% build consistency matrix
nxy = indx'*indy;
vxy = nxy(nxy >= 2);
w11 = 0;
for i = 1:length(vxy)
    w11 = w11 + vxy(i)*(vxy(i) - 1)/2;
end

% row and column sum of consistency matrix
nx = sum(nxy,2);
ny = sum(nxy);
M1 = sum(nx.*(nx - 1)/2);
M2 = sum(ny.*(ny - 1)/2);
M = n*(n - 1)/2;

c1 = n*(n.^2 - 3*n - 2) - 8*(n + 1)*M1 + 4*sum(nx.^3);
c2 = n*(n.^2 - 3*n - 2) - 8*(n + 1)*M2 + 4*sum(ny.^3);

a = M/16;
b = (4*M1 - 2*M).^2;
c = (4*M2 - 2*M).^2;
d = c1*c2/(16*n*(n - 1)*(n - 2));
e = (b - 4*c1 - 4*M)*(c - 4*c2 - 4*M)./(64*n*(n - 1)*(n - 2)*(n - 3));

s = a - (b*c)/(256*(M^2)) + d + e;
zr = (sqrt(s)^-1)*(w11 - (M1*M2/M));
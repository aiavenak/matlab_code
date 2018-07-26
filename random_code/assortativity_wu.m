function   r = assortativity_wu(CIJ)
%ASSORTATIVITY      Assortativity coefficient
%
%   r = assortativity_wu(CIJ);
%
%   The assortativity coefficient is a correlation coefficient between the
%   strength of all nodes on two opposite ends of a link. A positive
%   assortativity coefficient indicates that nodes tend to link to other
%   nodes with the same or similar degree.
%
%   Inputs:     CIJ,    weighted undirected connection matrix
%
%   Outputs:    r,      assortativity coefficient
%
%   Notes: the function computes the weighted assortativity coefficient 
%          described in Rubinov and Sporns (2010) NeuroImage.
%
%   Reference:  Leung st al. (2007) Physica A 378:591-602
%
clear all
clc
N = 20;
path2code = fullfile(pwd);
cd C:\Users\aiavenak\Projects\Rewire_Brain\data\Hg;
load('DSI','CIJ');
CIJ = CIJ(1:N,1:N);
cd(path2code);

CIJ = (CIJ - min(CIJ(:))) / max(CIJ(:));
stre = sum(CIJ);
[i,j] = find(triu(CIJ,1)>0);
K = length(i);
strei = stre(i);
strej = stre(j);
wij = zeros(1,K);
for count = 1:1:K
    wij(count) = CIJ(i(count),j(count)); 
end

% compute assortativity
r = ( sum(wij.*strei.*strej)/K - (sum(0.5*wij.*(strei+strej))/K)^2 ) / ...
    ( sum(0.5*wij.*(strei.^2+strej.^2))/K - (sum(0.5*wij.*(strei+strej))/K)^2 );


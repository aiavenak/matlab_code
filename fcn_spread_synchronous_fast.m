function [Y,Z] = fcn_spread_synchronous_fast(A,seeds,thr)
%FCN_SPREAD_SYNCHRONOUS_FAST    synchronous update spreading model
%
%   Y = FCN_SPREAD_SYNCHRONOUS_FAST(A,SEEDS,THR) runs the linear threshold
%       model on network A, starting with seeds specified in the SEEDS
%       array. A node adopts when a fraction of its neighbors, specified by
%       the input, THR, are also active.
%
%   Inputs:
%
%       A,      connectivity matrix; if directed, sums across rows/columns
%               are out/in degree, respectively.
%
%       SEEDS,  array of seeds; each row contains two columns, the first is
%               the seed node index and the second the seed node type,
%               indicated by a positive integer.
%
%       THR,    the threshold required for adoption - can be a scalar
%               applied uniformly to all nodes, or a vector applied to to
%               specific nodes.
%
%   Outputs:
%
%       Y,      each column is a time step and each row is a node; thus,
%               cells indicate the state of nodes at each time step.
%
%       Z,      N x 2 matrix
%               1st column: adoption times
%               2nd column: adoption colours
%
%
%   Richard Betzel, Indiana University, 2014

n = length(A);
k = sum(A,2);
P = bsxfun(@rdivide,A,k);
nseeds = size(seeds,1);

x = zeros(n,nseeds);
ind = (seeds(:,2) - 1)*n + seeds(:,1);
x(ind) = 1;

count = 0;
xupdate = [];
Y = zeros(n,25);

% Z     first column: time adopted
%       second column: color adopted
Z = inf(n,2);
Z(seeds(:,1),2) = seeds(:,2);
Z(seeds(:,1),1) = ones(nseeds,1);


vec = (1:nseeds)';
stopflag = true;
while stopflag
    %keyboard;
    count = count + 1;
    Y(:,count) = x*vec;
    Z(xupdate,1) = count;
    Z(xupdate,2) = x(xupdate,:)*vec;
    stopflag = false;
    hasadopted = Y(:,count) > 0;
    xprime = P*x;
    [xmax,indmax] = max(xprime,[],2);
    xupdate = xmax > thr & ~hasadopted;
    %keyboard;
    indupdate = (indmax(xupdate) - 1)*n + find(xupdate);
    x(indupdate) = 1;
    if sum(xupdate)
        stopflag = true;
    end
end
%keyboard;
Y = Y(:,1:count);
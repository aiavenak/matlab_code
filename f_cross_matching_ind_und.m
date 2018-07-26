function M = f_cross_matching_ind_und(CIJ1,CIJ2)
%MATCHING_IND_UND       matching index
%
%   M0 = CROSS_MATCHING_IND_UND(CIJ) computes matching index between 
%   node i in CIJ1 and all other nodes in CIJ2.
%   Matching index is a measure of similarity between two nodes'
%   connectivity profiles in two different matrices.
%
%   Inputs:     CIJ1, CIJ2    undirected adjacency matrices, must be the
%                             same size
%
%   Outputs:    M0,           cross-matching index matrix.
%
%   Richard Betzel, Indiana University, 2013
%   Modified by Andrea Avena-Koenigsberger, 2014
%

D1 = size(CIJ1);
D2 = size(CIJ2);
if sum(D1-D2) ~=0
    fprintf('matrix size does not match \n');
    return;
end

N = D1(1);
I = ~eye(N);
M = zeros(N,N);

for i = 1:N
    
    c1 = CIJ1(i,:);
    use = bsxfun(@or,c1,CIJ2);
    %use(:,i) = 0;
    %use = use.*I;
    
    ncon1 = bsxfun(@times,use,c1);
    ncon2 = bsxfun(@times,use,CIJ2);
    ncon = sum(ncon1 + ncon2,2);
    
    M(:,i) = 2*sum(ncon1 & ncon2,2)./ncon;
    
end

% M = M.*I;
M(isnan(M)) = 0;
%M0 = zeros(size(CIJ0));
%M0(R,R) = M;

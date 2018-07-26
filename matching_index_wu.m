function [Mall] = matching_index_wu(CIJ,i,j)
%MATCHING_IND       Matching index
%
%   [Mall] = matching_ind(CIJ);
%
%   For any two nodes u and v, the matching index computes the amount of
%   overlap in the connection patterns of u and v. Self-connections and
%   u-v connections are ignored. The matching index is a symmetric 
%   quantity, similar to a correlation or a dot product.
%
%   Input:      CIJ,    weighted undirected connection/adjacency matrix
%
%   Output:     Min,    matching index for incoming connections
%               Mout,   matching index for outgoing connections
%               Mall,   matching index for all connections
%
%   Notes:
%       Does not use self- or cross connections for comparison.
%       Does not use connections that are not present in BOTH u and v.
%       All output matrices are calculated for upper triangular only.
%
%
% Olaf Sporns, Indiana University, 2002/2007/2008/2013

N = size(CIJ,1);

if nargin==1

    % compare all connections
    Mall = zeros(N,N);
    for i=1:N-1
        c1 = CIJ(:,i);
        for j=i+1:N
            %c1 = [CIJ(:,i)' CIJ(i,:)];
            %c2 = [CIJ(:,j)' CIJ(j,:)];
            c2 = CIJ(:,j);
            use = ~(~c1&~c2);
            use(i) = 0;
            use(j) = 0;
            ncon = sum(c1(use))+sum(c2(use));
            if (ncon==0)
                Mall(i,j) = 0; %%% 0
            else
                c1u = c1(use); c2u = c2(use);
                ov = c1(use)&c2(use);
                Mall(i,j) = sum(c1u(ov)+c2u(ov))/ncon;
            end;
        end;
    end;

    % make full matrix
    Mall = Mall+Mall';

elseif nargin==3
    c1 = CIJ(:,i);
    c2 = CIJ(:,j);
    use = ~(~c1&~c2);
    use(i) = 0;
    use(j) = 0;
    ncon = sum(c1(use))+sum(c2(use));
    if (ncon==0)
        Mall = 0;
    else
        c1u = c1(use); c2u = c2(use);
        ov = c1(use)&c2(use);
        Mall = sum(c1u(ov)+c2u(ov))/ncon;
    end;
    
else
    error('invalid number of nargin!');
end


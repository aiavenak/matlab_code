function [ssd] = f_simm(G1,G2)

% note: G1 or G2 must not have nodes with zero degree

N = size(G1,1);
G1inv = full(inv(sparse(G1)));
G2inv = full(inv(sparse(G2)));

ssd = trace(G1*G2inv)+trace(G1inv*G2)-2*N;
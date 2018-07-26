function [eig_vals,eig_vect] = f_get_spectrum_lap(aij)

conn = find(sum(aij));
aij = aij(conn,conn);

[~,L]=f_lap(aij);
[eig_vect,eig_vals]=eig(L);

eig_vals = sort(diag(eig_vals));





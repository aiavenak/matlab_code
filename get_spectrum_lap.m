function eig_vals = get_spectrum_lap(aij)

conn = find(sum(aij));
aij = aij(conn,conn);

D = diag(sum(aij));
L = D - aij;

[kk eig_vals]=eig(L);

eig_vals = sort(diag(eig_vals));





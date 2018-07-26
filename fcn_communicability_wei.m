function C = fcn_communicability_wei(A)

D = (diag(sum(A,1)))^(-1/2);
C = expm(D*A*D);
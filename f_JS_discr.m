function JS = f_JS_discr(aij,bij,nbins)

[d1 d2] = size(aij);
[b1 b2] = size(bij);
if (d1 == d2) && (b1 == b2)
    % aij, bij are adyacency matrices
    % get laplacian and eigenvalues
    vaij = f_get_spectrum_lap(aij);
    vbij = f_get_spectrum_lap(bij);
else  % aij bij are vectors of eigenvalues
    vaij = aij;
    vbij = bij;
end

ul = max([max(vaij),max(vbij)]);
if nargin < 3
    nbins = mean(length(vaij),length(vbij));
    nbins = sqrt(nbins);
end
bins = linspace(0,ul,nbins);

pr_aij = hist(vaij,bins);
pr_aij = pr_aij./length(vaij);
pr_bij = hist(vbij,bins);
pr_bij = pr_bij./length(vbij);

pr = (pr_aij + pr_bij)./2;

KLaij = f_kullback_leibler(pr_aij,pr);
KLbij = f_kullback_leibler(pr_bij,pr);

JS = ((KLaij/2) + (KLbij/2))^(1/2);

function kl = f_kullback_leibler(p1,p2)

kl = p1.*log(p1./p2);
kl(isnan(kl)) = 0;
kl(isinf(kl)) = 0;
kl = sum(kl);
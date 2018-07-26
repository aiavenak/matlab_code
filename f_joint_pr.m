function [jpr,bins_i,bins_j] = f_joint_pr(M1,M2,nbins,lim1,lim2)

if nargin < 3
    bins_i = linspace(lim1(1),lim1(2),nbins(1)+1);
    bins_j = linspace(lim2(1),lim2(2),nbins(2)+1);
else
    deltaM1 = min(nonzeros(M1)).*0.01;
    deltaM2 = min(nonzeros(M2)).*0.01;
    bins_i = linspace(min(M1(:)),max(M1(:))+deltaM1,nbins(1)+1);
    bins_j = linspace(min(M2(:)),max(M2(:))+deltaM2,nbins(2)+1);
end

jpr = zeros(nbins(1),nbins(2));

for ii = 1:nbins(1)-1
    
   ind = M1 >= bins_i(ii) & M1 < bins_i(ii+1);
   M2ind = M2(ind);
   jpr(ii,:) = histcounts(M2ind,bins_j)./length(M2);
    
end

jpr(isnan(jpr)) = 0;


% imagesc(jpr,[0 .1]); colorbar;
% sum(jpr,1)
% sum(jpr,2)
% sum(sum(jpr))


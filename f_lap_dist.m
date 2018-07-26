function Dab = f_lap_dist(aij,bij)

% calculate the edit distance between two graphs A and B
% reference: Thomas Thorne and Michael P. H. Stumpf
%            evolution Graph spectral analysis of 
%            protein interaction network
%            J. R. Soc. Interface 2 May 2012

[d1 d2] = size(aij);
[b1 b2] = size(bij);
if (d1 == b1) && (d2 == b2)
    if d1 == d2  % aij, bij are matrices
        % get laplacian and eigenvalues
        [kk Lna] = f_lap(aij);
        [kk Lnb] = f_lap(bij);
        
        [kk Va] = eig(Lna);
        [kk Vb] = eig(Lnb);
        
        Va = diag(Va);
        Vb = diag(Vb);
    else  % aij, bij are vectors of eigenvalues
        Va = aij;
        Vb = bij;
    end
    
    Va = sort(Va);
    Vb = sort(Vb);
    
    Dab = sum((Va-Vb).^2);
    
else
    warning('sizes do not match');
    Dab = 0/0;
end
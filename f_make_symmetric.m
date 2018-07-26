function Asym = f_make_symmetric(A,flag_avr)

[dim1,dim2,dim3] = size(A);
mask = eye(dim1)>0;

mask_tril = find(bsxfun(@(X,Y) X.*Y, ones(dim1,dim2,dim3),tril(ones(dim1),-1)));
mask_triu = find(bsxfun(@(X,Y) X.*Y, ones(dim1,dim2,dim3),triu(ones(dim1),1)));

if all(isnan(A(mask_tril))) || all(isinf(A(mask_tril)))
    A(mask_tril) = 0;
elseif all(isnan(A(mask_triu))) || all(isinf(A(mask_triu)))
    A(mask_triu) = 0;
end
       
if dim3 == 1
    Asym = (A + A');
    Asym(mask) = A(mask);
else
    Asym = zeros([dim1,dim2,dim3]);
    for d = 1:dim3        
        Asym(:,:,d) = A(:,:,d) + (A(:,:,d))';
    end
    mask = (bsxfun(@(A1,A2) A1.*A2,ones([dim1,dim2,dim3]),mask))>0;
    Asym(mask) = A(mask);
end

if nargin == 2
    if flag_avr
        Asym = Asym./2;
    end
    Asym(mask) = A(mask);
end



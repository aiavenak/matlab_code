function J = f_jaccard_inx(M1,M2,edge_type)

if (nargin == 3) & (strcmp(edge_type,'wei'))
    
    intsct = sum(sum(min(M1,M2)));
    union = sum(sum(max(M1,M2)));
else
    
    intsct = sum(sum((M1>0) & (M2>0)));
    union  = sum(sum((M1>0) | (M2>0)));
end

J = intsct/union;


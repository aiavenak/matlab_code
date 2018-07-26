function fcn_mat2paj(aij,name)
%FCN_MAT2PAJ converts matlab matrix to pajek format
%
%   FCN_MAT2PAJ(AIJ,NAME) takes an [n x n] matrix AIJ of edge weights and
%   writes them to the file NAME in the same directory.
%

n = length(aij);
[x,y,wt] = find(triu(aij,1));
m = length(wt);

outname = sprintf('%s.net',name);
fid = fopen(outname,'w');
fprintf(fid,'*Vertices %i\n',n);
for i = 1:n
    fprintf(fid,'%i "%i"\n',i,i);
end

fprintf(fid,'*Edges %i\n',m);
for i = 1:m
    fprintf(fid,'%i %i %17.15f\n',x(i),y(i),wt(i));
end

fclose(fid);
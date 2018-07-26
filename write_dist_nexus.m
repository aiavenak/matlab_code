function write_dist_nexus(adj,filename)

nnodes = size(adj,1);

fid = fopen(sprintf('%s',filename),'wt');
fprintf(fid, '#nexus\n\n');
fprintf(fid, 'BEGIN taxa;\n');
fprintf(fid, 'DIMENSIONS ntax=%d;\n',nnodes);
fprintf(fid, 'TAXLABELS\n');

for i=1:nnodes
    fprintf(fid, '[%d] %d\n',i,i);
end
fprintf(fid, ';\n');
fprintf(fid, 'END [TAXA];\n\n');

fprintf(fid, 'BEGIN distances;\n');
fprintf(fid, 'DIMENSIONS ntax=%d;\n',nnodes);
fprintf(fid, 'FORMAT labels diagonal triangle=both;\n');
fprintf(fid, 'MATRIX\n');

for i=1:nnodes
    fprintf(fid, '[%d] %d\t',i,i);
    for j=1:nnodes
        fprintf(fid, '%f ',adj(i,j));
    end
    fprintf(fid, '\n');
end
fprintf(fid, ';\n');
fprintf(fid, 'END [distances];\n');

fclose(fid);        
end



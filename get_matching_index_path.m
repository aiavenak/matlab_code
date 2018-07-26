function MIsp = get_matching_index_path(adj,path)
MIsp = 0;

for i=1:length(path)-1
   adj(path(i),path(i+1))=0;
   adj(path(i+1),path(i))=0;
end
for i=1:length(path)
    for j=i+1:length(path)
        MIsp = MIsp + matching_index_wu(adj,path(i),path(j));
    end
end
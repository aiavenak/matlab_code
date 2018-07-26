function MIsp = get_shortest_path_matching_index_wei_und(adj,SPL,B)



N=size(adj,1);

MIsp = zeros(N,N);

for i=1:N-1
    
    for j=i+1:N
        

        path = retrieve_shortest_path(i,j,SPL,B);
        
        MIsp(i,j) = get_matching_index_path(adj,path);
        L = length(path);
        MIsp(i,j) = 2*MIsp(i,j)/(L*(L-1));
    end
end
MIsp = MIsp + MIsp';
function [ci_best q_best] = f_modularity_consensus(adj,q_reps,null_reps,path2results)

%addpath(genpath('/home/aiavenak/matlab/tools'));

n = size(adj,1);
if nargin == 4
    flag_save = 1;
else
    flag_save = 0;
end
% find modules

ci = zeros(n,q_reps);
for ii = 1:q_reps
    ci(:,ii) = modularity_louvain_und(adj);
end
if flag_save
    if exist(path2results,'file')
        save(fullfile(path2results),'ci','-append');
    else
        save(fullfile(path2results),'ci');
    end
end
% consensus

agr_mat = agreement(ci);                                            % bulid agreement matrix
agr_max = -inf;                                                         % set max null agreement to small value

for i = 1:null_reps                                                         % for each null model
    cinull = zeros(n,q_reps);                                            % initialize matrix for storing null modules
    
    for j = 1:q_reps                                                     % for each partition ...
        r = randperm(n);                                                % choose random order
        cinull(:,j) = ci(r,j);                                      % populate null matrix
    end
    
    agr_mat_null = agreement(cinull);                                   % build null agreement matrix
    agr_max = max(agr_max,max(agr_mat_null(:)));                        % find largest value in null agreement matrix
end

thr = max(agr_max);                                                     % threshold is the largest agreement value across all surrogate models

agr_mat_thr = agr_mat.*(agr_mat > thr);                                 % threshold the original agreement matrix

ci_new = zeros(n,q_reps);                                                % initialize matrix for storing reclustered communities

for i = 1:q_reps                                                         % for nreps ...
    ci_new(:,i) = modularity_louvain_und(agr_mat_thr);                  % recluster the communities
end

[ciu,~,mult] = fcn_unique_partitions(ci_new);                           % get unique partitions
[~,indmax] = max(mult);                                                 % find the partition with greatest multiplicity
ci_best = ciu(:,indmax);

%compute modularity
s = sum(adj(:));
m=max(ci_best);                                                         %number of modules
w=zeros(m);                                                             %new weighted matrix
for u=1:m
    for v=u:m
        wm=sum(sum(adj(ci_best==u,ci_best==v)));                        %pool weights of nodes in same module
        w(u,v)=wm;
        w(v,u)=wm;
    end
end
adj = w;
q_best = trace(adj)/s-sum(sum((adj/s)^2));                               %compute modularity

if flag_save
    save(fullfile(path2results),'ci_best','q_best','-append');
end



% 
% K = 8
% for rr = 1:20
%     load(fullfile(path2results,sprintf('SC_aal%d_%d',2^K,rr)),'ci_best','M_nf');
%     [X Y indsort] = grid_communities(ci_best);
%     figure; spy(M_nf(indsort,indsort))
%     hold on;
%     plot(X,Y,'k','linewidth',2)
%     shg
%     pause()  
% end






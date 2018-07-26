function [rand_mod] = random_modules_symmetric_nofeed(module,contig,coord)

% Input
% module = module assignments vector
% contig = contiguity matrix
% coords = centroid coordinates of parcellation, format = [x y z]
%
% Output
% rand_mod = randomized module assignments vector
%
% Frantisek Vasa
% April 2015

% checks
% check that length(contig) == size(module,1) == size(module,2)
if ~and(length(module) == size(contig,1),length(module) == size(contig,2))
    error('check dimensions of input')
end

% check that length(unique(module)) = max(module)
if length(unique(module)) ~= max(module)
    error('check module labeling') % or relabel modules here
end

% check that contiguity matrix is symmetric
if ~isequal(contig,contig')
    error('the contiguity matrix is asymmetric')
end

% check that nodes are not neighbours to themselves - if they are set diagonal = 0
if trace(contig) ~= 0
    contig(logical(eye(size(contig,1)))) = 0;
end

%%% store OVERALL initial size of modules for checking at the end
unique_mod = unique(module);
for i = 1:1:length(unique_mod)
    mod_size_input(i) = sum(module == unique_mod(i));
end
%%%

contig_in = contig; % store inputs
module_in = module;

% % check connectedness of contiguity matrix
% % if it is disconnected (eg: two hemispheres) run these separately
% if length(unique(get_components(contig))) > 2
%     error('the contiguity matrix has more than two connected components')
% elseif length(unique(get_components(contig))) == 2
%     disp('the contiguity matrix has two components - treating these as hemispheres')
% end

% identify larger component
for i = 1:2; lrg(i) = sum(get_components(contig) == i); end
[~,hi] = max(lrg); clear lrg;
hf = setdiff(1:2,hi);
ind_i = find(get_components(contig) == hi); % parcel IDs of 1st component (hemisphere)
ind_f = find(get_components(contig) == hf);

% run larger hemisphere first (then reflect + rearrange)
%for h = 1:2
        
hemi{hi} = (get_components(contig_in) == hi);
contig_hemi{hi} = contig_in(hemi{hi},hemi{hi}); % one hemisphere
mod_hemi{hi} = module_in(hemi{hi}); % input hemi 1

% prepare this for reflection at end
hemi{hf} = (get_components(contig_in) == hf);
contig_hemi{hf} = contig_in(hemi{hf},hemi{hf}); % one hemisphere
mod_hemi{hf} = module_in(hemi{hf}); % input hemi 1
        
module = mod_hemi{hi};
contig = contig_hemi{hi};

% check that each module forms a connected component and relabel split modules
%%% do this independently of hemisphere splitting?
mod_orig = module; % store input
c = max(module); % relabeling counter
mod_relabel = []; % keep track of which modules were relabeled to which
for i = 1:1:max(mod_orig)
    test_cont = contig(mod_orig == i,mod_orig == i);
    if length(unique(get_components(test_cont))) > 1
        ind = find(mod_orig == i); % find indices of nodes belonging to this module
        for j = 2:1:length(unique(get_components(test_cont))) % relabel components N = 2+
            c = c+1;
            module(ind(get_components(test_cont) == j)) = c;
            mod_relabel = [mod_relabel; [i c]];
        end
    end
end
clear ind

% store initial size of modules for checking at the end
unique_mod = unique(module);
for i = 1:1:length(unique_mod)
    mod_size(i) = sum(module == unique_mod(i));
end

% loop to ensure convergence
finished = false;
while ~finished
    
    % initialize output and 'working' variables
    rand_mod = zeros(size(module)); % random assignments (output)
    
    work_contig = contig; % "working" contiguity matrix (that gets reduced in loop)
    
    un = 1:1:length(module); % unassigned nodes
    un_mod = unique_mod; %1:1:max(module); % unassigned modules
    
    cnt = 0; % counter of assigned modules
    
    cnt_stuck1 = 0; % counter of stuck code 1 (lack of convergence)
    
    % while there are 2+ modules to be assigned
    while length(un_mod) > 1
        
        % evaluate degrees of contiguity matrix (to start assigning modules to
        % low-degree nodes, to prevent splitting the largest connected component)
        deg = sum(work_contig,2);
        
        % choose a node at random among lowest degree nodes
        % while loop to ensure convergence when lowest degree node breaks connected component
        while true    
            n = find(deg == min(deg)); % n = current node
            if length(n) > 1; n = n(randi(length(n))); end
            
            % check that assigning this node does not break connected component
            test_cont = work_contig; test_cont(n,:) = []; test_cont(:,n) = [];
            if length(unique(get_components(test_cont))) > 1
                deg(n) = nan; % maintain same number of elements in deg
                %disp('nan')
                continue
            elseif length(unique(get_components(test_cont))) == 1
                break
            end            
        end
        
        % choose a module at random (from unassigned modules)
        m = un_mod(randi(length(un_mod))); %%%%%
        
        cnt_stuck1 = cnt_stuck1+1; % ensure convergence
        if cnt_stuck1 > 5*length(unique_mod)
            %disp('not converging - reinitializing...');
            break
        end
        
        % number of nodes still to be assigned in current module (not counting initial node)
        n_remain = sum(module == m)-1; %sum(work_mod == m)-1;
        
        assign = n; % initialize (to-be-)assigned nodes
        curr_nei = find(work_contig(n,:)); % initialize neighbours of the initial node
        % could this be done within next while loop (to prevent doing this when n_remain = 0)?
        
        % while nodes remain in current module
        while n_remain > 0
            
            % choose a neighbor at random (new working node)
            n = curr_nei(randi(length(curr_nei)));
            
            % check that assigning this node does not break connected component
            test_cont = work_contig; test_cont([assign n],:) = []; test_cont(:,[assign n]) = [];
            if length(unique(get_components(test_cont))) > 1
                curr_nei = setdiff(curr_nei,n); % if it does remove this node from the neighbours
                if isempty(curr_nei); break; end % if neighbours are empty choose initial node again
                continue
            end
            
            % add node to list of nodes to be assigned to this module
            assign = [assign n];
            
            % update current neighbours with those of the newly assigned node
            curr_nei = [curr_nei find(work_contig(n,:))];
            curr_nei = unique(curr_nei); % remove duplicates
            curr_nei = setdiff(curr_nei,assign); % remove assigned nodes
            
            n_remain = n_remain-1;
            
        end
        
        if isempty(curr_nei); continue; end % if neighbours are empty choose initial node again
        
        % assign chosen nodes to current module
        rand_mod(un(assign)) = m;
        
        % remove assigned module from vector of unassigned modules
        un_mod = setdiff(un_mod,m); cnt = cnt+1;
        
        % remove assigned nodes from vector of unassigned nodes
        un(assign) = [];
        
        % update work_contig and work_mod
        work_contig(assign,:) = []; work_contig(:,assign) = [];
        %work_mod(assign) = [];
        
    end
    
    if length(un_mod) == 1
        
        % for the last module, assign all remaining unassigned nodes to it
        rand_mod(un) = un_mod; cnt = cnt+1;
        
        finished = true;
        
    end
    
end % while ~finished

% check that number of nodes in each module matches
for i = 1:1:length(unique_mod); rand_mod_size(i) = sum(rand_mod == unique_mod(i)); end
if ~isequal(mod_size,rand_mod_size); error('module size not preserved'); end

% check that each module forms a connected component
comp_num_check = [];
for i = 1:1:length(unique_mod)
    test_cont = contig(rand_mod == unique_mod(i),rand_mod == unique_mod(i));
    comp_num_check(i) = length(unique(get_components(test_cont)));
end
if ~isequal(comp_num_check,ones(1,length(unique_mod))); error('not all modules form connected components'); end

% relabel relabelled modules back
for i = 1:1:size(mod_relabel)
    rand_mod(rand_mod == mod_relabel(i,2)) = mod_relabel(i,1);
end

rand_out{hi} = rand_mod;

%end % loop over two hemishperes

%%% second hemisphere
% loop to ensure convergence
finished = false;
while ~finished
stuck = 0;

% label final (smaller) hemisphere based on initial (larger) hemisphere using coords
% for each left hemisphere parcel assign value of closest right parcel
ci = coord(ind_i,:); % initial hemisphere
cf = coord(ind_f,:); % final hemisphere

ni = length(ci);
nf = length(cf);

ci_r = ci; ci_r(:,1) = abs(ci_r(:,1));
cf_r = cf; cf_r(:,1) = abs(cf_r(:,1));

% for each lh coord, measure distance to all rh coords, find lowest
i_to_f = zeros(ni,nf);
for i = 1:1:ni
    i_to_f(i,:) = sqrt(sum((repmat(ci_r(i,:),nf,1)-cf_r).^2,2));
end

% find lowest 
[~,i_to_f_min] = min(i_to_f,[],1); % right
%[~,f_to_i_min] = min(i_to_f,[],2); % left

% reflect initial to final using i_to_f_min
rand_out{hf} = rand_out{hi}(i_to_f_min);

% initial number of members of each module on each side
mod_in_i = module_in(ind_i);
mod_in_f = module_in(ind_f);
% 1st col = initial / 2nd col = final
for i = 1:1:max(module_in)   
    mod_num_in(i,1) = sum(mod_in_i == i);
    mod_num_in(i,2) = sum(mod_in_f == i);
end

% evaluate (as)symmetry of input modules - deal with lateralized modules
% which modules have no members in init hemi but do in the final?
assym_to = find(and(mod_num_in(:,1) == 0,mod_num_in(:,2) > 0));

% which modules have no members in final hemi but do in the initial?
assym_from = find(and(mod_num_in(:,1) > 0,mod_num_in(:,2) == 0));

% flip some of these around
while ~or(isempty(assym_from),isempty(assym_to));
    temp_from = randi(length(assym_from)); % random assym_i module to relabel to
    temp_to = randi(length(assym_to)); % random assym_i module to relabel to
    rand_out{hf}(rand_out{hf} == assym_from(temp_from)) = assym_to(temp_to);
    assym_from(temp_from) = []; assym_to(temp_to) = []; 
    if or(isempty(assym_from),isempty(assym_to)); break; end
end
clear temp_from temp_to

% random number of members of each module on each side
mod_rand_i = module_in(ind_i);
mod_rand_f = module_in(ind_f);
% 1st col = initial / 2nd col = final
for i = 1:1:max(module_in)   
    mod_num_rand(i,1) = sum(rand_out{hi} == i);
    mod_num_rand(i,2) = sum(rand_out{hf} == i);
end

% for modules where numbers don't match, reassign parcels
% difference in numbers
mod_diff = mod_num_rand(:,2) - mod_num_in(:,2); % compare second columns (final hemispheres) of module numbers
mis = find(mod_diff ~= 0);

% if there are still modules in assym_from or assym_to, 
% reassign one of excess parcels to it
while ~isempty(assym_to)
    % random module from assym_to
    temp_to = randi(length(assym_to));
    % find modules with excess parcels, excluding ones which have no members in initial hemisphere
    temp_from = setdiff(find(mod_diff > 0),find(and(mod_num_rand(:,1) == 0,mod_num_rand(:,2) > 0)));
    % pick one of those modules
    if length(temp_from) > 1; temp_from = temp_from(randi(length(temp_from))); end
    % find a parcel from this module which has neighbours from other modules    
    temp = contig_hemi{hf};
    temp(rand_out{hf} ~= temp_from,:) = 0;
    temp(:,rand_out{hf} == temp_from) = 0; 
    [temp_parc,~] = find(temp); temp_parc = temp_parc(randi(length(temp_parc)));
    % relabel
    rand_out{hf}(temp_parc) = assym_to(temp_to);
    assym_to(temp_to) = []; 
end
clear temp_from temp_to temp_parc

%%% repeat for assym_from? shouldn't be an issue as nodes which have
%%% members in initial hemi but not in final should get relabeled in the
%%% final hemi to other modules

while ~isempty(mis)
        
    %%% ATTEMPTED "SYSTEMATIC" VERSION
    from_vect = find(mod_diff > 0); % modules with excess parcels
    to_vect = find(mod_diff < 0); % modules with lack of parcels
    
    from_to_mat = zeros(max(module_in));
    for i = 1:1:length(from_vect) % loop over modules with excess parcels
        for j = 1:1:length(to_vect) % loop over modules with lack of parcels
            % reduced contiguity matrix to identify neighbouring parcels
            % from "from" and "to" modules
            temp = contig_hemi{hf};
            temp(rand_out{hf} ~= from_vect(i),:) = 0;
            temp(:,rand_out{hf} ~= to_vect(j)) = 0;
            from_to_mat(from_vect(i),to_vect(j)) = sum(temp(:)); % identify cases where neighbours are present
        end
    end
    
    % if no more contiguous modules that need relabeling
    if sum(from_to_mat(:)) == 0; break; end
    
    % pick one of these entries
    [from_mod,to_mod] = find(from_to_mat);
    if length(from_mod) > 1;
        ind = randi(length(from_mod));
        from_mod = from_mod(ind); 
        to_mod = to_mod(ind); 
    end

%     %%% WORKING VERSION
%     % start with (one of) modules with excess parcels
%     from_mod = find(mod_diff > 0);
%     if length(from_mod) > 1; from_mod = from_mod(randi(length(from_mod))); end
%     
%     % redistribute to module
%     to_mod = find(mod_diff < 0);
%     if length(to_mod) > 1; to_mod = to_mod(randi(length(to_mod))); end
    
    % reduced contiguity matrix to identify neighbouring parcels
    % from "from" and "to" modules
    temp = contig_hemi{hf};
    temp(rand_out{hf} ~= from_mod,:) = 0;
    temp(:,rand_out{hf} ~= to_mod) = 0;
    
    if sum(temp(:)) == 0; continue; end % if these two modules have no neighbours
    
    %figure; imagesc(temp);
    
    [ind_from,ind_to] = find(temp);
    
    % pick a pair of nodes at random, relabel
    rel = randi(length(ind_from));
    
    rand_out{hf}(ind_from(rel)) = to_mod;
    %rand_out{hi}(ind_to(rel)) = from_mod;
    
    % update mod_diff and mis
    for i = 1:1:max(module_in)        
        mod_num_rand(i,2) = sum(rand_out{hf} == i);
        %mod_num_rand(i,1) = sum(rand_out{hi} == i);
    end
    
    % update mod_diff and mis
    mod_diff = mod_num_rand(:,2) - mod_num_in(:,2); % number of parcels to redistribute
    mis = find(mod_diff ~= 0); % modules where numbers don't match
    
    length(mis);
    
end

cnt = 0;
% deal with remaining modules that need relabeling
while ~isempty(mis)
    %disp('big loop')
    if cnt > 100;
        % find modules with excess parcels, excluding ones which have no members in initial hemisphere
        temp_from = setdiff(find(mod_diff > 0),find(and(mod_num_rand(:,1) == 0,mod_num_rand(:,2) > 0)));
        
        %%% extra condition - if no such cases
        if isempty(temp_from); break; end % disp('temp from empty'); 
        
        % pick one of those modules
        if length(temp_from) > 1; temp_from = temp_from(randi(length(temp_from))); end
        % find a parcel from this module which has neighbours from other modules
        temp = contig_hemi{hf};
        temp(rand_out{hf} ~= temp_from,:) = 0;
        temp(:,rand_out{hf} == temp_from) = 0;
        if sum(temp(:))==0; break; end % disp('temp empty'); 
        [temp_parc,~] = find(temp); temp_parc = temp_parc(randi(length(temp_parc)));
        % modules with lack of parcels
        to_mod = find(mod_diff < 0);
        if length(to_mod) > 1; to_mod = to_mod(randi(length(to_mod))); end
        
        % relabel
        rand_out{hf}(temp_parc) = to_mod;
        
        % update mod_diff and mis
        for i = 1:1:max(module_in)
            mod_num_rand(i,2) = sum(rand_out{hf} == i);
            %mod_num_rand(i,1) = sum(rand_out{hi} == i);
        end
        
        % update mod_diff and mis
        mod_diff = mod_num_rand(:,2) - mod_num_in(:,2); % number of parcels to redistribute
        mis = find(mod_diff ~= 0); % modules where numbers don't match
        %disp(['nonconvergence - relabel ' num2str(to_mod)])
        if isempty(mis); break; end % disp('break 1'); 
        cnt = 0;
        continue
    end % if cnt > 100
    
    % loop over all modules - for each module, find neighbour modules
    neigh_mod = {};
    for i = 1:1:max(module_in)
        % reduced contiguity matrix to identify neighbouring parcels
        % from "from" and "to" modules
        temp = contig_hemi{hf};
        temp(rand_out{hf} ~= i,:) = 0;
        temp(:,rand_out{hf} == i) = 0;
        [~,c] = find(temp); %%% STUCK HERE
        neigh_mod{i} = unique(rand_out{hf}(c));
    end
    
    %%% start with direct neighbours (recode from above?)
    
    % indirect neighbours - one step removed
    % pick a module which needs rewiring
    from_vect = find(mod_diff > 0); % modules with excess parcels
    to_mod = find(mod_diff < 0); % modules with lack of parcels
    if length(to_mod) > 1; to_mod = to_mod(randi(length(to_mod))); end
    
    % loop over its neigbouring modules, looking for modules with excess
    % parcels within their neighbours
    if isempty(neigh_mod{to_mod});
        % find modules with excess parcels, excluding ones which have no members in initial hemisphere
        temp_from = setdiff(find(mod_diff > 0),find(and(mod_num_rand(:,1) == 0,mod_num_rand(:,2) > 0)));
        
        %%% extra condition - if no such cases
        if isempty(temp_from); break; end % disp('temp from empty'); 
        
        % pick one of those modules
        if length(temp_from) > 1; temp_from = temp_from(randi(length(temp_from))); end
        % find a parcel from this module which has neighbours from other modules
        temp = contig_hemi{hf};
        temp(rand_out{hf} ~= temp_from,:) = 0;
        temp(:,rand_out{hf} == temp_from) = 0;
        [temp_parc,~] = find(temp); temp_parc = temp_parc(randi(length(temp_parc)));
        % relabel
        rand_out{hf}(temp_parc) = to_mod;
        
        % update mod_diff and mis
        for i = 1:1:max(module_in)
            mod_num_rand(i,2) = sum(rand_out{hf} == i);
            %mod_num_rand(i,1) = sum(rand_out{hi} == i);
        end
        
        % update mod_diff and mis
        mod_diff = mod_num_rand(:,2) - mod_num_in(:,2); % number of parcels to redistribute
        mis = find(mod_diff ~= 0); % modules where numbers don't match
        %disp(['relabel ' num2str(to_mod)])
        if isempty(mis); break; end % disp('break 2'); 
        continue
    end
    
    stuck = stuck+1;
    if stuck > 500; break; end %disp('unstuck (restarting)'); 
    %%% can get stuck here if doesn't enter below loop???
    
    for i = 1:1:length(neigh_mod{to_mod})
        if ~isempty(intersect(neigh_mod{neigh_mod{to_mod}(i)},from_vect))
            temp_fr = intersect(neigh_mod{neigh_mod{to_mod}(i)},from_vect);
            temp_fr = temp_fr(randi(length(temp_fr)));
            
            % find a parcel on boundary between modules temp_from and i
            temp = contig_hemi{hf};
            temp(rand_out{hf} ~= temp_fr,:) = 0;
            temp(:,rand_out{hf} ~= i) = 0;
            if sum(temp(:)) == 0;
                cnt=cnt+1;
                continue;
            end
            [ind_from,~] = find(temp);
            % pick a pair of nodes at random, relabel
            rel = randi(length(ind_from));
            rand_out{hf}(ind_from(rel)) = i;
            
            % find a parcel on boundary between modules i and to_mod
            % reduced contiguity matrix to identify neighbouring parcels from "from" and "to" modules
            temp = contig_hemi{hf};
            temp(rand_out{hf} ~= i,:) = 0;
            temp(:,rand_out{hf} ~= to_mod) = 0;
            if sum(temp(:)) == 0;
                cnt=cnt+1;
                continue;
            end
            [ind_from,~] = find(temp);
            % pick a pair of nodes at random, relabel
            rel = randi(length(ind_from));
            rand_out{hf}(ind_from(rel)) = to_mod;
            
            % update mod_diff and mis
            for i = 1:1:max(module_in)
                mod_num_rand(i,2) = sum(rand_out{hf} == i);
                %mod_num_rand(i,1) = sum(rand_out{hi} == i);
            end
            
            % update mod_diff and mis
            mod_diff = mod_num_rand(:,2) - mod_num_in(:,2); % number of parcels to redistribute
            mis = find(mod_diff ~= 0); % modules where numbers don't match
            if isempty(mis); break; end % disp('break 3'); 
            length(mis);
            
            continue
        end
    end
    
end % while ~isempty(mis)

if stuck > 500; continue; end % disp('restarting'); 
if exist('temp_from'); if isempty(temp_from); continue; end; end % disp('temp from empty - restarting'); 
if sum(temp(:))==0; continue; end % disp('temp empty - restarting'); 

finished = true;

end % while ~finished

% check that number of nodes in each module and each hemisphere matches
if ~isequal(mod_num_in,mod_num_rand); error('module size not preserved'); end

% combine output of two hemispheres
if size(module,1) == 1
    rand_mod = [rand_out{1} rand_out{2}]; % if input is a row vector
elseif size(module,2) == 1
    rand_mod = [rand_out{1}; rand_out{2}]; % if input is a column vector
end

%end % elseif two hemispheres

%%% store OVERALL initial size of modules for checking at the end
unique_rand = unique(rand_mod);
for i = 1:1:length(unique_rand)
    mod_size_output(i) = sum(rand_mod == unique_rand(i));
end
if ~isequal(mod_size_input,mod_size_output); error('overall module size not preserved'); end
%%%

%disp('finished!')

% optional condition - exclude cases where initial split module is contiguous

end
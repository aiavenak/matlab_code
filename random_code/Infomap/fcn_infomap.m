function [ci,L,R] = fcn_infomap(aij,iter)
% clear all
% close all
% 
% cd /Users/richardbetzel/Documents/MATLAB/_projects/_multiscale/_sub/sub01_108_rh/
% load sub01_108_rh
% 
% iter = 100;
% aij = fcn_flow_graph(aij,ones(108,1),time(40));

infodir = '/home/aiavenak/matlab/code/Infomap/infomap_undir';
fcn_mat2paj(aij,'dummy');

command = [infodir './infomap %i dummy.net 1'];

n = length(aij);
ci = zeros(n,iter);
L = zeros(iter,1);
R = zeros(iter,1);
for j = 1:iter
    
    r = randi(10000000);
    [~,~] = system(sprintf(command,r));
    R(j) = r;
    
    fid = fopen('dummy.clu','rt');
    com = textscan(fid,'%f','Headerlines',1);
    fclose(fid);
    com = com{1};
    
    fid = fopen('dummy.tree','rt');
    x = fgetl(fid);
    fclose(fid);
    
    for i = 1:4
        [str,x] = strtok(x,' ');
    end
    
    codeLength = str2double(str);
    
    ci(:,j) = com;
    L(j) = codeLength;
    
    delete('dummy_map.net');
    delete('dummy_map.vec');
    delete('dummy.clu');
    delete('dummy.map');
    delete('dummy.tree');
    
end

delete('dummy.net');
system('rm -rf ~/.Trash/*');
function [indsort,X,Y] = plot_communities_grid(ciu,aij,line_col,lims)

% input: 
% ciu - vector indicating RSN membership of nodes
% aij - connectivity matrix
% line_col - optional - color of the lines separating RSN
% limits - optional - max and min value of the matrix, eg [min(A(:)),max(A(:))];

if (nargin == 2) || ((nargin>2) && isempty(line_col))
    line_col = [1,1,1]
end

nc = max(ciu);
[ciu,indsort] = sort(ciu);
N = length(ciu)+0.5;

X = [];
Y = [];
for i = 1:nc
    ind = find(ciu == i);
    if ~isempty(ind)
        mn = min(ind) - 0.5;
        %mx = max(ind) + 0.5;
        %x = [mn mn mx mx mn NaN];
        x = [0 N NaN mn mn NaN]
        %y = [mn mx mx mn mn NaN];
        y = [mn mn NaN 0 N NaN];
        X = [X, x];
        Y = [Y, y];
    end
end

if nargin == 4
    imagesc(aij(indsort,indsort),lims); colorbar
else
   imagesc(aij(indsort,indsort)); colorbar 
end
hold on;
plot(X,Y,'color',line_col,'linewidth',1);
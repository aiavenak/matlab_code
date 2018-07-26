function [ordered_RGB, ordered_ind] = f_colors4scatterplot(node_attributes,num_clr_bins,clr_map,limits)

%   Input: node_atributes -- a list of numerical values; one value per node
%          num_clr_bins -- number of bins for the color-map
%          clr_map -- optional, the name of a matlab colormap or a list of
%          RGB colors (for example, from colormap.org); may or may not be
%          normalized. 

if nargin >= 3
    if isnumeric(clr_map)
        if ~(max(clr_map(:)) < 1)
            clr_map = clr_map./255;
        end
        if size(clr_map,1) ~= num_clr_bins
            ind_new = floor(linspace(1,size(clr_map,1),num_clr_bins));
            clr_map = clr_map(ind_new,:);
        end
    else            
        eval(sprintf('clr_map = %s(%d)',clr_map,num_clr_bins)); 
    end
else    
    clr_map = jet(num_clr_bins);    
end

[r,c] = size(node_attributes);
if (c > 1) && (r == 1)
    node_attributes = node_attributes';
elseif (c > 1) && (r > 1)
    warning('node attributes vector must be a single column vector \n\n');
end


if nargin == 4
    par = linspace(limits(1),limits(2),num_clr_bins);
else
    par = linspace(min(node_attributes),max(node_attributes),num_clr_bins);
end

[~,ordered_ind] = min(abs(bsxfun(@minus, par, node_attributes)),[],2);


ordered_RGB = clr_map(ordered_ind,:);





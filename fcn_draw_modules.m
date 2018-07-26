function [] = fcn_draw_modules(Aij,ci)

    %[ci, q] = modularity_louvain_und(Aij);
    [ci_ordered, indsort] = sort(ci); 
    Aij = Aij(indsort,indsort);
    imagesc(Aij); 
    hold;

    for i = 1:max(ci)
        [~, found] = find(ci_ordered == i);
        % find the first and last index of the module to outline
        x1 = min(found);
        x2 = max(found);
        % specify x and y coords of the corners of the square
        x = [x1-.5, x2+.5, x2+.5, x1-.5, x1-.5]; 
        % no need to find y coords because they are squares
        y = [x1-.5, x1-.5, x2+.5, x2+.5, x1-.5];
        % plot the square 
        plot(x,y,'w-', 'LineWidth', 2);
        shg;
    end
end

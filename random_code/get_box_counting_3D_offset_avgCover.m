function fractal = get_box_counting_3D_offset_avgCover(object,r,offset_param)

if isrow(r)
    r=r';
end
if any(mod(r,1))
    error('box sizes must be integers!')
end

nnz_object = nnz(object);
[width_i,width_j,width_k] = size(object);
        
D0num = zeros(numel(r),1);
D1num = zeros(numel(r),1);
D2num = zeros(numel(r),1);
Lnum = zeros(numel(r),1);

fractal = cell(length(offset_param),1);

%Lnum2 = zeros(numel(r),1);

num_offset = max(offset_param);

ind_offset = 1;

nnzNcum = zeros(numel(r),1);
entropycum = zeros(numel(r),1);
corrcum = zeros(numel(r),1);
laccum = zeros(numel(r),1);


for offset=1:num_offset % for each offset
    for index_size = 1:numel(r)
        cell_size = r(index_size);

        if (cell_size == 1) %max resolution: box is pixel-size
            D0num(index_size) = nnz_object;
            D1num(index_size) = -log(1/nnz_object); %simplified equation
            D2num(index_size) = 1/nnz_object; %num_points_object*((1/num_points_object)^2);
            Lnum(index_size) = std(object(:))^2/mean(object(:))^2; %heterogeneidad a la que el objeto esta esparcido a distintas escalas

            % N=object;

        else

            mod_i = mod(width_i,cell_size);
            mod_j = mod(width_j,cell_size);
            mod_k = mod(width_k,cell_size);

            %if (mod_i>0 || mod_j>0 || mod_k>0) || (offset==1) %if not perfect grid or is 1st run

                %% starting points for object in grid (monte-carlo offset strategy)
                offset_i = ceil(rand*cell_size);
                offset_j = ceil(rand*cell_size);
                offset_k = ceil(rand*cell_size);

                %% dimensions for grid covering object
                if offset_i > cell_size-mod_i
                    width_i_off = width_i + cell_size-mod_i + cell_size;
                else
                    width_i_off = width_i + cell_size-mod_i; %cell_size-mod_i + cell_size;
                end
                if offset_j > cell_size-mod_j
                    width_j_off = width_j + cell_size-mod_j + cell_size;
                else
                    width_j_off = width_j + cell_size-mod_j; %cell_size-mod_i + cell_size;
                end
                if offset_k > cell_size-mod_k
                    width_k_off = width_k + cell_size-mod_k + cell_size;
                else
                    width_k_off = width_k + cell_size-mod_k; %cell_size-mod_i + cell_size;
                end

                %% set object in grid at current corner
                object_offset = false(width_i_off,width_j_off,width_k_off);
                object_offset(...
                    offset_i : offset_i+width_i-1,...
                    offset_j : offset_j+width_j-1,...
                    offset_k : offset_k+width_k-1)...
                    = object;

                if nnz(object_offset)~=nnz_object
                    error('grid did not totally cover the object (%d out of %d at box-size %d)!',nnz(object_offset),nnz_object,cell_size)
                end
                Noff = box_count_3d_c(double(object_offset),cell_size); % c implementation of the box counting al

                if sum(Noff(:))~=nnz_object
                    error('grid did not totally cover the object (%d out of %d at box-size %d)!',sum(Noff(:)),nnz_object,cell_size)
                else
                    %% update results with new averaged grid solution
                    N = Noff; % N stores the minimum coverage obtained
                    nnzNcum(index_size) = nnzNcum(index_size) + nnz(N);
                    P = N./nnz_object;
                    P2=P;
                    P2(P2==1)=0;                    
                    entropycum(index_size) = entropycum(index_size) + sum(-P2(P2~=0).*log(P2(P2~=0))); %hay que ajustar dimensiones
                    corrcum(index_size) = corrcum(index_size) + sum(P(:).^2);
                    laccum(index_size) = laccum(index_size) + std(N(:))^2/(mean(N(:))^2);
                    D0num(index_size) = round(nnzNcum(index_size)/offset); %avg 
                    D1num(index_size) = entropycum(index_size)/offset; % avg 
                    D2num(index_size) = corrcum(index_size)/offset; % avg
                    Lnum(index_size) = laccum(index_size)/offset; % avg
                end
            %end
        end
    end
    %% if it is time to save fractal results
    if any(offset==offset_param)
        fractal{ind_offset} = struct('boxSizes',r,'D0num',D0num,'D1num',D1num,'D2num',D2num,'Lnum',Lnum,'V',nnz_object,...
            'offset',offset_param,'dimensions',[width_i,width_j,width_k]);
        ind_offset = ind_offset + 1;
    end

end

     
%figure, loglog(1./r,D0num,'s-'); 


% values = polyfit(log(1./r),log(D0num),1);
% D0 = values(1);
% RD0 = corr(log(1./r),log(D0num))^2; %R-square of the D0 fractal dimension
% values = polyfit(log(1./r),D1num,1);
% D1 = values(1);
% values = polyfit(log(1./r),log(D2num),1);
% D2 = -values(1);


%values = polyfit(log(1./r),log(Lnum),1);
%L = values(1);

%xdata = r;
%ydata = Lnum;
%L = lsqcurvefit(@(x,xdata) (x(2)./(xdata.^x(1))) + x(3), [1 1 1], xdata, ydata); %alpha, beta, gamma
%L = exp(values(2));




function plot_net(xyz,node_sz,color_indx,cij,flag_edges)

%FCN_PLOT_NET       plot of network
%
%   [F,H] = FCN_PLOT_NET(XYZ,NODE_SZ,COLOR_INDX,CIJ) generates plot of
%   network with connectivity matrix CIJ.
%
%   Inputs:             XYZ     xyz node coordinates
%                   NODE_SZ     vector of node sizes
%                COLOR_INDX     vector used to specify node color
%                       CIJ     weighted/binary symmetric connectivity
%                               matrix
%
%   Outputs:              F     figure handle
%                         H     axis handle
%
%   NOTE: for visualization, it maybe wise to threshold CIJ beforehand.
%
%   Richad Betzel, Indiana University, 2012

%nodes of interest NOI

p = size(xyz,2);
in1 = node_sz ~= 0;
in2 = color_indx ~= 0;
ind = in1 & in2;

%cij = cij(ind,ind);
cij  = (cij./max(cij(:))).*50;
noi = size(cij,1);

node_sz = node_sz(ind)+7;
color_indx = color_indx(ind);
xyz = xyz(ind,:);

clrs = jet(64);
color_indx = (color_indx - min(color_indx))./max(color_indx);
color_indx = floor(color_indx*63) + 1;

f = figure; 
h = axes;
hold(h,'on');

if flag_edges
    for i = 1:noi-1
        for j = i+1:noi
            if cij(i,j)
                xx = [xyz(i,1),xyz(j,1)];
                yy = [xyz(i,2),xyz(j,2)];
                zz = [xyz(i,3),xyz(j,3)];
                plot3(xx,yy,zz,'k','linewidth',.5);
            end
        end
    end
end

for i = 1:length(ind)
    switch p
        case 2
            plot(xyz(i,1),xyz(i,2),'ko','markersize',node_sz(i)...
                ,'markerfacecolor',clrs(color_indx(i),:));
        case 3
            plot3(xyz(i,1),xyz(i,2),xyz(i,3),'ko'...
                ,'markersize',node_sz(i)...
                ,'markerfacecolor',clrs(color_indx(i),:));
    end
end

axis image;
axis off;
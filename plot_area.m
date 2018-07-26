function hplot=plot_area(m1,m2,hfig,color)

figure(hfig);

if size(m1,1)>1
    m1_5=prctile(m1,2.5);
    m1_50=prctile(m1,50);
    m1_95=prctile(m1,97.5);
else
    m1_5=m1;
    m1_50=m1;
    m1_95=m1;
end
m2_5=prctile(m2,2.5);
m2_50=prctile(m2,50);
m2_95=prctile(m2,97.5);


X=[m1_5,m1_95(end:-1:1)];
Y=[m2_5,m2_95(end:-1:1)];

hfill = fill(X,Y,color);
set(hfill,'LineStyle','none');
alpha(hfill,0.25);
hold on
hplot = plot(m1_50,m2_50,'-','Color',color,'LineWidth',1);


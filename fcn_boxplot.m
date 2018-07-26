function ph = fcn_boxplot(Y,x,width)
% clear all
% close all
% clc
% 
% load Y
% x = 10;

miny = min(Y);
maxy = max(Y);
prcty = prctile(Y,[25,50,75]);

% width = 1;
boxx = [0, width, width, 0, 0] + x - (width/2);
boxy = [prcty(1), prcty(1), prcty(3), prcty(3), prcty(1)];

medx = [0, width] + x - (width/2);
medy = prcty(2)*ones(1,2);

xtop = (width/2)*ones(1,2) + x - (width/2);
ytop = [prcty(3), maxy];

xbot = xtop;
ybot = [prcty(1), miny];

% xvals = [boxx, nan, medx, nan, xtop, nan, xbot] + x - (width/2);
% yvals = [boxy, nan, medy, nan, ytop, nan, ybot];

ph = plot(boxx,boxy,xtop,ytop,xbot,ybot,medx,medy);
% ph = plot(xvals,yvals);
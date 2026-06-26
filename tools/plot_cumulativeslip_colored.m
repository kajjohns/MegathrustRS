function plot_cumulativeslip_colored(ss,t,V,slip)
% script to 2-d plot cumulative slip coloured as inter-,co- or post-seismic process along a transect/depth of the fault.
% INPUTS
% ss - fault structure or object (unicycle format preferred)
% ss.N - number of fault segments
% ss.W - vector containing width of each fault segment in the transect
% V,slip - velocity and slip as a matrix
% t - time vector
% adapted from James Moore's eqcycle2d repository (on github)
% Rishav Mallick, EOS, 2020

cspec=[cmap('steelblue',150,10,10);flipud(cmap('orange',100,57,10));flipud(cmap('orangered',100,20,25))];


dipdist = zeros(ss.N,1);
for i = 1:ss.N
    if i==1
        dipdist(i) = ss.W(1)/2;
    else
        dipdist(i) = dipdist(i-1) + ss.W(i);
    end
end

set(gcf,'Color','White','name','Time Evolution')
down = 2;
downt = 10;
% down = 1;
% downt = 1;

ydata=repmat(dipdist(1:down:end)./1e3,[1,size(t(1:downt:end))])';
xdata=slip(1:downt:end,1:down:end);
zdata=(V(1:downt:end,1:down:end));

xsample=linspace(min(min(xdata(:))),max(max(xdata(:))),1e3);

test1=griddata(xdata(:),ydata(:),zdata(:),xsample,dipdist(1:down:end)./1e3);

pcolor(xsample,dipdist(1:down:end)./1e3,test1), shading flat
alpha 0.7
set(gca,'YDir','reverse','Fontsize',15,'ColorScale','log')
h=colorbar('Location','NorthOutside');
% caxis(10.^[-10 0.5])
caxis(10.^[-11 log10(max(V(:)))])
colormap(cspec);
title(h,'V (m/s)')
xlabel('Accumulated slip (m)')
ylabel('Down-dip distance (km)');


end
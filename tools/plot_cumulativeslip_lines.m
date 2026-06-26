function index = plot_cumulativeslip_lines(ss,t,V,slip,skip)
% script to 2-d plot cumulative slip coloured as inter-,co- or post-seismic process along a transect/depth of the fault.
% INPUTS
% ss - fault structure or object (unicycle format preferred)
% ss.N - number of fault segments
% ss.W - vector containing width of each fault segment in the transect
% V,slip - velocity and slip as a matrix
% t - time vector
% Rishav Mallick, EOS, 2019

dipdist = zeros(ss.N,1);
for i = 1:ss.N
    if i==1
        dipdist(i) = ss.W(1)/2;
    else
        dipdist(i) = dipdist(i-1) + ss.W(i);
    end
end

j = 1;
% initialize arrays for Yn = integrated array; 
% Vn = velocity extracted from array of derivatives
index =[];
Veq = 1e-3;

for i = 1:length(t)-1
    if i == 1
        index(j) = i;
        j = j+1;
    elseif i > 1 && max(V(i,:)) >=Veq  
        if ((t(i) - t(index(j-1))) > 3) % number of seconds for EQ
            index(j) = i;
            j = j+1;
        end
    else
        if ((t(i) - t(index(j-1))) > 1/365*3.15e7) % number of days for slow process
            index(j) = i;
            j = j+1;
        end
    end
end
Vn = V(index',:);
slipn = slip(index',:);


count = 1;
for i = 1:length(Vn(:,1))
    if max(Vn(i,:)) >= Veq 
        plot(slipn(i,:),dipdist./1e3,'-','LineWidth',0.01,'Color',rgb('orangered')), hold on
    elseif (max(Vn(i,:))<Veq && max(Vn(i,:))>1e1*mean(ss.Vpl))
        if count>=skip/50
            plot(slipn(i,:),dipdist./1e3,'-','LineWidth',0.1,'Color',rgb('forestgreen')), hold on
            count=1;
        else
            count=count+1;
        end
    else
        if count>=skip
            plot(slipn(i,:),dipdist./1e3,'-','LineWidth',1,'Color',rgb('steelblue')), hold on
            count = 1;
        else
            count = count+1;
        end
    end    
end
ylabel('Distance along fault (km)')
xlabel('Cumulative Slip (m)')
axis tight, grid off
set(gca,'FontSize',15,'Color','none','YDir','reverse')


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ode output file name
savefn = 'ode_elastic_simple_sequence';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read the output file
% data = Read(filename,[cidxs],[ridxs],[isStride])
%   Read in the contents of a SSD file. If the optional argument cidxs is
%   provided, then data contains the columns listed in cidxs in ascending
%   order. cidxs can contain column numbers that exceed the number of rows in
%   the file. [ridxs] does the equivalent for rows. Set cidxs=[] if you want
%   to specifiy ridxs but not cidxs. If isStride=1 is given as the third
%   argument, cidxs is interpreted as a stride and should be a scalar.

stride=5;  %only read in every 'stride' time steps
x = SaveStreamData('Read',savefn,stride,[],1);
t = x(1,:)';  %times
x = x(2:end,:)';
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% unwrap x for velocity, slip, and stress, etc.
nel= size(FaultPatches,1);
vel = x(:,1:nel);
theta =  x(:,nel*1+1:nel*2);
tau = x(:,nel*2+1:nel*3);
phi =  x(:,nel*3+1:nel*4);
p =  x(:,nel*4+1:nel*5);
B = x(:,nel*5+1:nel*6);
slip = x(:,nel*6+1:nel*7);

dist = cumsum(lengths);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot b-value with time
figure
h=pcolor(dist/1000,t,B)
set(h, 'EdgeColor', 'none');
xlabel('distance along fault (km)')
ylabel('time, years')
colormap(jet)
colorbar
title('b')
set(gca,'fontsize',15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot log10 slip rates with time as color plots
figure
h=pcolor(dist/1000,t,log10(vel))
set(h, 'EdgeColor', 'none');
xlabel('distance along fault (km)')
ylabel('time, years')
colormap(jet)
colorbar
caxis([-4 1])
title('log_{10} slip rate (m/yr)')
set(gca,'fontsize',15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot log10 slip rates with solver time step as color plots
figure;hold on
imagesc(dist/1000,1:length(t),log10(vel))
colormap(jet)
colorbar
title('log_{10} slip rate (m/yr)')
xlabel('distance along fault (km)')
ylabel('solver time step')
caxis([-4 1])
set(gca,'fontsize',15)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot slip with time as color plots
figure
h=pcolor(dist/1000,t,slip)
set(h, 'EdgeColor', 'none');
xlabel('distance along fault (km)')
ylabel('time, years')
colormap(jet)
colorbar
title('cumulative slip (m)')
set(gca,'fontsize',15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot stress with solver time step as color plots
figure;hold on
imagesc(dist/1000,1:length(t),tau)
colormap(jet)
colorbar
title('stress(Pa)')
xlabel('distance along fault (km)')
ylabel('solver time step')
set(gca,'fontsize',15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot state with solver time step as color plots
figure;hold on
imagesc(dist/1000,1:length(t),theta)
colormap(jet)
colorbar
title('state (yr)')
xlabel('distance along fault (km)')
ylabel('solver time step')
set(gca,'fontsize',15)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot pore pressure with solver time step as color plots
figure;hold on
h=pcolor(dist/1000,t,p*1e-6)
set(h, 'EdgeColor', 'none');
colormap(jet)
colorbar
title('pore pressure (MPa)')
xlabel('distance along fault (km)')
ylabel('time')
set(gca,'fontsize',15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot slip evolution
figure;
yr2sec = 60*60*24*365;
ss.N = nel; ss.W=lengths;ss.Vpl = srate/yr2sec; plot_cumulativeslip_colored(ss,t/yr2sec,vel/yr2sec,slip)
%hold on; index = plot_cumulativeslip_lines(ss,t,vel,slip,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot velocity and displacement at points
distances = [140];   %distance along fault (km)
for k=1:length(distances)

    pos = find(dist>distances(k)*1000);
    pos = pos(1);
    figure
    semilogy(t,vel(:,pos))
    
   
    grid on
    title(['velocity at ' num2str(distances(k)) ' km along fault'])
    xlabel('time (yr)')
    ylabel('slip rate (m/yr)')
    set(gca,'fontsize',15)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOTE:  
% byou need to run the setup_simulation.m file associated 
% with ode file that is loaded below in order to load geometry 
% and all other parameters


%name of existing ode output file for extracting initial conditions  
savefn_existing = 'ode_elastic_dilatancy1';

%name of output file for new simulation results 
savename = 'ode_elastic_dilatancy1_stresschange';

%time at which to impose stress change (years)
t_impose = 20.8; 

%vector of induced stress change (Pa, must be size (nel,1) where nel is number of elements) 
delta_tau = 10^5*exp(-(depths-20000).^2/5000^2);  %gaussian-shaped stress change

%ramp time (time over which stress change ramps linearly from zero to full amount)
yr2sec = 60*60*24*365;
duration = 60/yr2sec; %seconds

%stop time in years
stopt = 10;  

%save every n-th time step
gamb.ode.saveEvery = 10;




% END input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%plot stress change

figure
scatter(centers(:,1),-centers(:,2),10,delta_tau,'fill')
axis ij
axis equal
colorbar
colormap(jet)
title('stress change, Pa')
ylabel('depth (km)')
xlabel('distance from trench (km)')
set(gca,'fontsize',15)
axis equal
grid on
drawnow

%specify file for outputs
dr = '.';
gamb.ode.savefn = [dr '/' savename ];


%load times and find nearest output time to t_timpose
t = SaveStreamData('Read',savefn_existing,[],1);
ind = find(t>t_impose); ind = ind(1);

%load simulation outputs at requested time
x = SaveStreamData('Read',savefn_existing,ind,[]);
x0 = x(2:end);  %initial values for simulation, discard time



%elastic stressing rate
nel= size(FaultPatches,1);
[gamb.hmat.id nnz] = hm_mvp('init',gamb.hmat.shear_savefn,4);
fprintf(1,'Compression factor = %3.1f\n',(nel)^2/nnz);
v = srate*ones(nel,1);
stressing_rate = -hm_mvp('mvp',gamb.hmat.id,v);


[t,x] = get_ratestate_dilatancy_initial(x0,delta_tau,duration,srate,dc,stopt,a,b,sigma,mu_0,stressing_rate,c_star,eps,beta,pinf,a_b_slope,a_b_cut,logVcut,use_ikari);

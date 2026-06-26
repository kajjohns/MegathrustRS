%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%save every n-th time step
gamb.ode.saveEvery = 20;

%name of output file for simulation results 
savename = 'ode_elastic_large_Linf';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



gamb.ode.savefn = [gamb.ode.savefn savename];




%elastic stressing rate
nel= size(FaultPatches,1);
[gamb.hmat.id nnz] = hm_mvp('init',gamb.hmat.shear_savefn,4);
fprintf(1,'Compression factor = %3.1f\n',(nel)^2/nnz);
v = srate*ones(nel,1);
stressing_rate = -hm_mvp('mvp',gamb.hmat.id,v);


%compute slip evolution -- ODE solver

[t,x] = get_ratestate_dilatancy(srate,dc,stopt,a,b,sigma,mu_0,stressing_rate,c_star,eps,beta,pinf,a_b_slope,a_b_cut,logVcut,use_ikari);

vel = x(:,1:nel);
theta = x(:,nel+1:nel*2);
slip = x(:,nel*6+1:nel*7);
tau = x(:,nel*2+1:nel*3);

for loop=1:size(tau,2)
    tau_rate(:,loop) = gradient(tau(:,loop),t);
end

hm_mvp('cleanup',gamb.hmat.id);


figure
h=pcolor(depths,t,log10(vel))
set(h, 'EdgeColor', 'none');
xlabel('depth (m)')
ylabel('time, years')
colormap(jet)
colorbar
title('log_{10} slip rate (m/yr)')

figure
h=pcolor(depths,1:length(t),log10(vel))
set(h, 'EdgeColor', 'none');
xlabel('depth (m)')
ylabel('solver time step')
colormap(jet)
colorbar
title('log_{10} slip rate (m/yr)')

figure
time_slice = 12000;
plot(depths,vel(time_slice,:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%save every n-th time step

gamb.ode.saveEvery = 10;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%specify file for outputs
dr = '.';
gamb.ode.savefn = [dr '/' savename ];







%elastic stressing rate
nel= size(FaultPatches,1);
[gamb.hmat.id nnz] = hm_mvp('init',gamb.hmat.shear_savefn,4);
fprintf(1,'Compression factor = %3.1f\n',(nel)^2/nnz);
v = srate*ones(nel,1);
stressing_rate = -hm_mvp('mvp',gamb.hmat.id,v);


%compute slip evolution -- ODE solver

[t,x] = get_ratestate_dilatancy(srate,dc,stopt,a,b,sigma,mu_0,stressing_rate,c_star,eps,beta,pinf,a_b_slope,a_b_cut,logVcut,use_ikari);

vel = x(:,1:nel);


slip = x(:,nel*3+1:nel*4);
tau = x(:,nel*2+1:nel*3);

hm_mvp('cleanup',gamb.hmat.id);


figure
h=pcolor(depths,t,log10(vel))
set(h, 'EdgeColor', 'none');
xlabel('depth (m)')
ylabel('time, years')
colormap(jet)
colorbar
title('log_{10} slip rate (m/yr)')

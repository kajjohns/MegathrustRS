function [t,x] = get_ratestate_elastic(srate,d_c,stopt,A,B,sigma,mu_0,stressing_rate,coarse2fine,fine2coarse,visco_times,G_tau_rate_visco)
                

%number of patches
nps=length(d_c);


   
%pre-quake slip rate
v_ss = srate*ones(nps,1);  %m/yr
v_ss(B>A) = 10^-3;
v_star = v_ss;  %Note: v_star has no effect on the solution!!

%pre-earthquake stress
%tau0    =    sigma.*(mu_0 + (A-B).*log(v_ss./v_star));
tau0    =    (sigma-pinf).*(mu_0 + (A-B).*log(v_ss./v_star));

%pre-earthquake state variable
theta0 = d_c./v_ss;
sigma0=sigma;

phi0 = zeros(size(theta0));  %solution only depends on dphi_dt
p0 = pinf*ones(size(theta0));

% SEISMIC RADIATION DAMPING TERM
yr2sec = 60*60*24*365;
vs = 5*1000*yr2sec;  %5 km/s --> m/yr
eta = 3*10^10/(2*vs);  %eta = mu/2*vs  Pa-yr/m

eta = 100*eta;  %increase damping

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%RUNGE KUTTA INTEGRATION
% START AND STOP TIME
t0      =  0;
tf      =  stopt;     % time in years

   
x0      =   [v_ss; theta0; tau0; phi0; p0; u0];



opt = Ode('Options');
global gamb;

 OutputFcn('Init',gamb.ode.savefn,'saveEvery',gamb.ode.saveEvery);
 opt = Ode_visco('Options');
 opt.OutputFcn = @(varargin)OutputFcn('Fn',varargin{:});
 
 
 Ode_elastic('Pair23','ratestate_3D',[t0 tf],x0,opt,...
     srate,coarse2fine,fine2coarse,visco_times,G_tau_rate_visco,stressing_rate,d_c,A,B,sigma,eta,nps);
 OutputFcn('Cleanup');
 x = SaveStreamData('Read',gamb.ode.savefn);

 t = x(1,:)';
 x = x(2:end,:)';
 

 
%times=linspace(t0,tf,5000);
%[t,x] = ode23('ratestate_3D', times,x0,opt,d_c,A,B,sigma,eta,v_ss,nps,stressing_rate);

 


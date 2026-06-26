function [t,x] = get_ratestate_dilatancy_initial(x0,delta_tau,duration,sliprate,d_c,stopt,A,B,sigma,mu_0,stressing_rate,c_star,eps,beta,pinf,a_b_slope,a_b_cut,logVcut,use_ikari)


%number of patches
nps=length(d_c);


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

   
opt = Ode('Options');
global gamb;


 OutputFcn('Init',gamb.ode.savefn,'saveEvery',gamb.ode.saveEvery);
 opt = Ode('Options');
 opt.OutputFcn = @(varargin)OutputFcn('Fn',varargin{:});
 
 
 Ode('Pair23','ratestate_3D_dilatancy_ramp',[t0 tf],x0,opt,...
     stressing_rate,delta_tau,duration,d_c,A,sigma,eta,c_star,eps,beta,pinf,mu_0,nps,a_b_slope,logVcut,use_ikari);
 OutputFcn('Cleanup');

 if nargout==2
 x = SaveStreamData('Read',gamb.ode.savefn);
 t = x(1,:)';
 x = x(2:end,:)';
 end
 

 


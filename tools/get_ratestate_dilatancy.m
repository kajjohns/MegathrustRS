function [t,x] = get_ratestate_dilatancy(sliprate,d_c,stopt,A,B,sigma,mu_0,stressing_rate,c_star,eps,beta,pinf,a_b_slope,a_b_cut,logVcut,use_ikari)



%number of patches
nps=length(d_c);

%compute velocity-dependent b value


   
%pre-quake slip rate
v_ss = sliprate*ones(nps,1);  %m/yr
v_ss(B>A | use_ikari) = 0.1*v_ss(B>A | use_ikari);
v_star = v_ss;  


%ikari a-b
maxb = A - a_b_cut;
B(use_ikari) = A(use_ikari) - (a_b_slope*(log10(v_ss(use_ikari))-logVcut)+a_b_cut);
B(use_ikari) = min(B(use_ikari),maxb(use_ikari));



%pre-earthquake stress
%tau0    =    sigma.*(mu_0 + (A-B).*log(v_ss./v_star));
tau0    =    (sigma-pinf).*(mu_0 + (A-B).*log(v_ss./v_star));

%pre-earthquake state variable
theta0 = d_c./v_ss;

phi0 = zeros(size(theta0));  %solution only depends on dphi_dt
p0 = pinf;


%initiate displacements
u0 = zeros(size(v_ss));

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

   
x0      =   [v_ss; theta0; tau0; phi0; p0; B; u0];



opt = Ode('Options');
global gamb;


 OutputFcn('Init',gamb.ode.savefn,'saveEvery',gamb.ode.saveEvery);
 opt = Ode('Options');
 opt.OutputFcn = @(varargin)OutputFcn('Fn',varargin{:});
 
 
 Ode('Pair23','ratestate_3D_dilatancy',[t0 tf],x0,opt,...
     stressing_rate,d_c,A,sigma,eta,c_star,eps,beta,pinf,mu_0,nps,a_b_slope,logVcut,use_ikari);
 OutputFcn('Cleanup');
 x = SaveStreamData('Read',gamb.ode.savefn);

 t = x(1,:)';
 x = x(2:end,:)';
 

 


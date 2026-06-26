function [out, serr] = ratestate_3D_dilatancy(t,x,srate,deltau,duration,d_c,A,sigma,eta,c_star,eps,beta,pinf,mu,nps,a_b_slope,logVcut,use_ikari)

%   t = times
%   x = model vector

global gamb;

%  unwrap parameters
v     = x(1:nps);
theta = x(nps+1:2*nps);
tau   = x(2*nps+1:3*nps);
phi   =  x(3*nps+1:4*nps);
p     =  x(4*nps+1:5*nps);
B     =  x(5*nps+1:6*nps);
u     =  x(6*nps+1:7*nps);
    

%return rates

dtheta_dt = 1-theta.*v./d_c; %Aging Law
dphi_dt = -eps.*dtheta_dt./theta;  %equation 17, Segall and Rice, 1995
dp_dt = -dphi_dt./beta + c_star.*(pinf-p);

if t<duration
    ramp_rate = deltau/duration;
else
    ramp_rate = 0*deltau;
end

dtau_dt = hm_mvp('mvp',gamb.hmat.id,v) + srate + ramp_rate;

dB_dv = zeros(nps,1);
dB_dv(use_ikari) =  -a_b_slope/log(10)./v(use_ikari);
dB_dv(log10(v)<logVcut) = 0;


dv_dt =  (eta./(sigma-p) + A./v + log(theta).*dB_dv).^-1.*((dtau_dt+tau.*dp_dt./(sigma-p))./(sigma-p) - B.*dtheta_dt./theta);
du_dt = v;  

dB_dt = dB_dv.*dv_dt;



out(1:nps) = dv_dt;
out(1+nps:2*nps) = dtheta_dt;
out(1+2*nps:3*nps) = dtau_dt;
out(1+3*nps:4*nps) = dphi_dt;
out(1+4*nps:5*nps) = dp_dt;
out(1+5*nps:6*nps) = dB_dt;
out(1+6*nps:7*nps) = du_dt;



out = out';
serr=0;

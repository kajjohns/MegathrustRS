%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT section

clear all


%file name for elastic Greens Functions
elastic_name = 'elastic_GFs';
%load GFs
load(elastic_name)


%name of output file for simulation results 
savename = 'ode_elastic_simple_interface3';


%specify long-term slip rate (m/yr)
srate = 0.02;  

%nominal coefficient of friction
mu_0=0.2;


%depth-dependent friction parameters
%NOTE:  depth means depth below seafloor/ground, not sea level, because
%the seafloor is the top of the elasic medium

%specify a, b, d_c at points at depth (km)
%values will be assigned to patches by linear interpolation of the
%specified values
%param_depth -- depth (km) at which friction paramter values are assigned
%assignments at depth=0 and depth at bottom edge of the model fault required
%you may add more depths 
param_depth(1) = 0;  %km, top edge (required)
param_depth(2) = 4;  
param_depth(3) = 4.5; 
param_depth(4) = 10.5; 
param_depth(5) = 11; 
param_depth(6) = 15; 
param_depth(7) = 15.01;
param_depth(8) = 25;
param_depth(9) = 25.01;
param_depth(10) = -FaultPatches(end,4);  % bottom edge (required, set to bottom edge of interface) 




%values corresponding with intervals above

%switch between velocity-dendent a-b and fixed a,b values
use_ikari(1) = false;  %true for velocity-dependent a-b formulation, false for fixed values
use_ikari(2) = false;
use_ikari(3) = false;
use_ikari(4) = false;
use_ikari(5) = false;
use_ikari(6) = false;
use_ikari(7) = true;
use_ikari(8) = true;
use_ikari(9) = false;
use_ikari(10) = false;

%if use_ikari=true, specify parameters below
a_b_slope =   0.001333;  %slope of a-b vs. log10 slip rate
a_b_cut = -0.0005;  %minimum a-b
logVcut = -1.5;   %slip rate at minimum a-b value (m/yr)

%specif a, b within intervals below
a(1) = 0.01;  %unitless
a(2) = 0.01;
a(3) = 0.01;
a(4) = 0.01;
a(5) = 0.01;
a(6) = 0.01;
a(7) = 0.01;
a(8) = 0.01;
a(9) = 0.01;
a(10) = 0.01;

%note, b values below for fixed a-b. If use_ikari=true, then b values are a
%function of slip rate
b(1) = 0.003; %unitless
b(2) = 0.005;
b(3) = 0.012;
b(4) = 0.012;
b(5) = 0.012;
b(6) = 0.002;
b(7) = 0.0105;
b(8) = 0.0105;
b(9) = 0.003;
b(10) = 0.003;


%critical slip distance
d_c(1) = 0.01;  %m
d_c(2) = 0.01;
d_c(3) = 0.01;
d_c(4) = 0.01;
d_c(5) = 0.001;
d_c(6) = 0.001;
d_c(7) = 0.001;
d_c(8) = 0.001;
d_c(9) = 0.001;
d_c(10) = 0.01;





%nominal pore pressure (as fraction of lithostatic pressure)
pore_pressure_fraction(1) = 0.7;
pore_pressure_fraction(2) = 0.7;
pore_pressure_fraction(3) = 0.7;
pore_pressure_fraction(4) = 0.7;
pore_pressure_fraction(5) = 0.7;
pore_pressure_fraction(6) = 0.7;
pore_pressure_fraction(7) = 0.7;
pore_pressure_fraction(8) = 0.7;
pore_pressure_fraction(9) = 0.9;
pore_pressure_fraction(10) = 0.9;

%maximum effective normal stress (effective normal stress, litho-pore pressure, will not exceed this value)
max_eff = 5e7;


%fluid parameters

c_star(1) = 1/2; %diffusivity, 1/yrs
c_star(2) = 1/2;
c_star(3) = 1/2;
c_star(4) = 1/2;
c_star(5) = 1/2;
c_star(6) = 1/2;
c_star(7) = 1/2;
c_star(8) = 1/2;
c_star(9) = 1/2;
c_star(10) = 1/2;

beta(1) = 1e-3*10^-6;  % bulk compressibility, 1/Pa
beta(2) = 1e-3*10^-6; 
beta(3) = 1e-3*10^-6; 
beta(4) = 1e-3*10^-6; 
beta(5) = 1e-3*10^-6; 
beta(6) = 1e-3*10^-6; 
beta(7) = 1e-3*10^-6; 
beta(8) = 1e-3*10^-6; 
beta(9) = 1e-3*10^-6; 
beta(10) = 1e-3*10^-6; 

eps(1) = 0;  %dilatancy coefficient,  1.7e-4 (Segall and Rice, quartz), Samualsen (4.7e-5 to 3e-4)
eps(2) = 0;
eps(3) = 0;   %Gap
eps(4) = 0*.1*10^-4;  %ETS
eps(5) = 0;
eps(6) = 0;
eps(7) = 0;
eps(8) = 0;
eps(9) = 0;
eps(10) = 0;

%stop time in years
stopt = 6000;  

% END input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ispc
    addpath hmmvp0.16_win
else
    addpath hmmvp0.16
end

addpath tools




%compute stressing rate due to detachment slip
pm = make_pm(FaultPatches);

%confining stress (lithostatic minus hydrostatic)
depths = -.5*(FaultPatches(:,2) + FaultPatches(:,4))*1000;  %convert to m
centers = (FaultPatches(:,1:2) + FaultPatches(:,3:4))/2;
lengths = 1000*sqrt( (FaultPatches(:,1)-FaultPatches(:,3)).^2 + (FaultPatches(:,2)-FaultPatches(:,4)).^2 );



%fix paramter designation on bottom edge
if isinf(param_depth(end))
    param_depth(end) = -FaultPatches(end,4);
end


%friction parameters -- interpolate 
a = interp1(param_depth*1000,a,depths);
b = interp1(param_depth*1000,b,depths);
dc = interp1(param_depth*1000,d_c,depths);
pore_pressure_fraction = interp1(param_depth*1000,pore_pressure_fraction,depths);
beta = interp1(param_depth*1000,beta,depths);
c_star = interp1(param_depth*1000,c_star,depths);
eps = interp1(param_depth*1000,eps,depths);

use_ikari_i = interp1(param_depth*1000,single(use_ikari),depths);
use_ikari = use_ikari_i>0.5;  %make logical


%lithostatic normal stress
sigma = 2700*9.81*depths;
pinf = pore_pressure_fraction.*sigma;  %nominal pore pressure outside of deforming fault zone

%implement maximum effective normal stress
ind = (sigma-pinf)>max_eff;
pinf(ind) = sigma(ind) - max_eff;

maxb = a - a_b_cut;


figure;
plot(a,depths/1000,'b','linewidth',3)
hold on
plot(b,depths/1000,'r','linewidth',3)
plot(maxb(use_ikari),depths(use_ikari)/1000,'g','linewidth',3)
legend('a','b','ikari maxb')
axis ij
set(gca,'fontsize',15)
ylabel('depth, km')

figure
scatter(centers(:,1),-centers(:,2),10,a-b,'fill')
axis ij
axis equal
colorbar
title('a-b')
ylabel('depth (km)')
xlabel('distance from trench (km)')
set(gca,'fontsize',15)


% h*, from Rice 1993, h*=2*mu*Dc/[pi*(B-A)], A=sigma*a, B=sigma*B
h_star = 2*3*10^10*dc./(pi*(b-a).*(sigma-pinf));
ind = b<a;
h_star(ind) = 3*10^10*dc(ind)./(b(ind).*(sigma(ind)-pinf(ind)));

%for velocity-dependent a-b, compute at maxb
% h*, from Rice 1993, h*=2*mu*Dc/[pi*(B-A)], A=sigma*a, B=sigma*B
h_star_ikari = 2*3*10^10*dc./(pi*(maxb-a).*(sigma-pinf));
ind = maxb<a;
h_star_ikari(ind) = 3*10^10*dc(ind)./(maxb(ind).*(sigma(ind)-pinf(ind)));


figure
plot(sigma*1e-6,depths)
hold on
plot(pinf*1e-6,depths)
plot((sigma-pinf)*1e-6,depths)
axis ij
ylabel('depth')
legend('lithostatic','pore pressure', 'effective stress')
xlabel('pressure, MPa')
set(gca,'fontsize',15)




figure
plot(log10(h_star),depths)
hold on
plot(log10(h_star_ikari(use_ikari)),depths(use_ikari))
plot(log10(lengths),depths)
title('log_{10} h-star and patch lengths (meters)')
ylabel('depth')
legend('h^*','h^* (ikari)','patch length')
axis ij
set(gca,'fontsize',15)




figure
logv=linspace(-5,5);
a_ikari = mean(a(use_ikari));
b_ikari = a_ikari - (a_b_slope*(log10(10.^logv)-logVcut)+a_b_cut);
b_ikari = min(b_ikari,a_ikari-a_b_cut);
plot(log10(10.^logv),a_ikari*ones(size(b_ikari)),'linewidth',3)
hold on
plot(log10(10.^logv),b_ikari,':','linewidth',3)
plot(log10(10.^logv),a_ikari-b_ikari,'linewidth',3)
set(gca,'fontsize',15)
grid on
legend('a (mean of assigned)','b','a-b')
ylabel('parameter value (unitless)')
xlabel('log10 slip rate (m/yr)')
title('velocity dependent a-b, only for use\_ikari=true')




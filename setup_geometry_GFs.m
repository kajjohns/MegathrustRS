
clear all

addpath tools
addpath viscoGFs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT PARAMETERS


%geometry of subduction interface
faultdip_trench=10;  %dip of interface at trench
faultdip_bottom=20;  %dip of interface at bottom 
x_trench = 0;   %x-position of trench
x_bottom = 238;   %x-position of bottom of interface (bottom of cold wedge)
z_bot = 45;    %depth of bottom of interface

%file name for saving elastic green's functions
name_elastic = 'elastic_GFs';


%patch length for elastic caclulations (km)
%maximum size -- actual patch size is smaller
pL = .2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




addpath tools


%build geometry -- note: 

[topx_interface,topz_interface,botx_interface,botz_interface] = make_geometry_interface(faultdip_trench,x_trench,x_bottom,faultdip_bottom,z_bot,pL);


figure
hold on
plot(topx_interface,-topz_interface,'r.')
axis ij
grid on
axis equal
drawnow
title('interface geometry (top edge of patches)')
ylabel('depth (km)')
xlabel('distance from trench (km)')
set(gca,'fontsize',15)




FaultPatches = [topx_interface' topz_interface' botx_interface' botz_interface'];

    

% ---- Software configuration
% Make hmmvp available
addpath hmmvp0.16;




global gamb;
gamb.use = 1;
% Optional H-matrix compression of the BEM matrix:
% Use H-matrix compression as implemented in hmmvp?
gamb.hmat.use = 1;
% Element-wise relative error of approximation. Should be ~ RelTol in the
% ODE integration.
gamb.hmat.rerr = 1e-3;
% Directory to which to save H-mat for future use:
dr = '.';



%generate shear stress matrix

% ODE software parameters:
% File to which to save ODE solution
ne = size(FaultPatches,1);
fn_dec = sprintf('ne%dhmat%drerr%3.2f',...
       ne,gamb.hmat.use,log10(gamb.hmat.rerr));
gamb.ode.savefn = [dr '/ode_' fn_dec '_' name_elastic ];

% Save at every k'th ODE time step
gamb.ode.saveEvery = 10; 

% H-matrix file name. 
gamb.hmat.shear_savefn = sprintf('%s/shear_Hmat_%s.dat',dr,fn_dec);
gamb.hmat.normal_savefn = sprintf('%s/normal_Hmat_%s.dat',dr,fn_dec);


%make stress-slip matrix, stressG, such that stress = stressG*slip;
CalcShear_edgestress('Make',FaultPatches,gamb.hmat.rerr,gamb.hmat.shear_savefn);
%CalcNormal('Make',FaultPatches,gamb.hmat.rerr,gamb.hmat.normal_savefn);


%prepare geometry parameters for saving
Geometry.faultdip_trench = faultdip_trench;
Geometry.x_trench = x_trench;
Geometry.x_bottom = x_bottom;
Geometry.faultdip_bottom = faultdip_bottom;
Geometry.z_bot = z_bot;
Geometry.pL = pL;


%save outptuts to mat file
save(name_elastic,'gamb','FaultPatches','Geometry')



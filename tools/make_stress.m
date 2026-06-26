function [Sxx,Sxz,Szz]=make_stress(pm,pm_small,mu,nu)

%convert km to m
pm(:,1:3)=1000*pm(:,1:3);
pm(:,6:7)=1000*pm(:,6:7);
%dip = pm(:,4);

pm_small(:,1:3)=1000*pm_small(:,1:3);
pm_small(:,6:7)=1000*pm_small(:,6:7);


npatches=size(pm_small,1);  %number of patches

%calculate patch centers
angle = 180-pm(:,5);
centerx = pm(:,6)+.5*pm(:,2).*cos(pm(:,4)*pi/180).*cos(angle*pi/180);
centery = pm(:,7)+.5*pm(:,2).*cos(pm(:,4)*pi/180).*sin(angle*pi/180);
centerz = -pm(:,3)+.5*pm(:,2).*sin(pm(:,4)*pi/180);

Xcenter=[centerx';centery';centerz'];


%calculate patch normals
%strike=90-pm(:,5);
%strikevec=[cos(strike*pi/180) sin(strike*pi/180) zeros(size(strike))];
%dipvec=[pm(:,6) pm(:,7) -pm(:,3)]-[centerx centery centerz];
%vecnorm=sqrt(dipvec(:,1).^2+dipvec(:,2).^2+dipvec(:,3).^2);
%vecnorm=repmat(vecnorm,1,3);

%need to be very careful about signs
%if dip<90, dipvec should point up
%if dip>90, dipvec should point down
%index = dip<90;
%dipvec=dipvec./vecnorm;  %change sign so that dipvec points up
%dipvec(index,:) = -dipvec(index,:);

%normal=cross(strikevec,dipvec,2);

Sxx=zeros(size(pm,1),size(pm_small,1));
Sxz=Sxx;
Szz=Sxx;

%calculate coseismic stress changes
for k=1:npatches
  % mss=[pm(k,:) 1 0 0]';
   %mds=[pm_small(k,:) 0 1 0]';
   
   %[U,D,Sss,flag]=disloc3d(mss,Xcenter,mu,.25);
   %[U,D,Sds,flag]=disloc3d(mds,Xcenter,mu,nu);
   

   m = [pm_small(k,3) pm_small(k,6) pm_small(k,2) pm_small(k,4)-180 1000];  %1 km of slip 
   [Sxx(:,k),Sxz(:,k),Szz(:,k)]=EdgeStress(m,Xcenter(1,:),Xcenter(3,:),nu,mu);
  
   %Sxx(:,k)=Sds(1,:)';
   %Sxz(:,k)=Sds(3,:)';
   %Szz(:,k)=Sds(6,:)';
   
   
      
end



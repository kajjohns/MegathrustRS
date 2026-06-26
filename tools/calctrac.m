function [tau,sig]=calctrac(pm,Sxx,Sxz,Szz)


dip = pm(:,4);



%calculate patch centers
angle = 180-pm(:,5);
centerx = pm(:,6)+.5*pm(:,2).*cos(pm(:,4)*pi/180).*cos(angle*pi/180);
centery = pm(:,7)+.5*pm(:,2).*cos(pm(:,4)*pi/180).*sin(angle*pi/180);
centerz = -pm(:,3)+.5*pm(:,2).*sin(pm(:,4)*pi/180);



%calculate patch normals
strike=90-pm(:,5);
strikevec=[cos(strike*pi/180) sin(strike*pi/180) zeros(size(strike))];
dipvec=[pm(:,6) pm(:,7) -pm(:,3)]-[centerx centery centerz];
vecnorm=sqrt(dipvec(:,1).^2+dipvec(:,2).^2+dipvec(:,3).^2);
vecnorm=repmat(vecnorm,1,3);

%need to be very careful about signs
%if dip<90, dipvec should point up
%if dip>90, dipvec should point down
index = dip<90;
dipvec=dipvec./vecnorm;  %change sign so that dipvec points up
dipvec(index,:) = -dipvec(index,:);

normal=cross(strikevec,dipvec,2);



%traction vector
T1ds = Sxx.*normal(:,1) + 0*normal(:,2) + Sxz.*normal(:,3);
%T2ds = 0;
T3ds = Sxz.*normal(:,1) + 0*normal(:,2) + Szz.*normal(:,3);


tau = T1ds.*dipvec(:,1) + T3ds.*dipvec(:,3);
sig = T1ds.*normal(:,1) + T3ds.*normal(:,3);


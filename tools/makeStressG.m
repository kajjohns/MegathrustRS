function G = makeStressG(patches)
        
pm = make_pm(patches);



%convert km to m
pm(:,1:3)=1000*pm(:,1:3);
pm(:,6:7)=1000*pm(:,6:7);
dip = pm(:,4);

%npatches=size(pm_small,1);  %number of patches

%calculate patch centers
angle = 180-pm(:,5);
centerx = pm(:,6)+.5*pm(:,2).*cos(pm(:,4)*pi/180).*cos(angle*pi/180);
centery = pm(:,7)+.5*pm(:,2).*cos(pm(:,4)*pi/180).*sin(angle*pi/180);
centerz = -pm(:,3)+.5*pm(:,2).*sin(pm(:,4)*pi/180);



centers=[centerx';centery';centerz'];


%calculate patch normals
strike = 90-pm(:,5);
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

   


 G = zeros(size(pm,1),size(pm,1));

 mu = 3e10;
 
 for k =1:size(pm,1)
    
     m=[pm(k,:) 0 1 0]';
     xloc = centers;
     
     
     [U, D, S] = disloc3d(m, xloc, mu, .25);

     
    
%     %calculate tractions
    T1=S(1,:)'.*normal(:,1)+S(2,:)'.*normal(:,2)+S(3,:)'.*normal(:,3);
    T2=S(2,:)'.*normal(:,1)+S(4,:)'.*normal(:,2)+S(5,:)'.*normal(:,3);
    T3=S(3,:)'.*normal(:,1)+S(5,:)'.*normal(:,2)+S(6,:)'.*normal(:,3);

    
    %component along dipvec
     S_ds = T1.*dipvec(:,1)+T2.*dipvec(:,2)+T3.*dipvec(:,3);

     
     
     %S_ten = T1.*normal(:,1)+T2.*normal(:,2)+T3.*normal(:,3);

     G(:,k) = S_ds;
     %Gten(:,i) = S_ten;
       
   
     
  %  end
  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
  end
  
  
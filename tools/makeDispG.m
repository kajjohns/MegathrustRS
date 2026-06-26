function [Geast,Gnorth,Gup] = makeDispG(patches,xpos)
        
pm = make_pm(patches);

%npatches=size(pm_small,1);  %number of patches

Geast = zeros(length(xpos),size(pm,1));
Gnorth = zeros(length(xpos),size(pm,1));
Gup = zeros(length(xpos),size(pm,1));

 
 for k =1:size(pm,1)
    
     m=[pm(k,:) 0 1 0]';
     xloc = [xpos';0*xpos';0*xpos'];
     
     
     [U, D, S] = disloc3d(m, xloc, 1, .25);

     
    Geast(:,k) = U(1,:)';
    Gnorth(:,k) = U(2,:)';
    Gup(:,k) = U(3,:)';
    
     
    
  end
  
  
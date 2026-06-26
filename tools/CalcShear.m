function varargout = CalcShear(fn,varargin)
% Rearranged code from make_stressG_ss.m.
  [varargout{1:nargout}] = feval(fn,varargin{:});

function Make(patches,rerr,savefn)
  g = CalcShear('Init',patches);
  CalcShear('Compress',g,rerr,savefn);
  
function g = Init(patches)
    
    

        
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

   

        g.pm = pm;
        g.centers=centers;
        g.strikevec=strikevec;
        g.normal=normal;

        g.mu = 3e10;  %everything is in m
        g.dipvec = dipvec;
         
     
  
function G = Fn(g,rs,cs)



  G = zeros(length(rs),length(cs));
  for(i = 1:length(cs))
    k = cs(i);
    
     m=[g.pm(k,:) 0 1 0]';
     xloc = g.centers;
     
     
     
     [U, D, S] = disloc3d(m, xloc(:,rs), g.mu, .49);

     
     normal = g.normal(rs,:);
%     %calculate tractions
    T1=S(1,:)'.*normal(:,1)+S(2,:)'.*normal(:,2)+S(3,:)'.*normal(:,3);
    T2=S(2,:)'.*normal(:,1)+S(4,:)'.*normal(:,2)+S(5,:)'.*normal(:,3);
    T3=S(3,:)'.*normal(:,1)+S(5,:)'.*normal(:,2)+S(6,:)'.*normal(:,3);

    
    %component along dipvec
     dipvec = g.dipvec(rs,:);
     S_ds = T1.*dipvec(:,1)+T2.*dipvec(:,2)+T3.*dipvec(:,3);

     
     
     %S_ten = T1.*normal(:,1)+T2.*normal(:,2)+T3.*normal(:,3);

     G(:,i) = S_ds;
     %Gten(:,i) = S_ten;
       
   
     
  %  end
  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
  end
  
  
function bs = Compress(g,rerr,savefn)


  G = @(rs,cs)Fn(g,rs,cs);
  
   %X is 3xN point cloud
   X = (g.centers);

 
    [bs p q] = hmcc_hd(X,X);

    t = tic();
    n = length(p);
    eBfro = hm('EstBfroLowBndCols',G,p,q,1:n,1:n);
    savefn_tmp = [savefn '.tmp'];
    hm('Compress',bs,p,q,G,1e-5*rerr,savefn_tmp,'Bfro',eBfro);
    et = toc(t);
    [id hnnz] = hm_mvp('init',savefn_tmp);
    in = hm('HmatInfo',savefn_tmp);
    cf = in.m*in.n/hnnz;
    Bfro = sqrt(hm_mvp('fronorm2',id));
    [ce re] = hm('CondEstFro',id,struct('reltol',0.05));
    hm_mvp('cleanup',id);
    fprintf(1,'et: %1.1f  cf: %1.1f  ||B||_F: %1.1e %1.1e  ce: %1.1e\n',...
	    et,cf,eBfro,Bfro,ce);
    
    if (false)
      t = tic();
      hm('Compress',bs,p,q,G,rerr,savefn,...
	 'old_hmat_fn',savefn_tmp,'tol_method','BinvFro','BinvFro',ce/Bfro);
      [id hnnz] = hm_mvp('init',savefn);
      hm_mvp('cleanup',id);
      cf = in.m*in.n/hnnz;
      fprintf(1,'et: %1.1f  cf: %1.1f\n',toc(t),cf);
    else
      system(sprintf('mv %s %s',savefn_tmp,savefn));
    end
 

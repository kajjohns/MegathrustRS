function varargout = CalcShear_edgestress(fn,varargin)

[varargout{1:nargout}] = feval(fn,varargin{:});

function Make(patches,rerr,savefn)
  g = CalcShear_edgestress('Init',patches);
  CalcShear_edgestress('Compress',g,rerr,savefn);
  
function g = Init(patches)
    
    %EdgeStress
    %The form of m:   
    %   
    %     m(1) = Depth of updip end    
    %     m(2) = Horizontal position of updip end    
    %     m(3) = Length downdip   
    %     m(4) = Dip (degrees)   
    %     m(5) = Slip   

      Length = sqrt( (patches(:,1)-patches(:,3)).^2 + (patches(:,2)-patches(:,4)).^2 );
      Dip = 180/pi*atan(-(patches(:,2)-patches(:,4))./(patches(:,1)-patches(:,3)));

       m(1,:) = -patches(:,2)'; %negative to convert to depth
       m(2,:) = patches(:,1)';
       m(3,:) = Length';
       m(4,:) = Dip';
       m(5,:) = ones(size(m(1,:)));

      
        
        %convert km to m
        m(1,:)=1000*m(1,:);
        m(2,:)=1000*m(2,:);
        m(3,:)=1000*m(3,:);

     
        %calculate patch centers
        centers_x = (patches(:,1)+patches(:,3))/2;
        centers_z = (patches(:,2)+patches(:,4))/2;
        
      
       
        %calculate patch normals
        
        strikevec=[zeros(size(centers_x)) ones(size(centers_x)) zeros(size(centers_x))];
        dipvec=[patches(:,3) 0*patches(:,3) patches(:,4)]-[centers_x 0*centers_x centers_z];
        vecnorm=sqrt(dipvec(:,1).^2+dipvec(:,2).^2+dipvec(:,3).^2);
        vecnorm=repmat(vecnorm,1,3);

        %need to be very careful about signs
        %if dip<90, dipvec should point up
        %if dip>90, dipvec should point down
        index = Dip<90;
        dipvec=dipvec./vecnorm;  %change sign so that dipvec points up
        dipvec(index,:) = -dipvec(index,:);

        normal=cross(strikevec,dipvec,2);

  


        g.m = m;
        g.strikevec = strikevec;
        g.centers = 1000*[centers_x';0*centers_x';centers_z']; %convert to meters
        g.dipvec = dipvec;
        g.normal = normal;

        g.mu = 3e10;  %everything is in m
        g.nu = 0.25;
        
         
     
  
function G = Fn(g,rs,cs)



  G = zeros(length(rs),length(cs));
  for(i = 1:length(cs))
    k = cs(i);
    
     %m=[g.pm(k,:) 0 1 0]';
     %xloc = g.centers;
     
     m = g.m(:,k);
     [Sxx,Sxz,Szz]=EdgeStress(m,g.centers(1,rs),g.centers(3,rs),g.nu,g.mu);




     normal = g.normal(rs,:);
%     %calculate tractions
    T1=Sxx'.*normal(:,1)+Sxz'.*normal(:,3);
    %T2=S(2,:)'.*normal(:,1)+S(4,:)'.*normal(:,2)+S(5,:)'.*normal(:,3);
    T2 = zeros(size(T1));
    T3=Sxz'.*normal(:,1)+Szz'.*normal(:,3);

    
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
 

function [topx_interface,topz_interface,botx_interface,botz_interface] = make_geometry_interface(faultdip_trench,x_trench,x_bottom,faultdip_bottom,wedge_bot,pL)



%make subduction interface


%top
  %represent curved surface as cubic
     %need to shift dislocation so that top and bottom lines up with patch end
       fault_bot_x = x_bottom;
    
         
       %get coefficients of cubic 
       top_fault = 0;
       b=[tan(-faultdip_bottom*pi/180);tan(-faultdip_trench*pi/180);-top_fault;-wedge_bot];
       A=[ [3*fault_bot_x^2 2*fault_bot_x 1 0];[3*x_trench^2 2*x_trench 1 0];[x_trench^3 x_trench^2 x_trench 1];[fault_bot_x^3 fault_bot_x^2 fault_bot_x 1] ];
       c=inv(A)*b;
       
       
       depth = 0;
       dist = 0;
       x_coord = x_trench;
       cnt=1;
       z_coord = depth;
       while abs(depth)<wedge_bot
           
           cnt = cnt+1;
           dip = atan(c(1)*3*dist^2 + c(2)*2*dist + c(3));
           
           dx = pL*cos(dip);
           dist = dist + dx;
           z_coord(cnt) = c(1)*dist^3 + c(2)*dist^2 + c(3)*dist + c(4);
           x_coord(cnt) = x_trench + dist;
           
           depth = z_coord(end);
           
           
       end
      
       
topx_interface = x_coord(1:end-1);
topz_interface = z_coord(1:end-1);
botx_interface = x_coord(2:end);
botz_interface = z_coord(2:end);



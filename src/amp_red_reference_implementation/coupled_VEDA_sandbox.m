function [final, time, Force, h_0, h_0_dot, u0]=coupled_VEDA(meshnum,time_steps,Z,Amp,TR,MR,f1,Ham,z_0,nu,E_inf,E_0,tau)

    E_inf_s= E_inf/(1-nu*nu);
    E_0_s= E_0/(1-nu*nu);

    h0_r_t=zeros(meshnum,time_steps);
    h = zeros(meshnum, time_steps);
    
    u_dot_drk = zeros(meshnum, time_steps);
    u_dot_bahram = zeros(meshnum, time_steps);
    u_dot = zeros(meshnum, time_steps);
    
    u = zeros(meshnum, time_steps);

    %% Func def
    profilefun=@(r) r.^2/2/TR; % tip profile
    
%Way Bahram coded it
    Pressure = @(h,r) Ham/6/pi./h.^3 .* ((z_0 ./ h).^6 - 1);     %==Lennard-Jones adhesion==        
    Preprime = @(h,r) Ham/2/pi./h.^4 .* (1-3.*(z_0 ./ h).^6);     %3-9 Lennard Jones derivitive

% this is the way I coded it originally in fortran. tiny differences in roundoff errors.
    %Pressure_drk = @(h,r)  (Ham  ./ (6. * pi * h.^3)) .* ( (z_0^6) ./(h.^6) - 1.)
    %Preprime_drk = @(h,r)  ( 0.5 * Ham *  (-3. * z_0^6 + h.^6)   ./ (pi * h.^10))
    
    %% radial disc
    radius = (linspace(0, MR, meshnum))';
    dr = radius(2);
    radius = radius + dr/2;

    %% time disc
    h0min=Z-Amp;
    h0max=5e-9;
    
    if ((h0max-Z)/Amp > 1)
      %drk 7/31/2020 handle the case when the tip is always within h0max of the sample for the entire cycle
      start_time = 0
    else     
      start_time = 1/ (2*pi*f1)*acos((h0max-Z)/Amp);
    end 
    
    end_time = (1/f1)-start_time;
    time=linspace(start_time,end_time,time_steps);
    h_0=    Z+Amp*        cos(2*pi*f1*time);
    h_0_dot= -Amp*2*pi*f1*sin(2*pi()*f1*time);
    q=Amp*cos(2*pi*f1*time);
    dt=time(2)-time(1);

    for i=1:meshnum
        for j=1:time_steps
            h0_r_t(i,j)=h_0(j)+ profilefun(radius(i));
        end
    end

    %% Kernel calc

    kernel = zeros(meshnum, meshnum);
    for r = 1:meshnum
        for s = 1:r
          if s==r
            Iminus = (dr/pi - dr^2/4/pi/radius(r)) ...
                     * (1 + log(16*radius(r)^2/(radius(r)*dr - dr^2/4)) );
            Iplus  = dr/pi*log(16*(radius(r)+dr/2)^2/(radius(r)*dr + dr^2/4))...
                     + 4*radius(r)/pi*log((2*radius(r)+dr)/(2*radius(r)+dr/2));
            kernel(r,s) = (Iminus+Iplus)/radius(r)/dr;
            break
          else
            kernel(r,s) = 4/pi/radius(r) * ellipke((radius(s)/radius(r))^2);
            kernel(s,r) = kernel(r,s);
          end
        end
    end


    %% Main Computation part
    
%    B_drk = zeros(meshnum, 1);
    
    for t=1:time_steps-1

        h(:,t)=h0_r_t(:,t)-u(:,t);

        A=(radius.*Preprime(h(:,t)).*kernel*dr)'/E_0_s;                        
        
        B=kernel*(radius.*Pressure(h(:,t),radius))* dr/E_inf_s/tau;
         
         %  p = (radius.*Pressure(h(:,t),radius))/E_inf_s/tau;
         %  for j = 1:meshnum
         %      B_drk(j,1) = trapz(radius, kernel(j,:)'.* p );
         %  end
%           %B is -uinfinity/tau 
%          

        
        AB=A*ones(meshnum,1);
        
                  
    %       for j = 1:meshnum
    %           AB_drk(j,1) = trapz( A(j,:) );
    %       end
         
        
        AA=A-eye(size(A));
 
   %       AA_drk = A;
   %       AA_drk(:,1) = AA_drk(:,1) /2;
   %       AA_drk(:,end) = AA_drk(:,end) /2;
   %       AA_drk=(AA_drk-eye(size(A)));
%         
        b = u(:,t)/tau + B + h_0_dot(t)*AB;
        
  %       b_drk = u(:,t)/tau + B_drk + h_0_dot(t)*AB_drk;
%         
 %        u_dot_drk(:,t) = AA_drk \ b_drk;        
        u_dot_bahram(:,t)=(AA\b);
        
        u_dot(:,t) = u_dot_bahram(:,t);
        
        u(:,t+1)=u(:,t)+ dt*u_dot(:,t);  
    end

    %I'm doing trapezoidal for this whereas he is doing a sum.  is that it?
    Force_r= ones(1,meshnum)*(radius.*Pressure(h(:,1:time_steps-1),radius))*2*pi*dr;

%    Force_t= trapz(radius, radius.* Pressure(h(:,1:time_steps-1),radius) *2*pi);
    
    Force = Force_r;
    
    u0 = u(1,:);
    
    dissipation = -trapz(time(1:length(time)-1),Force.*h_0_dot(1:length(time)-1));    % unit N.nm (nJ)
    virial = f1*trapz(time(1:length(time)-1),Force.*q(1:length(time)-1));         % unit N.nm (nJ)
    
    final=[dissipation,virial];

%    close all
%    figure(1);
%    plot(h_0(1:length(time)/2),Force(1:length(time)/2)*1e9,'b','LineWidth',1);
%    hold on;
%    plot(h_0(length(time)/2:length(time)-1),Force(length(time)/2:length(time)-1)*1e9,'r','LineWidth',1);
%    set(get(gca,'XLabel'),'String','Separation (h0 nm)','LineWidth',1)
%    set(get(gca,'YLabel'),'String','Force (nN)');
    
    

end

function final=u_Fourier_VEDA(m,meshnum,time_steps,Z,Amp,TR,MR,f1,Ham,z_0,nu,E_inf,E_0,tau)

    E_inf_s= E_inf/(1-nu*nu);
    E_0_s= E_0/(1-nu*nu);

    h0_r_t=zeros(meshnum,time_steps);
    cos_r=zeros(meshnum,m+1);
    a=zeros(m+1,time_steps);
    a_dot=zeros(m+1,time_steps-1);
    h = zeros(meshnum, time_steps);


    %% Func def
    profilefun=@(r) r.^2/2/TR; % tip profile
    Pressure = @(h) Ham/6/pi./h.^3 .* ((z_0 ./ h).^6 - 1);     %==Lennard-Jones adhesion==
    Preprime = @(h) Ham/2/pi./h.^4 .* (1-3.*(z_0 ./ h).^6);     %3-9 Lennard Jones derivitive

    %% radial disc
    radius = (linspace(0, MR, meshnum))';
    dr = radius(2);
    radius = radius + dr/2;

    %% time disc
    h0max=5e-9;

    if ((h0max-Z)/Amp > 1)
      %drk 7/31/2020 handle the case when the tip is always within h0max of the sample for the entire cycle
      start_time = 0
    else     
      start_time = 1/ (2*pi*f1)*acos((h0max-Z)/Amp);
    end 

    end_time = (1/f1)-start_time;
    time=linspace(start_time,end_time,time_steps);
    h_0=Z+Amp*cos(2*pi*f1*time);
    h_0_dot=-Amp*2*pi*f1*sin(2*pi()*f1*time);
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
    
    for i=1:m+1
        cos_r(:,i)=(cos((i-1)*pi*radius/MR));
    end

    gamma =radius.*(kernel*cos_r*dr)*dr;
    Coef1=MR*E_0_s*ones(m+1,1)/2;
    Coef1(1)=MR*E_0_s;
    Coef2=diag(Coef1);

    for t=1:time_steps-1

        h(:,t)=h0_r_t(:,t)-cos_r*a(:,t);

        A=gamma'*(cos_r.*Preprime(h(:,t)));
        B=gamma'*(Pressure(h(:,t)));

        AA=(a(:,t).*Coef1/tau)+(E_0_s/tau/E_inf_s*B)+(A(:,1)*h_0_dot(t));
        A=A-Coef2;

        a_dot(:,t)=A\AA;
        a(:,t+1)=a(:,t)+dt*a_dot(:,t);
    end

    Force= ones(1,meshnum)*(radius.*Pressure(h(:,1:time_steps-1)))*2*pi*dr;
    dissipation = -trapz(time(1:length(time)-1),Force.*h_0_dot(1:length(time)-1));    % unit N.nm (nJ)
    virial = f1*trapz(time(1:length(time)-1),Force.*q(1:length(time)-1));         % unit N.nm (nJ)
    
    final=[dissipation,virial];
    
    close all
    figure(1);
    plot(h_0(1:length(time)/2),Force(1:length(time)/2)*1e9,'b','LineWidth',1);
    hold on;
    plot(h_0(length(time)/2:length(time)-1),Force(length(time)/2:length(time)-1)*1e9,'r','LineWidth',1);
    set(get(gca,'XLabel'),'String','Separation(h0) (nm)','LineWidth',1)
    set(get(gca,'YLabel'),'String','Force (nN)');

end

function [F]=F_BM_AEM_VEDA(casenum, cyc_numb, Ar1, h0,time_step, profileype, num_time_steps, R, maximum_radius, meshnum, A, z0, nu, tau_j, E0, Einf , Z_iter_num, A2_iter_num)

    address=strcat(pwd,'\');  
    arstep=1;
    
    CC_0=1/Einf;
    CC_j=CC_0-1/E0;
    ncct=1;
    ODEsolver=1;
    
    tic

    %%
    Pressure = @(h,r) A/6/pi./h.^3 .* ((z0 ./ h).^6 - 1);     %==Lennard-Jones adhesion==
    Preprime = @(h,r) A/2/pi./h.^4 .* (1-3.*(z0 ./ h).^6);     %6-12 Lennard Jones derivitive

    %% Radial discretization and tip profile
    radius = linspace(0, maximum_radius, meshnum);
    dr = radius(2);
    radius = radius + dr/2;
    if profileype==1 
        profilefun=@(r) r.^2/2/R;
    elseif profileype==2
        profilefun = @(r) roundcone(r, 60, R)-R;
    elseif profileype==3
        profilefun = @(r) -sqrt(R^2-r.^2)+R;
    elseif profileype==4
        profilefun = @(r) (-1.329920e-11)*r.^6 + (-2.161565e-6)*r.^4 + (3.175362e-2)*r.^2;
    end
    profile_ = profilefun(radius');


    %% INITIALISE AND CALCULATE NEEEDED VARIABLES
    u = zeros(meshnum, num_time_steps);
    ul = zeros(ncct* meshnum, num_time_steps);
    uldot = zeros(ncct* meshnum, num_time_steps);
    h = zeros(meshnum, num_time_steps);
    F = zeros(1, num_time_steps);
    Jacob0= zeros(meshnum*ncct, meshnum*ncct);
    kernel = zeros(meshnum, meshnum);
    conmat=zeros(meshnum,meshnum*ncct);
    dtr= zeros(num_time_steps*20, 6);
    b1= zeros(meshnum*ncct, 1);
    b2= zeros(meshnum*ncct, 1);
    b3= zeros(meshnum*ncct, 1);
    h0dot= zeros(1, num_time_steps);

    %% COMPUTE THE KERNEL/Jacob/b matrices
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

    for i=1:meshnum*ncct		% jacob and b matrix definition 
        for j = 1:meshnum*ncct 
            if mod(i,ncct)==0
                ii=floor(i/ncct);
                iii=ncct;
            elseif mod(i,ncct)~=0
                ii=floor(i/ncct)+1;
                iii=mod(i,ncct);
            end
            if mod(j,ncct)==0
                jj=floor(j/ncct);
            elseif mod(j,ncct)~=0
                jj=floor(j/ncct)+1;
            end

            Jacob0(i,j)=kernel(ii,jj)*radius(jj)*(1-nu^2)*dr*(CC_0-CC_j(iii));
        end

        if mod(i,ncct)==0
            iiii=ncct;
        elseif mod(i,ncct)~=0
            iiii=mod(i,ncct);
        end
        b1(i)=1/tau_j(iiii);
        b2(i)=CC_0*(1-nu^2)*dr/tau_j(iiii);
        b3(i)=(CC_0-CC_j(iiii))*(1-nu^2)*dr;

    end

    for ii=1:meshnum
        conmat(ii,(ii-1)*ncct+1:ii*ncct)=1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % IC
    h(:,1) = h0(1) + profile_;
    v=0;

    %% MAIN SIMULATION
    for t = 2:length(h0)		% loop for each time

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %clearvars kpr1 kpr2 kpr kppr1 kppr2 kppr Jacob1 Jacob2 Jacob3 Jacob b;
        clear kpr1 kpr2 kpr kppr1 kppr2 kppr Jacob1 Jacob2 Jacob3 Jacob b;

        kpr1=transpose(kernel*(Pressure(h(:,t-1),radius(:)).*radius'));
        kpr2=repmat(kpr1,ncct,1);
        kpr=reshape(kpr2,meshnum*ncct,1);

        kppr1=transpose(kernel*(Preprime(h(:,t-1),radius(:)).*radius'));
        kppr2=repmat(kppr1,ncct,1);
        kppr=reshape(kppr2,meshnum*ncct,1);

        Jacob2=transpose(Preprime(h(:,t-1),radius(:)));
        Jacob3=repmat(Jacob2,ncct,1);
        Jacob1=Jacob0.*repmat(reshape(Jacob3,1,meshnum*ncct),meshnum*ncct,1);

        Jacob=Jacob1-eye(size(Jacob1));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%        if and(v>0,exist('Ar_guess')>0) 
%             dtr(v,:)=[Ar_guess Ar(arstep) dzold h0(t)];
%            [arstep Ar_guess Ar1(arstep) dzold h0(t)]
%        elseif and(v>0,exist('Ar_guess')==0)
%             dtr(v,:)=[0 Ar1(arstep) 0 h0(t)];
%            [arstep 0 Ar1(arstep) 0 h0(t)]
%        end  
%        v=v+1;

        %%%%%%%%%%%%%%%%%%%%%% Mat. Inversion method
        
        h0dot(t-1)=(h0(t)-h0(t-1))/time_step(t-1);
        b = ul(:,t-1).*b1 + kpr.*b2 + h0dot(t-1)*kppr.*b3;
        uldot(:,t-1)=(Jacob\b);

        if and(ODEsolver==5,t>5)
            ul(:,t)=ul(:,t-1)+ time_step(t-1)/720*(1901*uldot(:,t-1)-2774*uldot(:,t-2)+2616*uldot(:,t-3)-1274*uldot(:,t-4)+251*uldot(:,t-5)) ;
        elseif or(and(ODEsolver==4,t>4),and(ODEsolver>4,t==5))
            ul(:,t)=ul(:,t-1)+ time_step(t-1)/24*(55*uldot(:,t-1)-59*uldot(:,t-2)+37*uldot(:,t-3)-9*uldot(:,t-4)) ;
        elseif or(and(ODEsolver==3,t>3),and(ODEsolver>3,t==4))
            ul(:,t)=ul(:,t-1)+ time_step(t-1)/12*(23*uldot(:,t-1)-16*uldot(:,t-2)+5*uldot(:,t-3)) ;
        elseif or(and(ODEsolver==2,t>2),and(ODEsolver>2,t==3))
            ul(:,t)=ul(:,t-1)+ time_step(t-1)/2*(3*uldot(:,t-1)-uldot(:,t-2)) ;
        elseif or(and(ODEsolver==1,t>1),and(ODEsolver>1,t==2)) 
            ul(:,t)=ul(:,t-1)+ time_step(t-1)*uldot(:,t-1);  
        end

        u(:,t)= conmat*ul(:,t);
        h(:,t) = h0(t) + profile_ - u(:,t);
        F(t) = 2 * pi * sum(Pressure(h(:,t),radius(:)) .* radius(:)) * dr; % (N)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    %% 
%     close all
%     figure(1);
%     plot(h0,F*1e9,'b','LineWidth',1);
%     title('Force v.s. Separation');
%     set(get(gca,'XLabel'),'String','Separation (h0)  (nm)','LineWidth',1);
%     set(get(gca,'YLabel'),'String','Force (nN)');
%     xlim([-inf 5]) ;
%     
%     filename4=[address 'results\' num2str(casenum) '_' num2str(Ar1) '_' num2str(Z_iter_num) '_' num2str(A2_iter_num) '_' num2str(cyc_numb) '.jpg'];
%     ii=1;

    output1=zeros(1000,2);
    for i=1:max(1,floor(num_time_steps/1000)):length(F)
        output1(ii,:)=[h0(i), F(i)*1e9];
        ii=ii+1;
    end
    
%    print('-f1','-dtiff','-r200',filename4)
 
end



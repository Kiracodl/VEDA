function final=zsingle_VEDA(casenum,Ar,Af,profileype,A, z0, nu, ncct, tau_j, CC_j, CC_0...
    ,ODEsolver,nstep,R,maximum_radius,meshnum, tol,k,Q,per)
    address=strcat(pwd,'\');  

    amp=Ar*Af;
    sign=1;
    dz=1;       %nm 
    Zguess=zeros(size(amp));
    
    Zguess(1)=amp(1)+dz;
    
    results=zeros(length(Ar),10);
    Ar_guess=NaN
    
    arstep=1;
    jjj=1; 
    while (arstep==1)    
      
        %% ==Lennard-Jones adhesion==
        Pressure = @(h,r) A/6/pi./h.^3 .* ((z0 ./ h).^6 - 1);     %==Lennard-Jones adhesion==
        Preprime = @(h,r) A/2/pi./h.^4 .* (1-3.*(z0 ./ h).^6);     %3-9 Lennard Jones derivitive

        %% Waveform

        h0min_tapping=Zguess(arstep)-amp(arstep);
        if (Zguess(arstep)+amp(arstep))<5
            h0max_tapping=Zguess(arstep)+amp(arstep);
        else
            h0max_tapping=5; %nm
        end

        start_time = per/360*acosd((h0max_tapping-h0min_tapping-amp(arstep))/amp(arstep));
        time=linspace(start_time,per-start_time,nstep);
        time1=time(1:floor(nstep/2));
        time2=time(floor(nstep/2)+1:end);    
        ApproachLength=length (time1);
        RetractLength=length (time2);
        num_time_steps = length(time);
        time_step = diff(time);
        h0 = (amp(arstep)+h0min_tapping)+ amp(arstep)*cosd(time/per*360);
        h0dot= -2*pi*amp(arstep)/per*sind(360*time/per);
        h0_step = -diff(h0);
        q_tapp=h0-Zguess(arstep);
        qdot_tapp=h0dot;

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
        profile = profilefun(radius');


        %% INITIALISE AND CALCULATE NEEEDED VARIABLES
        u = zeros(meshnum, num_time_steps);
        ul = zeros(ncct* meshnum, num_time_steps);
        uldot = zeros(ncct* meshnum, num_time_steps);
        h = zeros(meshnum, num_time_steps);
        F = zeros(1, num_time_steps);
        b= zeros(meshnum*ncct, 1);
        Jacob= zeros(meshnum*ncct, meshnum*ncct);
        Jacob0= zeros(meshnum*ncct, meshnum*ncct);
        Jacob1= zeros(meshnum*ncct, meshnum*ncct);
        q= zeros( num_time_steps,1);
        %dtr= zeros(num_time_steps*20, 7);  %this is just iteration statistics?
        kpr= zeros(meshnum, 1);
        kppr= zeros(meshnum, 1);
        conmat=zeros(meshnum,meshnum*ncct);
        uu = zeros(meshnum, num_time_steps);
        uul = zeros(ncct* meshnum, num_time_steps);
        b1= zeros(meshnum*ncct, 1);
        b2= zeros(meshnum*ncct, 1);
        b3= zeros(meshnum*ncct, 1);

        %% COMPUTE THE KERNEL/Jacob/b matrices
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
        h(:,1) = h0(1) + profile;
        %v=0;

        %% MAIN SIMULATION
        for t = 2:num_time_steps		% loop for each time

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%            clearvars kpr1 kpr2 kpr kppr1 kppr2 kppr Jacob1 Jacob2 Jacob3 Jacob b;
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

            %just iteration statistics? 
            
            %if and(v>0,exist('Ar_guess')>0) 
            %    dtr(v,:)=[Ar_guess Ar(arstep) dz dzold Zguess(arstep) q_tapp(t) h0(t)];
            %    [arstep Ar_guess Ar(arstep) dz dzold Zguess(arstep) q_tapp(t) h0(t)]
            %elseif and(v>0,exist('Ar_guess')==0)
            %    dtr(v,:)=[0 Ar(arstep) dz 0 Zguess(arstep) q_tapp(t) h0(t)];
            %    [arstep 0 Ar(arstep) dz 0 Zguess(arstep) q_tapp(t) h0(t)]
            %end  
            %v=v+1;

            %%%%%%%%%%%%%%%%%%%%%% Mat. Inversion method

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
            h(:,t) = h0(t) + profile - u(:,t);
            F(t) = 2 * pi * sum(Pressure(h(:,t),radius(:)) .* radius(:)) * dr; % (N)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end



    %%    

        dissipation = -trapz(time,F.*qdot_tapp);    % unit N.nm (nJ)
        virial=1/per*trapz(time,F.*q_tapp);         % unit N.nm (nJ)

        Ar_guess=1/Q/sqrt(((1/Q)+(dissipation/pi/k/amp(arstep)/amp(arstep)))^2+...
            (2*virial/k/amp(arstep)/amp(arstep))^2);
        Phase_lag(arstep)=atan(((1/Q)+(dissipation/pi/k/amp(arstep)/amp(arstep)))...
            /(2*virial/k/amp(arstep)/amp(arstep)))*180/pi;

        ['Z single, Zguess: ' num2str(Zguess(arstep)), ' error: ' num2str(abs(Ar_guess-Ar(arstep)))]
    
        signold=sign;
        dzold=dz;

        if (abs(Ar_guess-Ar(arstep))>tol)
            if (Ar_guess>Ar(arstep))
                sign=1;
                if sign*signold>0
                    dz=dzold;
                elseif sign*signold<0
                    dz=dzold/2;
                end    
                Zguess(arstep)=Zguess(arstep)-dz;
            elseif (Ar_guess<Ar(arstep))
                sign=-1;
                if sign*signold>0
                    dz=dzold;
                elseif sign*signold<0
                    dz=dzold/2;
                end    
                Zguess(arstep)=Zguess(arstep)+dz;
            end
        elseif abs(Ar_guess-Ar(arstep))<tol
            if arstep<length(Ar)
                Zguess(arstep+1)=amp(arstep+1)+Zguess(arstep)-amp(arstep);
                dz=5*dz;
            end
            %%
         %   close all
         %   figure(1);
         figure
            plot(h0(1:ApproachLength),F(1:ApproachLength)*1e9,'b','LineWidth',1);
            hold on;
            plot(h0(ApproachLength:ApproachLength+ RetractLength),F(ApproachLength:ApproachLength+ RetractLength)*1e9,'r','LineWidth',1);
            title('Force vs. Separation');
            set(get(gca,'XLabel'),'String','Separation (h0 nm)','LineWidth',1)
            set(get(gca,'YLabel'),'String','Force (nN)');
            %%
            filename4=[address 'results\' num2str(casenum) '_' num2str(Ar(arstep)) '_tapping.jpg'];
           %saveas(figure(1),filename4)

            phase=atan(((1/Q)+(dissipation/pi/k/amp(arstep)/amp(arstep)))/(-2*virial/k/amp(arstep)/amp(arstep)))*180/pi;
            if phase<0
                phase=phase+180;
            end

            results(jjj,:)=[arstep Ar(arstep) Ar_guess Zguess(arstep)  amp(arstep) h0min_tapping h0max_tapping dissipation*6.242e9 virial*6.242e9 phase]; 
            jjj=jjj+1;
            arstep=arstep+1;    
        end
    end
    
    final=[Zguess(arstep-1),phase,dissipation,virial];
    
end
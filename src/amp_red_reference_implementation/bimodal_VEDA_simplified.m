%opengl('save', 'software');
   
clear
%clc
%close all

ex_opt=11

Z_single_guess_override = []; % for debugging/testing speed up iteration by pre-specifying the initial value
E0 = []; 
Einf = [];
CC_0 = [];
CC_j = [];

if (ex_opt==11)
  
  k2=1103.2 * 1e-9;          % Cantilever 2nd equivalent stiffness (N/nm)	
  k=28.09 * 1e-9;        	% Cantilever equivalent stiffness (N/nm)	

  per=1/(277.7778e3);      	% Time Period of the first mode (sec) 
  per2=1/(0.1736111E+07);        % Time Period of the 2nd mode (sec)
  
  E0=2e-9;               	% N/(nm)^2 = 1e-18 GPa 
                          % typo? Meant 1e-18 Pa? not GPa?
                          % 2 GPa = E1
  Einf=0.2e-9;           	% N/(nm)^2 = 1e-18 GPa
                          % 0.2 GPa
                          %E2 = (Einf E0 ) / (E0 - Einf) = (0.2 * 2)/ (2-0.2)=0.222222 GPa
  
  %in original version, F_BM_AEM_VEDA_markedup is using E0 & Einf, but  zsingle_VEDA_markedup is using CC_0 & CC_j.  
  
  Sol_type=1;             % 2: Simplified, 2: Coupled
  nstep=1.00E+04;         % Number of time steps	
  R=10;                   % Radius of the tip (nm)	
  maximum_radius=15;   	% computation domain (nm)
  meshnum=70;         	% number of radial discretization
  
  Q=429.5;            	% Quality factor
  Q2=Q;            % Quality factor of 2nd mode
  
 
  Armax=0.9;          	% max Ar range
  Armin=0.1;              % min Ar range
  npoints=9;        		% number of the Ar selected on the range
  
  Af=75;              	% free Amp (nm)
  Af2=3.5;                  % free Amp of 2nd mode(nm)
  
  A=1.00E-11;         	% Hamaker's constant (N/nm)	
  z0=0.3;             	% Equilibrium position (nm)
  nu=0.49;            	% Poisson's ratio
  tau_j=1.00E-07;       	% Creep time (Sec)	
    
  tol=0.01;           	% predefined tolerance for Ar predictions
  %basis_func_number=25; 	% number of basis functions, only needed for simplified model
  n_c=4;                 % max number of the 1st mode oscilations to study the convergence
  
  Z_single_guess_override = [ 6.3707896E+001   5.5253437E+001   4.7260971E+001  3.9485657E+001  3.1875355E+001  2.4412828E+001  1.7114560E+001   1.0026368E+001   3.2763681E+000];
  A2_override =             [ 0.3214267E+01    0.3118547E+01    0.3048926E+01    0.2982637E+01   0.2915743E+01  0.2848869E+01   0.2776804E+01  0.2682047E+01   0.2522870E+01];

elseif (ex_opt==12)
  %exampled as delivered
  k2=5.0982e-07;          % Cantilever 2nd equivalent stiffness (N/nm)	
  Q2=599.7000;            % Quality factor of 2nd mode
  per2=6.3000e-07;        % Time Period of the 2nd mode (sec)
  Af2=5;                  % free Amp of 2nd mode(nm)
  CC_j=1.5000e+09;        % creep compliance terms_j GPa
  CC_0=2.0000e+09;        % creep compliance terms_0 GPa
  Sol_type=1;             % 2: Simplified, 2: Coupled
  nstep=1.00E+04;         % Number of time steps	
  R=10;                   % Radius of the tip (nm)	
  maximum_radius=15;   	% computation domain (nm)
  meshnum=70;         	% number of radial discretization
  k=2.809E-08;        	% Cantilever equivalent stiffness (N/nm)	
  Q=429.5;            	% Quality factor
  per=0.0000036;      	% Time Period of the first mode (sec)
  Armax=0.5;          	% max Ar range
  Armin=0.5;              % min Ar range
  npoints=1;        		% number of the Ar selected on the range
  Af=75;              	% free Amp (nm)
  A=1.00E-10;         	% Hamaker's constant (N/nm)	
  z0=0.4;             	% Equilibrium position (nm)
  nu=0.49;            	% Poisson's ratio
  tau_j=1.00E-07;       	% Creep time (Sec)	
  %E0=2e-9;               	% N/(nm)^2 = 1e-18 GPa !material double specified!
  %Einf=0.2e-9;           	% N/(nm)^2 = 1e-18 GPa
  tol=0.01;           	% predefined tolerance for Ar predictions
  %basis_func_number=25; 	% number of basis functions, only needed for simplified model
  n_c=11;                 % max number of the 1st mode oscilations to study the convergence
elseif (ex_opt == 13)
  %stripped down to run faster
  k2=5.0982e-07;          % Cantilever 2nd equivalent stiffness (N/nm)	
  Q2=599.7000;            % Quality factor of 2nd mode
  per2=6.3000e-07;        % Time Period of the 2nd mode (sec)
  Af2=5;                  % free Amp of 2nd mode(nm)
  CC_j=1.5000e+09;        % creep compliance terms_j GPa
  CC_0=2.0000e+09;        % creep compliance terms_0 GPa
  Sol_type=1;             % 2: Simplified, 2: Coupled
  nstep=1000;         % Number of time steps	
  R=10;                   % Radius of the tip (nm)	
  maximum_radius=15;   	% computation domain (nm)
  meshnum=25;         	% number of radial discretization
  k=2.809E-08;        	% Cantilever equivalent stiffness (N/nm)	
  Q=429.5;            	% Quality factor
  per=0.0000036;      	% Time Period of the first mode (sec)
  Armax=0.5;          	% max Ar range
  Armin=0.5;              % min Ar range
  npoints=1;        		% number of the Ar selected on the range
  Af=75;              	% free Amp (nm)
  A=1.00E-10;         	% Hamaker's constant (N/nm)	
  z0=0.4;             	% Equilibrium position (nm)
  nu=0.49;            	% Poisson's ratio
  tau_j=1.00E-07;       	% Creep time (Sec)	
  %E0=2e-9;               	% N/(nm)^2 = 1e-18 GPa !material double specified
  %Einf=0.2e-9;           	% N/(nm)^2 = 1e-18 GPa
  tol=0.01;           	% predefined tolerance for Ar predictions
  %basis_func_number=25; 	% number of basis functions, only needed for simplified model
  n_c=3;                 % max number of the 1st mode oscilations to study the convergence  
  
  Z_single_guess_override = [33.75];
end
 
%per Bahram e-mail 8/14/2020
%The zsingle_VEDA is using the creep compliance data, and the F_BM_AEM_VEDA is using the stress relaxation data!!!!
if (isempty(E0))
  E0=1/(CC_0-CC_j);                       % N/(nm)^2 = 1e-18 GPa
  Einf=1/CC_0;                    % N/(nm)^2 = 1e-18 GPa
elseif (isempty(CC_0))
  CC_0 = 1/Einf;
  CC_j = CC_0 - 1/E0
end
 
casenum=1;
profileype=1;
ncct=1;
%pl=0;
ODEsolver=1;
address = strcat(pwd,'\results\');  

dA2_initial=Af2/10;
dZ_initial=Af2/30;

tol_Ar_1=0.025;
tol_Ar_2=0.025;
%AR2_tol=tol_Ar_2;
%PL2_tol=0.1;
%conv_lim=20; % number of cycles with less than Ar2 and PL2 tolerance deviation
Ar1_d=fliplr(linspace(Armin,Armax,npoints));  
max_height=5;

%%
Ar2m_Ext=zeros(npoints,1000);
Ar2m_Ext_ave =zeros(npoints,1000);
Phase_lag2m_Ext =zeros(npoints,1000);
Phase_lag2m_Ext_ave =zeros(npoints,1000);
Phase_lag2m_Ext_cyc =zeros(npoints,1000);
Ar2m_Ext_cyc=zeros(npoints,1000);
Ar2=zeros(npoints,1000);
Ar1=zeros(npoints,1000);
A2=zeros(npoints,1000);
Zguess=zeros(npoints,1000);
Phase_lag_sm_2=zeros(1,npoints);
Phase_lag_sm_1=zeros(1,npoints);
Ets_sm_1=zeros(1,npoints);
Vts_sm_1=zeros(1,npoints);
z_sm_1=zeros(1,npoints);
z_sm_2=zeros(1,npoints);
Phase_lag2=zeros(1,npoints);
Phase_lag1=zeros(1,npoints);   
Z_final=zeros(1,npoints);
A2_final=zeros(1,npoints);
Ar2_final=zeros(1,npoints);
Phase_lag1_final=zeros(1,npoints);
Phase_lag2_final=zeros(1,npoints);
E_ts_2_final = zeros(1,npoints);
V_ts_2_final = zeros(1,npoints);
E_ts_1_final = zeros(1,npoints);
V_ts_1_final = zeros(1,npoints);
errorfinal_1=zeros(1,npoints);
errorfinal_2=zeros(1,npoints);
it_i=zeros(1,npoints);
it_j=zeros(1,npoints);
ave_diss_2 = zeros(npoints,n_c);
ave_virial_2= zeros(npoints,n_c);
ave_Ar2_2= zeros(npoints,n_c);
ave_Phase_lag2_2= zeros(npoints,n_c);
ave_Ar2= zeros(npoints,n_c);
ave_Phase_lag2 = zeros(npoints,n_c);
ave_diss_mean_2 = zeros(npoints,n_c);
ave_virial_mean_2 = zeros(npoints,n_c);
Ar2sin_single = zeros(npoints,(floor(n_c*per/per2)+1));
Ar2cos_single = zeros(npoints,(floor(n_c*per/per2)+1));
ave_A2sin_2 = zeros(npoints,n_c);
ave_A2cos_2 = zeros(npoints,n_c);
ave_sin = zeros(npoints,n_c);
ave_cos = zeros(npoints,n_c);
tan_ave_Phase_lag2_2 = zeros(npoints,n_c);
sin_ave_Phase_lag2_2 = zeros(npoints,n_c);
cos_ave_Phase_lag2_2 = zeros(npoints,n_c);
dissipation1= zeros(1,nstep);
virial1= zeros(1,nstep);
dissipation1m= zeros(npoints,n_c);
virial1m= zeros(npoints,n_c);
dissipation2m= zeros(npoints,n_c);
virial2m= zeros(npoints,n_c);  
dissipation1m_ave= zeros(npoints,n_c);
virial1m_ave= zeros(npoints,n_c);
dissipation2m_ave= zeros(npoints,n_c);
virial2m_ave= zeros(npoints,n_c);        
Ar1m= zeros(npoints,n_c);
Phase_lag1m= zeros(npoints,n_c);
Ar2m= zeros(npoints,n_c);
Phase_lag2m= zeros(npoints,n_c);
time_dec= zeros(n_c,nstep);
qf_dec= zeros(n_c,nstep);
mem= zeros(1,nstep);
h_dec= zeros((floor(n_c*per/per2)+1),nstep);
min_dec= zeros((floor(n_c*per/per2)+1),1);
%initial_t =zeros(n_c,nstep);
%final_t =zeros(n_c,nstep);
start_end= zeros(n_c,2);
ii =zeros(1,n_c);
mm =zeros(1,n_c);
%Ar2m_fit=zeros(npoints,n_c);
%Phase_lag2m_fit=zeros(npoints,n_c);
time_min_q =zeros(n_c,1);
t_min_max =zeros(n_c,2);
time_sim =zeros(n_c,nstep*2);  
time_tot =zeros(n_c,nstep*2);  
n_times=zeros(n_c,3); 

for r=1:npoints

    %% first guess for Z  

    j=1;
    %arstep=1;
    %tip_amp(arstep)=Ar1_d(r)*Af;         %wasn't being used!
    %Zguess(r,1)=tip_amp(arstep)+.1;  %this wasn't even being used!
        
    %sign=1;  %wasn't being used!
    %signchange=1; %wasn't being used!
    dZ=0.1;
    %time_sm(1,:)=linspace(0,per,nstep); !wasnt being used!
    %qf1_norm = cosd(time_sm(1,:)/per*360); !wasnt being used!
    %qdf1_norm = -2*pi/per*sind(time_sm(1,:)/per*360); !wasnt being used!

    if (isempty(Z_single_guess_override))
        funresults=zsingle_VEDA_markedup(casenum,Ar1_d(r),Af,profileype,A, z0, nu, ncct, tau_j, CC_j, CC_0,ODEsolver,nstep,R,maximum_radius,meshnum, tol,k,Q,per);
        z_sm_1(r)=funresults(1);
        Phase_lag_sm_1(r)=funresults(2);
        Ets_sm_1(r)=funresults(3)*6.242e9;
        Vts_sm_1(r)=funresults(4)*6.242e9;
    end
    
    %%
    Zguess(r,:)=zeros(size(Zguess(r,:)));
    Ar1(r,:)=zeros(size(Ar1(r,:)));

    
    %% Main part   

    i=1;
    A2(r,j)=0.1;
    signZ=1;
    signchangeZ=1;
    error=inf;
    condition1=1;

    while condition1

        sinZold=signZ;

        if (i==1)
            signZ=+1;
        elseif (condition1==2)   
            signZ=+1;    
        elseif (Ar1(r,i-1)<Ar1_d(r))
            signZ=+1;
        elseif (Ar1(r,i-1)>Ar1_d(r))   
            signZ=-1;
        end

        if (i>2)&& (signchangeZ==1)&& (condition1~=2)
            signchangeZ=sinZold*signZ;
        end

        if signchangeZ==1
            if (i==1)||(condition1==2)
                dZ=dZ_initial;
            else
                dZ=dZ_initial/2;
            end
        elseif signchangeZ==-1 
            if i>2
                dZ=(Ar1(r,i-1)-Ar1_d(r))*abs(Zguess(r,i-1)-Zguess(r,i-2))/(abs(Ar1(r,i-1)-Ar1_d(r))+abs(Ar1(r,i-2)-Ar1_d(r)));
            else
                dZ=dZ/2;
            end
        end

        if i>1
            Zguess(r,i)=Zguess(r,i-1)+ signZ*dZ;
        elseif i==1
            Zguess(r,i)=z_sm_1(r)+ dZ;
        end

        if (~isempty( Z_single_guess_override))
            Zguess(r,i)=Z_single_guess_override(r);
            condition1=0;
        end
        
        if (i==1)&&(r>1)
            A2_old=A2_final(r-1); 
        elseif (i~=1)
            A2_old=A2_final(r); 
        else
            A2_old=Af2/2;
        end

        j=1;
        A2(r,:)=zeros(size(A2(r,:)));
        A2(r,j)=A2_old;
        sign2=1;
        signchangeA2=1;
        dA2=dA2_initial;

        condition2=1;

        while condition2

            if (j>1)

                signold2=sign2;

                if Ar2(r,j-1)>(A2(r,j-1)/Af2)
                    sign2=+1;
                else
                    sign2=-1;
                end

                if (signchangeA2==1)&&(j>2)
                    signchangeA2= sign2*signold2;
                end

                if signchangeA2==1
                    dA2=dA2_initial;
                elseif (signchangeA2==-1)&&(j>2)
                    dA2=dA2*abs(Ar2(r,j-1)-(A2(r,j-1)/Af2))/(abs(Ar2(r,j-1)-(A2(r,j-1)/Af2))+abs(Ar2(r,j-2)-(A2(r,j-2)/Af2))); 
                else
                    dA2=dA2/2;
                end

                A2(r,j)=A2(r,j-1)+sign2*dA2;

                if A2(r,j)<=0
                    dA2=dA2/2;
                    A2(r,j)=A2(r,j-1)+sign2*dA2;
                end

            end

            if (~isempty( A2_override))
                A2(r,j) = A2_override(j);
                condition2=0;
            end
                
            
        %%      
            time_tot =zeros(n_c*nstep,1);  
            force=zeros(size(time_tot));
            qf1 =zeros(size(time_tot)); 
            qf2 =zeros(size(time_tot)); 
            qdf1 =zeros(size(time_tot)); 
            qdf2 =zeros(size(time_tot)); 
            
            
            %qf1_tot =zeros(size(time_tot)); 
            %qf2_tot =zeros(size(time_tot)); 
            %qdf1_tot =zeros(size(time_tot)); 
            %qdf2_tot =zeros(size(time_tot));                 
            %dissipation1m(r,:) = zeros(size(dissipation1m(r,:)));
            %virial1m(r,:)= zeros(size(virial1m(r,:)));
            %dissipation2m(r,:) = zeros(size(dissipation2m(r,:)));
            %virial2m(r,:)= zeros(size(virial2m(r,:)));
            dissipation1m_ave(r,:) = zeros(size(dissipation1m_ave(r,:)));
            virial1m_ave(r,:)= zeros(size(virial1m_ave(r,:)));
            dissipation2m_ave(r,:) = zeros(size(dissipation2m_ave(r,:)));
            virial2m_ave(r,:)= zeros(size(virial2m_ave(r,:)));
            %Ar1m(r,:)=zeros(size(Ar1m(r,:)));
            %Phase_lag1m(r,:)=zeros(size(Phase_lag1m(r,:)));
            %Ar2m(r,:)=zeros(size(Ar2m(r,:)));
            %Phase_lag2m(r,:)= zeros(size(Phase_lag2m(r,:)));
            %Ar2m_Ext(r,:)=zeros(size(Ar2m_Ext(r,:)));
            %Ar2m_Ext_ave(r,:)=zeros(size(Ar2m_Ext_ave(r,:)));
            %Ar2m_Ext_cyc(r,:)=zeros(size(Ar2m_Ext_cyc(r,:)));
            %Phase_lag2m_Ext(r,:)=zeros(size(Phase_lag2m_Ext(r,:)));
            %Phase_lag2m_Ext_ave(r,:)=zeros(size(Phase_lag2m_Ext_ave(r,:)));
            %Phase_lag2m_Ext_cyc(r,:)=zeros(size(Phase_lag2m_Ext_cyc(r,:)));
            %Ar2m_fit(r,:)=zeros(size(Ar2m_fit(r,:)));
            %Phase_lag2m_fit(r,:)=zeros(size(Phase_lag2m_fit(r,:)));
            

            condition3=1; 
            cyc_numb=1;
            tt=1;
            b=1;bb=1;
            n_Ar2m_Ext_ave_old=0;
            n_Phase_lag2m_Ext_ave_old=0;

            for y=1:n_c
                clear mem
                time_dec(y,:)=linspace((y-1)*per,y*per,nstep);
                qf_dec(y,:)=Zguess(r,i)+(Ar1_d(r)*Af)*cosd(time_dec(y,:)/per*360)+A2(r,j)*cosd((time_dec(y,:)/per2*360));    
                [qw,mem]=find(qf_dec(y,:)<max_height);
                
                %initial_t(y,1:min(mem)-1)=time_dec(y,1:min(mem)-1);
                %final_t(y,1:nstep-max(mem))=time_dec(y,max(mem)+1:nstep);
                start_end(y,1)= time_dec(y,min(mem));
                start_end(y,2)= time_dec(y,max(mem));
            
                ndx = ((y-1)*nstep+1) : (y * nstep);
                time_tot(ndx)=linspace(start_end(y,1), start_end(y,2), nstep);

                qf1(ndx) =(Ar1_d(r)*Af)*cosd(time_tot(ndx)/per*360);
                qf2(ndx) =      A2(r,j)*cosd(time_tot(ndx)/per2*360);
                qdf1(ndx)=(Ar1_d(r)*Af)*(-2*pi/per *sind( time_tot(ndx)/per*360     ));
                qdf2(ndx)=      A2(r,j)*(-2*pi/per2*sind((time_tot(ndx)/per2*360)));
            end


            %bahram only used the 1st cycle as a way to extend the second cycle backwards in time.  didn't use directly
            %so if we want to match his, would need to start the calculation at the second cycle.
            for cyc_numb=1:n_c
                disp(['cycle number', num2str(cyc_numb)]);
                ndx = ((cyc_numb-1)*nstep+1) : (cyc_numb * nstep);
                
                force(ndx) = F_BM_AEM_VEDA_markedup(casenum,cyc_numb,Ar1_d(r),Zguess(r,i)+ qf1(ndx)+qf2(ndx), diff(time_tot(ndx)),profileype, nstep, R, ...
                    maximum_radius, meshnum, A, z0, nu, tau_j, E0, Einf, i , j);                
                
                %this is a kludge.  need zero at ends of each cycle so we don't get interpolation errors in the integration.  do more cleanly fixme
                force(ndx(1))=0;
                force(ndx(end))=0;
             end
                              

            dissipation1m_ave(r) =         -trapz(time_tot,force.*qdf1) / n_c;    % unit N.nm (nJ)                    
            virial1m_ave(r)      = (1/per)* trapz(time_tot,force.* qf1) / n_c;    % unit N.nm (nJ)
            
            Ar1(r,i)=(((1/Q)./sqrt(((1/Q)+(dissipation1m_ave(r,cyc_numb)/pi/k/(Ar1_d(r)*Af)/(Ar1_d(r)*Af))).^2+(2*virial1m_ave(r,cyc_numb)/k/(Ar1_d(r)*Af)/(Ar1_d(r)*Af)).^2)));
            Phase_lag1(r,i)=(atan( ((1/Q)+(dissipation1m_ave(r,cyc_numb)/pi/k/(Ar1_d(r)*Af)/(Ar1_d(r)*Af))) ./(-2*virial1m_ave(r,cyc_numb)/k/(Ar1_d(r)*Af)/(Ar1_d(r)*Af)))*180/pi);
            if Phase_lag1(r,i)<0
              Phase_lag1(r,i)=Phase_lag1m(r)+180;
            end                    
            
            dissipation2m_ave(r) =         -trapz(time_tot,force.*qdf2) / n_c / (per/per2);    % unit N.nm (nJ)                    
            virial2m_ave(r)      = (1/per2)* trapz(time_tot,force.* qf2) / n_c / (per/per2);    % unit N.nm (nJ)
            
            Ar2(r,j)=(1/Q2)/sqrt(((1/Q2)+(dissipation2m_ave(r)/pi/k2/A2(r,j)/A2(r,j))).^2+(-2*virial2m_ave(r)/k2/A2(r,j)/A2(r,j)).^2);
             Phase_lag2(r,j)= atan( ((1/Q2)+(dissipation2m_ave(r)         /pi/k2/A2(r,j)/A2(r,j)))            ./(-2*virial2m_ave(r)         /k2/A2(r,j)/A2(r,j))           )*180/pi;            
            if Phase_lag2(r,j)<0
              Phase_lag2(r,j)=Phase_lag2(r,j)+180;
            end  
            
            if (isempty(A2_override))      
                condition2= (abs(Ar2(r,j)-(A2(r,j)/Af2))>tol_Ar_2)&&(dA2>tol_Ar_2^2);
            end

            %table(converged_cycle, r, i, j, A2(r,j), Ar2(r,j)-(A2(r,j)/Af2), signchangeA2, Zguess(r,i),(Ar1(r,i)-Ar1_d(r)), Ar1_d(r)*Af, signchangeZ)
            'r, i, j, A2(r,j), Ar2(r,j)-(A2(r,j)/Af2), signchangeA2, Zguess(r,i),(Ar1(r,i)-Ar1_d(r)), Ar1_d(r)*Af, signchangeZ'
            [r, i, j, A2(r,j), Ar2(r,j)-(A2(r,j)/Af2), signchangeA2, Zguess(r,i),(Ar1(r,i)-Ar1_d(r)), Ar1_d(r)*Af, signchangeZ]
            
            j=j+1;
        end

        i=i+1;
        errold=error;

        error= abs(Ar1(r,i-1) - Ar1_d(r))+ abs(Ar2(r,j-1)-(A2(r,j-1)/Af2));

        if error<errold
            Z_final(r)=Zguess(r,i-1);                    
            Phase_lag1_final(r)=Phase_lag1(r,i-1);
            A2_final(r)=A2(r,j-1);
            Ar2_final(r)=A2(r,j-1)/Af2;
            Phase_lag2_final(r)=Phase_lag2(r,j-1);
            E_ts_2_final(r)=dissipation2m_ave(r)*6.242e9;
            V_ts_2_final(r)=virial2m_ave(r)*6.242e9;
            E_ts_1_final(r)=dissipation1m_ave(r)*6.242e9;
            V_ts_1_final(r)=virial1m_ave(r)*6.242e9;
            errorfinal_1(r)=Ar1(r,i-1)-Ar1_d(r);
            errorfinal_2(r)=Ar2(r,j-1)-(A2(r,j-1)/Af2);
            it_i(r)=i-1;
            it_j(r)=j-1;
        end

        if (isempty(Z_single_guess_override))
            condition1=(abs(Ar1(r,i-1)-Ar1_d(r))>tol_Ar_1) && (dZ>tol_Ar_1^2);
        end

        %%
        dissipation2m_1=zeros(1,n_c);
        virial2m_1=zeros(1,n_c);
% 
%         if ishandle(3)
%             print('-f3','-dpng','-r1000',[address 'diss.jpg'])
%             print('-f4','-dpng','-r1000',[address 'Virial.jpg'])
%         end

    end %while condition1

      filename1=['case' num2str(ex_opt)  '_th' num2str(r) '.csv'];
      thresults = [time_tot, force, qf1+qf2+Z_final(r), qdf1+qdf2 ];
      csvwrite(filename1,thresults,1)
        
end %for r = 1:npoints

 
results = [1:r; Ar1_d; errorfinal_1+Ar1_d; Z_final;  Phase_lag1_final; A2_final; Phase_lag2_final;  E_ts_2_final;  V_ts_2_final; E_ts_1_final; V_ts_1_final]';
filename1=['case_' num2str(ex_opt) 'simplified_matlab.csv'];
csvwrite(filename1,results,1)

            


% %% plot area
% S_M_P=0;
% 
% figure(2)
% subplot(2,2,1)
% plot(Z_final,Ar1_d,'LineWidth',2)
% if S_M_P==1
%     hold on
%     plot(z_sm_1,Ar1_d,'o--r')
%     legend('Bimodal','Single mode');
% end
% set(get(gca,'XLabel'),'String','Z (nm)','LineWidth',2)
% set(get(gca,'YLabel'),'String','A_1/A_0_1','LineWidth',2);
% 
% subplot(2,2,3)
% plot(Ar1_d,Phase_lag1_final,'LineWidth',2)
% if S_M_P==1
%     hold on
%     plot(Ar1_d,Phase_lag_sm_1,'o--r')
%     legend('Bimodal','Single mode');
% end
% set(get(gca,'XLabel'),'String','A_1/A_0_1','LineWidth',2)
% set(get(gca,'YLabel'),'String','{\phi}_1 (degree)','LineWidth',2);
% subplot(2,2,2)
% plot(Ar1_d, A2_final/Af2,'LineWidth',2)
% set(get(gca,'XLabel'),'String','A_1/A_0_1','LineWidth',2)
% set(get(gca,'YLabel'),'String','A_2/A_0_2','LineWidth',2);
% 
% subplot(2,2,4)
% plot(Ar1_d,Phase_lag2_final,'LineWidth',2)
% set(get(gca,'XLabel'),'String','A_1/A_0_1','LineWidth',2)
% set(get(gca,'YLabel'),'String','{\phi}_2 (degree)','LineWidth',2);

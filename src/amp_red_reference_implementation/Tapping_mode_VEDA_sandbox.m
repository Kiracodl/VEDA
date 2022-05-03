% % Tapping mode investigation by amplitude reduction algorithm %%
% % developed by BAHRAM RAJABIFAR
% % PURDUE UNIVERSITY
% % RAMAN GROUP

%close all
clc
clear
tic


ex_opt = 8;

if (ex_opt ==5)

  %% input parameters
  Sol_type=1;             % 1: Simplified (fourier), 2: Coupled (spatial)
  nstep=1.00E+04;         % Number of time steps	
  %nstep = 300

  R=10;                   % Radius of the tip (nm)	
  maximum_radius=15;   	% computation domain (nm) = 1.5 tip radii

  %meshnum=10
  meshnum=70;         	% number of radial discretization

  k=2.809E-08;        	% Cantilever equivalent stiffness (N/nm)	
                          %28.09 N/m
  Q=429.5;            	% Quality factor
  per=0.0000036;      	% Time Period of the first mode (sec)
                          %277.7778 kHz
  Armax=0.9;          	% max Ar range
  Armin=0.1;              % min Ar range
  npoints=9;        		% number of the Ar selected on the range
  Af=75;              	% free Amp (nm)
  A=1.00E-10;         	% Hamaker's constant (N/nm)	
                          % has to be a typo... must have meant N-nm = 1e-9 N m = 1e-9 J

  z0=0.3;             	% Equilibrium position (nm)
  nu=0.49;            	% Poisson's ratio

  tau_j=1.00E-07;       	% Creep time (Sec)	
                          % eta = T * (E1+E2) = 218.1818                        
                          % props%threeelm_t1 = props%etasample / ( G1 + G2))
                          % G = E / 2 / (1+nu) = 0.33557 E in this case
                          % tau = t1 * E0 / Einf
                          
                          % eta = 0.33557 (E1 + E2) * t1 * Einf / E0 
                          
                          
                          
  E0=2e-9;               	% N/(nm)^2 = 1e-18 GPa 
                          % typo? Meant 1e-18 Pa? not GPa?
                          % 2 GPa = E1
  Einf=0.2e-9;           	% N/(nm)^2 = 1e-18 GPa
                          % 0.2 GPa
                          %E2 = (Einf E0 ) / (E0 - Einf) = (0.2 * 2)/ (2-0.2)=0.222222 GPa
  tol=0.01;           	% predefined tolerance for Ar predictions
  basis_func_number=70; 	% number of basis functions, only needed for simplified model
elseif (ex_opt == 6)
  %% input parameters
  Sol_type=1;             % 1: Simplified (fourier), 2: Coupled (spatial)
  nstep=1.00E+04;         % Number of time steps	
  %nstep = 300

  R=10;                   % Radius of the tip (nm)	
  maximum_radius=2*R;   	% computation domain (nm) = 1.5 tip radii

  %meshnum=10
  meshnum=100;         	% number of radial discretization

  k=5E-09;        	% Cantilever equivalent stiffness (N/nm)	
     
  Q=300;            	% Quality factor
  per=1/100e3;      	% Time Period of the first mode (sec)             
  Armax=0.9;          	% max Ar range
  Armin=0.1;              % min Ar range
  npoints=9;        		% number of the Ar selected on the range
  Af=20;              	% free Amp (nm)
  A=1.00E-10;         	% Hamaker's constant (N/nm)	
                          % has to be a typo... must have meant N-nm 

  z0=0.3;             	% Equilibrium position (nm)
  nu=0.49;            	% Poisson's ratio

  tau_j=2.2219774000000000E-007;       	% Creep time (Sec)	                         
                          
  E0starPa = 1315962626.6614029 % Pa
  E0_Pa = E0starPa*(1-nu*nu);  
  E0= E0_Pa * 1e-18;               	% N/(nm)^2 = 1e-18 GPa . 
                                    %typo? meant 1e-18 Pa? 
          
  Einf_star_Pa=119632966.06012751;
  Einf_Pa = Einf_star_Pa*(1-nu^2);
  Einf    = Einf_Pa * 1e-18;  % N/(nm)^2 = 1e-18 GPa
                              %typo? meant 1e-18 Pa? %
                          
  tol=0.01;           	% predefined tolerance for Ar predictions
  basis_func_number=20; 	% number of basis functions, only needed for simplified model
elseif (ex_opt == 6.1)
  %stripped down version of 6
  
  Sol_type=2;             % 1: Simplified (fourier), 2: Coupled (spatial)
  nstep=150;         % Number of time steps	
  R=10;                   % Radius of the tip (nm)	
  maximum_radius=2*R;   	% computation domain (nm) = 1.5 tip radii
  meshnum=20;         	% number of radial discretization
  basis_func_number=15; 	% number of basis functions, only needed for simplified model
  k=5E-09;        	% Cantilever equivalent stiffness (N/nm)	
  Q=300;            	% Quality factor
  per=1/100e3;      	% Time Period of the first mode (sec)             
  Armax=0.9;          	% max Ar range
  Armin=0.1;              % min Ar range
  npoints=9;        		% number of the Ar selected on the range
  Af=20;              	% free Amp (nm)
  A=1.00E-10;         	% Hamaker's constant (N/nm)	
  z0=0.3;             	% Equilibrium position (nm)
  nu=0.49;            	% Poisson's ratio
  tau_j=2.2219774000000000E-007;       	% Creep time (Sec)	                                                  
  E0starPa = 1315962626.6614029 % Pa
  E0_Pa = E0starPa*(1-nu*nu);  
  E0= E0_Pa * 1e-18;               	% N/(nm)^2 = 1e-18 GPa . 
  Einf_star_Pa=119632966.06012751;
  Einf_Pa = Einf_star_Pa*(1-nu^2);
  Einf    = Einf_Pa * 1e-18;  % N/(nm)^2 = 1e-18 GPa
  tol=0.01;           	% predefined tolerance for Ar predictions
elseif (ex_opt == 7)
  %fixme, rerun, hamaker was wrong
  
  %% input parameters
  Sol_type=1;             % 1: Simplified (fourier), 2: Coupled (spatial)
  
  nstep=1.00E+04;         % Number of time steps	
  
  R=10;                   % Radius of the tip (nm)	
  maximum_radius=3*R;   	% computation domain (nm) = 1.5 tip radii

  %meshnum=10
  meshnum=50;         	% number of radial discretization

  k=10E-09;        	% Cantilever equivalent stiffness (N/nm)	
     
  Q=300;            	% Quality factor
  per=1/100e3;      	% Time Period of the first mode (sec)             
  Armax=0.99;          	% max Ar range
  Armin=0.09;              % min Ar range
  npoints=9;        		% number of the Ar selected on the range
  Af=20;              	% free Amp (nm)
  A=1.00E-11;         	% Hamaker's constant (N/nm)	
                          % has to be a typo... must have meant N-nm 

  z0=0.3;             	% Equilibrium position (nm)
  nu=0.49;            	% Poisson's ratio

  tau_j=2.2219774E-006;       	% Creep time (Sec)	                         
                          
  E0starPa = 1315962626.6614029 % Pa
  E0_Pa = E0starPa*(1-nu*nu);  
  E0= E0_Pa * 1e-18;               	% N/(nm)^2 = 1e-18 GPa . 
                                    %typo? meant 1e-18 Pa? 
          
  Einf_star_Pa=119632966.06012751;
  Einf_Pa = Einf_star_Pa*(1-nu^2);
  Einf    = Einf_Pa * 1e-18;  % N/(nm)^2 = 1e-18 GPa
                              %typo? meant 1e-18 Pa? %
                          
  tol=0.001;           	% predefined tolerance for Ar predictions
  basis_func_number=20; 	% number of basis functions, only needed for simplified model
elseif (ex_opt == 7.5)
  %fixme, rerun, hamaker was wrong
  
  %same as 7, but spatial instead of fourier.  
  %fails to run.  ill conditioned matrices
  %% input parameters
  Sol_type=2;             % 1: Simplified (fourier), 2: Coupled (spatial)
  
  nstep=1.00E+04;         % Number of time steps	
  
  R=10;                   % Radius of the tip (nm)	
  maximum_radius=3*R;   	% computation domain (nm) = 1.5 tip radii

  %meshnum=10
  meshnum=50;         	% number of radial discretization

  k=10E-09;        	% Cantilever equivalent stiffness (N/nm)	
     
  Q=300;            	% Quality factor
  per=1/100e3;      	% Time Period of the first mode (sec)             
  Armax=0.99;          	% max Ar range
  Armin=0.09;              % min Ar range
  npoints=9;        		% number of the Ar selected on the range
  Af=20;              	% free Amp (nm)
  A=1.00E-11;         	% Hamaker's constant (N/nm)	
                          % has to be a typo... must have meant N-nm 

  z0=0.3;             	% Equilibrium position (nm)
  nu=0.49;            	% Poisson's ratio

  tau_j=2.2219774E-006;       	% Creep time (Sec)	                         
                          
  E0starPa = 1315962626.6614029 % Pa
  E0_Pa = E0starPa*(1-nu*nu);  
  E0= E0_Pa * 1e-18;               	% N/(nm)^2 = 1e-18 GPa . 
                                    %typo? meant 1e-18 Pa? 
          
  Einf_star_Pa=119632966.06012751;
  Einf_Pa = Einf_star_Pa*(1-nu^2);
  Einf    = Einf_Pa * 1e-18;  % N/(nm)^2 = 1e-18 GPa
                              %typo? meant 1e-18 Pa? %
                          
  tol=0.001;           	% predefined tolerance for Ar predictions
  %basis_func_number=20; 	% number of basis functions, only needed for simplified model
elseif (ex_opt == 8)
  %same material properties as 7, but larger amplitude and stiffer k.  
  Sol_type=1;             % 1: Simplified (fourier), 2: Coupled (spatial)
  
  nstep=1.00E+04;         % Number of time steps	
  
  R=10;                   % Radius of the tip (nm)	
  maximum_radius=3*R;   	% computation domain (nm) = 1.5 tip radii

  %meshnum=10
  meshnum=50;         	% number of radial discretization

  k=40E-09;        	% Cantilever equivalent stiffness (N/nm)	
     
  Q=300;            	% Quality factor
  per=1/100e3;      	% Time Period of the first mode (sec)             
  Armax=0.99;          	% max Ar range
  Armin=0.09;              % min Ar range
  npoints=9;        		% number of the Ar selected on the range
  Af=50;              	% free Amp (nm)
  A=1.00E-11;         	% Hamaker's constant (N/nm)	
                          % has to be a typo... must have meant N-nm 

  z0=0.3;             	% Equilibrium position (nm)
  nu=0.49;            	% Poisson's ratio

  tau_j=2.2219774E-006;       	% Creep time (Sec)	                         
                          
  E0starPa = 1315962626.6614029 % Pa
  E0_Pa = E0starPa*(1-nu*nu);  
  E0= E0_Pa * 1e-18;               	% N/(nm)^2 = 1e-18 GPa . 
                                    %typo? meant 1e-18 Pa? 
          
  Einf_star_Pa=119632966.06012751;
  Einf_Pa = Einf_star_Pa*(1-nu^2);
  Einf    = Einf_Pa * 1e-18;  % N/(nm)^2 = 1e-18 GPa
                              %typo? meant 1e-18 Pa? %
                          
  tol=0.001;           	% predefined tolerance for Ar predictions
  basis_func_number=20; 	% number of basis functions, only needed for simplified model

  end
  
save_all_figures = 0

%% Initialization
Ar=linspace(Armax,Armin,npoints);
amp=Ar*Af;
sign=1;
dz_0=2.5;                   %nm 
Zguess=zeros(size(amp));
Zguess(1)=amp(1)+dz_0;
signchange=0;
results=zeros(length(Ar),8);
address=strcat(pwd,'\'); 
casenum=1;
arstep=1;
jjj=1; 
dz=dz_0;

%% Main computation part 



while (arstep<(length(Ar)+1)) 
    
    %    %take the iteration out of it for debugging
    %this is for case 5.
    %Zguess = [63.125,54.53125,46.640625,38.896484375,31.243896484375,23.743896484375,16.556396484375,9.447021484375,2.435302734375];
    %this is case 6    
    %Zguess = [  18.898438,    16.825195,    14.764160,    12.687866,    10.687866,     8.747471,     6.821977,     4.915109,     3.002420];
    
    if Sol_type==1   
        [final, time, Force, gap_to_undeform, dpr, u0] =u_Fourier_VEDA_markedup(basis_func_number,meshnum,nstep,Zguess(arstep)*1e-9,amp(arstep)*1e-9,R*1e-9,maximum_radius*1e-9,1/per,A*1e-9,z0*1e-9,nu,Einf*1e18,E0*1e18,tau_j);        
    elseif Sol_type==2
        [final, time, Force, gap_to_undeform, dpr, u0] =coupled_VEDA_sandbox(meshnum,nstep,Zguess(arstep)*1e-9,amp(arstep)*1e-9,R*1e-9,maximum_radius*1e-9,1/per,A*1e-9,z0*1e-9,nu,Einf*1e18,E0*1e18,tau_j);
    end
%%    

    dissipation = final(1)*1e9;    % unit N.nm (nJ)
    virial= final(2)*1e9;         % unit N.nm (nJ)

     Ar_guess=1/Q/sqrt(((1/Q)+(dissipation/pi/k/amp(arstep)/amp(arstep)))^2+...
         (2*virial/k/amp(arstep)/amp(arstep))^2);

     disp(['Zguess ', num2str(Zguess(arstep)), ' error ',  num2str(Ar_guess-Ar(arstep)) ])

    signold=sign;
    dzold=dz;

%    if (0)
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
%    elseif (1)
    elseif abs(Ar_guess-Ar(arstep))<tol
        if arstep<length(Ar)
            Zguess(arstep+1)=amp(arstep+1)+Zguess(arstep)-amp(arstep);
            dz=min(dz_0,5*dz);
        end

        %%
        
%         columnheaders={'time', 'force', 'gap_to_undeform', 'dpr', 'surf_coord'};
%         thsave = array2table([time(1:end-1); Force; gap_to_undeform(1:end-1); dpr(1:end-1); u0(1:end-1)]', 'VariableNames', columnheaders);
%         if (Sol_type == 1)
%             tag = 'fourier';
%         else
%             tag = 'spatial';
%         end
%         
%         filename1=['case5_' tag '_th' num2str(arstep) '_meshnum' num2str(meshnum) '_nstep' num2str(nstep) '.csv'];
%             
%         writetable(thsave,filename1)
        filename1=['case' num2str(ex_opt)  '_th' num2str(arstep) '_meshnum' num2str(meshnum) '_nstep' num2str(nstep) '_octave.csv'];
        thresults = [time(1:end-1); Force; gap_to_undeform(1:end-1); dpr(1:end-1); u0(1:end-1)]';
        csvwrite(filename1,thresults,1)
        
        if (save_all_figures)
            fig = gcf;
            filename4=[address 'results\' num2str(casenum) '_' num2str(arstep) '_tapping.jpg'];
            saveas(fig,filename4)
        end

        phase=atan(((1/Q)+(dissipation/pi/k/amp(arstep)/amp(arstep)))/(-2*virial/k/amp(arstep)/amp(arstep)))*180/pi;
        if phase<0
            phase=phase+180;
        end

        results(jjj,:)=[arstep Ar(arstep) Ar_guess Zguess(arstep) amp(arstep) dissipation*6.242e9 virial*6.242e9 phase]; 
        jjj=jjj+1;
        arstep=arstep+1;    
    end
    
%     phase0=atan(((1/Q)+(dissipation/pi/k/amp(arstep)/amp(arstep)))/(-2*virial/k/amp(arstep)/amp(arstep)))*180/pi;
%     if phase0<0
%         phase0=phase0+180;
%     end
%     
%     [ Ar(arstep) Ar_guess Zguess(arstep) phase0 dz]

end

toc

%% results output generation

figure(2)
subplot(2,2,1)
plot(results(:,2),results(:,6),'.','MarkerSize',25)
set(get(gca,'XLabel'),'String','Amplitude ratio', 'FontSize', 12);
set(get(gca,'YLabel'),'String','Dissipation (ev/cycle)', 'FontSize', 12);
subplot(2,2,2)
plot(results(:,2),-results(:,7),'.','MarkerSize',25)
set(get(gca,'XLabel'),'String','Amplitude ratio', 'FontSize', 12);
set(get(gca,'YLabel'),'String','Virial (ev/cycle)', 'FontSize', 12);
subplot(2,2,3:4)
plot(results(:,2),results(:,8),'.','MarkerSize',25)
set(get(gca,'XLabel'),'String','Amplitude ratio', 'FontSize', 12);
set(get(gca,'YLabel'),'String','Phase (deg)', 'FontSize', 12);


%% Saving, Bahram's versions
% 
% resultssave(1,:)={'arstep' 'Ar(arstep)' 'Ar_guess' 'Zguess(arstep)'  'amp(arstep)' 'dissipation*6.242e9' 'virial*6.242e9' 'Phase'};
% for i=1:length(results(:,1))
%     resultssave(i+1,:)=num2cell(results(i,:));
% end
% 
% filename1=[address 'results\' num2str(casenum) '_tapping.xlsx'];
% filename2=[address 'results\' num2str(casenum) '_tapping.jpg'];
% 
% print('-f2',filename2,'-djpeg','-r500')
% xlswrite(filename1,resultssave,1)

%%
% 
filename1=['case_' num2str(ex_opt) '_octave.csv'];
csvwrite(filename1,results,1)

%columnheaders={'arstep' 'Ar(arstep)' 'Ar_guess' 'Zguess(arstep)'  'amp(arstep)' 'dissipation' 'virial' 'Phase'};
%resultssave = array2table(results, 'VariableNames', columnheaders)
%filename1=['case5_' tag '.csv'];
%writetable(resultssave,filename1)

% % Tapping mode investigation by amplitude reduction algorithm %%
% % developed by BAHRAM RAJABIFAR
% % PURDUE UNIVERSITY
% % RAMAN GROUP

close all
clc
clear
tic

%% input parameters
Sol_type=1;             % 1: Simplified, 2: Coupled
nstep=1.00E+04;         % Number of time steps	
R=10;                   % Radius of the tip (nm)	
maximum_radius=15;   	% computation domain (nm)
meshnum=70;         	% number of radial discretization
k=2.809E-08;        	% Cantilever equivalent stiffness (N/nm)	
Q=429.5;            	% Quality factor
per=0.0000036;      	% Time Period of the first mode (sec)
Armax=0.9;          	% max Ar range
Armin=0.1;              % min Ar range
npoints=9;        		% number of the Ar selected on the range
Af=75;              	% free Amp (nm)
A=1.00E-10;         	% Hamaker's constant (N/nm)	
z0=0.3;             	% Equilibrium position (nm)
nu=0.49;            	% Poisson's ratio
tau_j=1.00E-07;       	% Creep time (Sec)	
E0=2e-9;               	% N/(nm)^2 = 1e-18 GPa
Einf=0.2e-9;           	% N/(nm)^2 = 1e-18 GPa
tol=0.01;           	% predefined tolerance for Ar predictions
basis_func_number=25; 	% number of basis functions, only needed for simplified model


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

    if Sol_type==1   
        final=u_Fourier_VEDA(basis_func_number,meshnum,nstep,Zguess(arstep)*1e-9,amp(arstep)*1e-9,R*1e-9,maximum_radius*1e-9,1/per,A*1e-9,z0*1e-9,nu,Einf*1e18,E0*1e18,tau_j);
    elseif Sol_type==2
        final=coupled_VEDA(meshnum,nstep,Zguess(arstep)*1e-9,amp(arstep)*1e-9,R*1e-9,maximum_radius*1e-9,1/per,A*1e-9,z0*1e-9,nu,Einf*1e18,E0*1e18,tau_j);
    end
%%    

    dissipation = final(1)*1e9;    % unit N.nm (nJ)
    virial= final(2)*1e9;         % unit N.nm (nJ)

    Ar_guess=1/Q/sqrt(((1/Q)+(dissipation/pi/k/amp(arstep)/amp(arstep)))^2+...
        (2*virial/k/amp(arstep)/amp(arstep))^2);

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
            dz=min(dz_0,5*dz);
        end

        %%
        fig = gcf;
        filename4=[address 'results\' num2str(casenum) '_' num2str(arstep) '_tapping.jpg'];
        saveas(fig,filename4)

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


%% Saving

resultssave(1,:)={'arstep' 'Ar(arstep)' 'Ar_guess' 'Zguess(arstep)'  'amp(arstep)' 'dissipation*6.242e9' 'virial*6.242e9' 'Phase'};
for i=1:length(results(:,1))
    resultssave(i+1,:)=num2cell(results(i,:));
end

filename1=[address 'results\' num2str(casenum) '_tapping.xlsx'];
filename2=[address 'results\' num2str(casenum) '_tapping.jpg'];

print('-f2',filename2,'-djpeg','-r500')
xlswrite(filename1,resultssave,1)


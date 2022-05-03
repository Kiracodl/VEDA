opengl('save', 'software');
   
clear
clc
close all

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
E0=2e-9;               	% N/(nm)^2 = 1e-18 GPa
Einf=0.2e-9;           	% N/(nm)^2 = 1e-18 GPa
tol=0.01;           	% predefined tolerance for Ar predictions
basis_func_number=25; 	% number of basis functions, only needed for simplified model
n_c=11;                 % max number of the 1st mode oscilations to study the convergence

casenum=1;
profileype=1;
ncct=1;
pl=0;
ODEsolver=1;
address = strcat(pwd,'\results\');  

dA2_initial=Af2/10;
dZ_initial=Af2/30;

tol_Ar_1=0.025;
tol_Ar_2=0.025;
AR2_tol=tol_Ar_2;
PL2_tol=0.1;
conv_lim=20; % number of cycles with less than Ar2 and PL2 tolerance deviation
Ar1_d=fliplr(linspace(Armin,Armax,npoints));  
beta=5^(2/nstep); %mesh refinement factor
mesh_refine=0;  % 1:Yes, 0:No
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
initial_t =zeros(n_c,nstep);
final_t =zeros(n_c,nstep);
start_end= zeros(n_c,2);
ii =zeros(1,n_c);
mm =zeros(1,n_c);
Ar2m_fit=zeros(npoints,n_c);
Phase_lag2m_fit=zeros(npoints,n_c);
time_min_q =zeros(n_c,1);
t_min_max =zeros(n_c,2);
time_sim =zeros(n_c,nstep*2);  
time_tot =zeros(n_c,nstep*2);  
n_times=zeros(n_c,3); 

for r=1:npoints

    %% first guess for Z  

    i=1;
    j=1;
    arstep=1;
    tip_amp(arstep)=Ar1_d(r)*Af;
    Zguess(r,i)=tip_amp(arstep)+.1;
    sign=1;
    signchange=1;
    dZ=0.1;
    time_sm(1,:)=linspace(0,per,nstep);
    qf1_norm = cosd(time_sm(1,:)/per*360);
    qdf1_norm = -2*pi/per*sind(time_sm(1,:)/per*360);

    funresults=zsingle_VEDA(casenum,Ar1_d(r),Af,profileype,A, z0, nu, ncct, tau_j, CC_j, CC_0...
         ,ODEsolver,nstep,R,maximum_radius,meshnum, tol,k,Q,per);
    z_sm_1(r)=funresults(1);
    Phase_lag_sm_1(r)=funresults(2);
    Ets_sm_1(r)=funresults(3)*6.242e9;
    Vts_sm_1(r)=funresults(4)*6.242e9;

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

        %%      
            qf1 =zeros(size(time_tot)); 
            qf2 =zeros(size(time_tot)); 
            qdf1 =zeros(size(time_tot)); 
            qdf2 =zeros(size(time_tot)); 
            time_tot =zeros(size(time_tot)); 
            time_tot_per2=zeros(size(time_tot)); 
            force_tot=zeros(size(time_tot)); 
            qf1_tot =zeros(size(time_tot)); 
            qf2_tot =zeros(size(time_tot)); 
            qdf1_tot =zeros(size(time_tot)); 
            qdf2_tot =zeros(size(time_tot));                 
            dissipation1m(r,:) = zeros(size(dissipation1m(r,:)));
            virial1m(r,:)= zeros(size(virial1m(r,:)));
            dissipation2m(r,:) = zeros(size(dissipation2m(r,:)));
            virial2m(r,:)= zeros(size(virial2m(r,:)));
            dissipation1m_ave(r,:) = zeros(size(dissipation1m_ave(r,:)));
            virial1m_ave(r,:)= zeros(size(virial1m_ave(r,:)));
            dissipation2m_ave(r,:) = zeros(size(dissipation2m_ave(r,:)));
            virial2m_ave(r,:)= zeros(size(virial2m_ave(r,:)));
            Ar1m(r,:)=zeros(size(Ar1m(r,:)));
            Phase_lag1m(r,:)=zeros(size(Phase_lag1m(r,:)));
            Ar2m(r,:)=zeros(size(Ar2m(r,:)));
            Phase_lag2m(r,:)= zeros(size(Phase_lag2m(r,:)));
            Ar2m_Ext(r,:)=zeros(size(Ar2m_Ext(r,:)));
            Ar2m_Ext_ave(r,:)=zeros(size(Ar2m_Ext_ave(r,:)));
            Ar2m_Ext_cyc(r,:)=zeros(size(Ar2m_Ext_cyc(r,:)));
            Phase_lag2m_Ext(r,:)=zeros(size(Phase_lag2m_Ext(r,:)));
            Phase_lag2m_Ext_ave(r,:)=zeros(size(Phase_lag2m_Ext_ave(r,:)));
            Phase_lag2m_Ext_cyc(r,:)=zeros(size(Phase_lag2m_Ext_cyc(r,:)));
            Ar2m_fit(r,:)=zeros(size(Ar2m_fit(r,:)));
            Phase_lag2m_fit(r,:)=zeros(size(Phase_lag2m_fit(r,:)));
            time_3 =zeros(n_c,nstep);
            qf_dec_2 =zeros(n_c,nstep);
            force=zeros(size(time_sim));
            t_shared=zeros(n_c,4);
            tp_shared=zeros(n_c,2);
            m_shared=zeros(n_c,4);
            n_m_shared=zeros(n_c,2);
            force_tot_exp=zeros(n_c,2*nstep);
            coef_tot_exp=zeros(n_c,2*nstep);
            time_tot_exp=zeros(n_c,2*nstep);
            qf1_tot_exp=zeros(n_c,2*nstep);
            qf2_tot_exp=zeros(n_c,2*nstep);
            qdf1_tot_exp=zeros(n_c,2*nstep);
            qdf2_tot_exp=zeros(n_c,2*nstep);

            condition3=1; 
            cyc_numb=1;
            tt=1;
            b=1;bb=1;
            n_Ar2m_Ext_ave_old=0;
            n_Phase_lag2m_Ext_ave_old=0;

            for y=1:n_c
                clear mem
                time_dec(y,:)=linspace((y-1)*per,y*per,nstep);
                qf_dec(y,:)=Zguess(r,i)+(Ar1_d(r)*Af)*cosd(time_dec(y,:)/per*360)+...
                    A2(r,j)*cosd((time_dec(y,:)/per2*360)-pl);    
                [qw,mem]=find(qf_dec(y,:)<max_height);
                initial_t(y,1:min(mem)-1)=time_dec(y,1:min(mem)-1);
                final_t(y,1:nstep-max(mem))=time_dec(y,max(mem)+1:nstep);
                start_end(y,1)= time_dec(y,min(mem));
                start_end(y,2)= time_dec(y,max(mem));

                if mesh_refine==1
                    [mm,ii(y)]=min(transpose(qf_dec(y,:)));
                    time_min_q(y,1)=time_dec(y,ii(y));
                    t_min_max(y,:)=sort(abs(start_end(y,:)-time_min_q(y)));
                    fun=@(k) ((beta.^(k+1)-1)/(beta.^(nstep-k+1)-1))-(t_min_max(y,2)/t_min_max(y,1));
                    x = fzero(fun,t_min_max(y,2)/(t_min_max(y,2)+ t_min_max(y,1))*nstep);
                    delta_t=t_min_max(y,2)*(beta-1)/((beta^(x+1))-1);
                    time_sim(y,floor(x))=time_min_q(y,1);

                    permit=1;
                    for g=1:floor(x)
                        if permit==1
                            time_sim(y,floor(x)+g)=time_sim(y,floor(x)+g-1)+ delta_t*(beta^g);
                            if time_sim(y,floor(x)+g-1)>start_end(y,2)
                                permit=0;
                            end
                        end
                        if floor(x)-g>0
                            time_sim(y,floor(x)-g)=time_sim(y,floor(x)-g+1)-delta_t*(beta^g);
                        end
                    end
                elseif mesh_refine==0
                    time_sim(y,1:nstep)=linspace(start_end(y,1), start_end(y,2), nstep);
                end                  

                n_times(y,:)=[length(initial_t(y,1:min(mem)-1)) length(initial_t(y,1:min(mem)-1))+length(time_sim(y,1:nnz(time_sim(y,:)))) ...
                    length(initial_t(y,1:min(mem)-1))+length(time_sim(y,1:nnz(time_sim(y,:))))+length(final_t(y,1:nstep-max(mem)))];
                time_tot(y,1:n_times(y,3))= [initial_t(y,1:min(mem)-1) time_sim(y,1:nnz(time_sim(y,:))) final_t(y,1:nstep-max(mem))];

                qf1(y,1:n_times(y,2)-n_times(y,1))=(Ar1_d(r)*Af)*cosd(time_tot(y,n_times(y,1)+1:n_times(y,2))/per*360);
                qf2(y,1:n_times(y,2)-n_times(y,1))=A2(r,j)*cosd((time_tot(y,n_times(y,1)+1:n_times(y,2))/per2*360)-pl);
                qdf1(y,1:n_times(y,2)-n_times(y,1))=(Ar1_d(r)*Af)*(-2*pi/per*sind(time_tot(y,n_times(y,1)+1:n_times(y,2))/per*360));
                qdf2(y,1:n_times(y,2)-n_times(y,1))= A2(r,j)*(-2*pi/per2*sind((time_tot(y,n_times(y,1)+1:n_times(y,2))/per2*360)-pl));

                qf1_tot(y,1:n_times(y,3))=(Ar1_d(r)*Af)*cosd(time_tot(y,1:n_times(y,3))/per*360);
                qf2_tot(y,1:n_times(y,3))=A2(r,j)*cosd((time_tot(y,1:n_times(y,3))/per2*360)-pl);
                qdf1_tot(y,1:n_times(y,3))=(Ar1_d(r)*Af)*(-2*pi/per*sind(time_tot(y,1:n_times(y,3))/per*360));
                qdf2_tot(y,1:n_times(y,3))= A2(r,j)*(-2*pi/per2*sind((time_tot(y,1:n_times(y,3))/per2*360)-pl));

                time_tot_per2(y,1:n_times(y,3)) = floor(time_tot(y,1:n_times(y,3))/per2)+1;

                if y>2
                    clear a1 a2 a3 a4 a5 a6 a7 a8

                    [a1,a2]=find(time_tot_per2(y-1,1:n_times(y-1,3)) == time_tot_per2(y-2,n_times(y-2,3)));
                    [a3,a4]=find(time_tot_per2(y-2,1:n_times(y-2,3)) == time_tot_per2(y-1,1));
                    [a5,a6]=find(time_tot_per2(y-1,1:n_times(y-1,3)) == time_tot_per2(y,1));
                    [a7,a8]=find(time_tot_per2(y,1:n_times(y,3)) == time_tot_per2(y-1,n_times(y-1,3)));

                    if y==3
                        m_shared(1,:)= [1 1 min(a4) max(a4)];
                        n_m_shared(1,2)= m_shared(1,4)-m_shared(1,3)+1;
                    end

                    m_shared(y-1,:)=[min(a2) max(a2) min(a6) max(a6)];
                    t_shared(y-1,:)=[(time_tot(y-2,max(a4))-time_tot(y-2,min(a4))) (time_tot(y-1,max(a2))-time_tot(y-1,min(a2))) ...
                        (time_tot(y-1,max(a6))-time_tot(y-1,min(a6))) (time_tot(y,max(a8))-time_tot(y,min(a8)))];
                    tp_shared(y-1,:)=[t_shared(y-1,2)/(t_shared(y-1,1)+t_shared(y-1,2)) t_shared(y-1,3)/(t_shared(y-1,3)+t_shared(y-1,4))];
                    n_m_shared(y-1,:)= [m_shared(y-1,2)-m_shared(y-1,1)+1  m_shared(y-1,4)-m_shared(y-1,3)+1];
                    
                end

            end

            tp_shared(1,:)=[0 1-tp_shared(2,1)];

            while condition3==1 

                if time_sim(cyc_numb,1)==0
                    n_i=nnz(time_sim(cyc_numb,:))+1;
                else
                    n_i=nnz(time_sim(cyc_numb,:));
                end

                force(cyc_numb,1:n_i) = F_BM_AEM_VEDA(casenum,cyc_numb,Ar1_d(r),Zguess(r,i)+ qf1(cyc_numb,1:n_i)...
                    + qf2(cyc_numb,1:n_i), diff(time_sim(cyc_numb,1:n_i)),profileype, nstep, R, ...
                    maximum_radius, meshnum, A, z0, nu, tau_j, E0, Einf, i , j);                

                force_tot(cyc_numb,1:n_times(cyc_numb,3))=...
                    [zeros(1,n_times(cyc_numb,1)) force(cyc_numb,1:(n_times(cyc_numb,2)-n_times(cyc_numb,1))) zeros(1,(n_times(cyc_numb,3)-n_times(cyc_numb,2)))];

                if cyc_numb>1
                    cyc_numb=cyc_numb-1;

                    if cyc_numb==1
                        force_tot_exp(cyc_numb,1:n_times(cyc_numb,3)+n_times(cyc_numb+1,1))=[force_tot(cyc_numb,1:n_times(cyc_numb,3)) force_tot(cyc_numb+1,1:n_times(cyc_numb+1,1))];
                        coef_tot_exp(cyc_numb,1:nnz(time_tot(cyc_numb,:))+n_m_shared(cyc_numb+1,1))=[ones(1,nnz(time_tot(cyc_numb,:))) tp_shared(cyc_numb,2)*ones(1,n_m_shared(cyc_numb+1,1))];
                        time_tot_exp(cyc_numb,1:nnz(time_tot(cyc_numb,:))+n_m_shared(cyc_numb+1,1)+1)=[time_tot(cyc_numb,m_shared(cyc_numb,1):m_shared(cyc_numb,4)) time_tot(cyc_numb+1,m_shared(cyc_numb+1,1):m_shared(cyc_numb+1,2))];
                        qf1_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:))+1)=(Ar1_d(r)*Af)*cosd(time_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:))+1)/per*360);
                        qf2_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:))+1)=A2(r,j)*cosd(time_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:))+1)/per2*360-pl);
                        qdf1_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:))+1)=(Ar1_d(r)*Af)*(-2*pi/per*sind(time_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:))+1)/per*360));
                        qdf2_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:))+1)= A2(r,j)*(-2*pi/per2*sind(time_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:))+1)/per2*360-pl));                            
                    else
                        force_tot_exp(cyc_numb,1:m_shared(cyc_numb-1,4)-m_shared(cyc_numb-1,3)+m_shared(cyc_numb,4)+m_shared(cyc_numb+1,2)+1)=...
                            [force_tot( cyc_numb-1, m_shared(cyc_numb-1,3):m_shared(cyc_numb-1,4)) force_tot(cyc_numb,1:m_shared(cyc_numb,4)) force_tot(cyc_numb+1,1:m_shared(cyc_numb+1,2))];

                        coef_tot_exp(cyc_numb,1:m_shared(cyc_numb-1,4)-m_shared(cyc_numb-1,3)+m_shared(cyc_numb,4)+m_shared(cyc_numb+1,2)+1)=...
                            [tp_shared(cyc_numb,1)*ones(1,m_shared(cyc_numb-1,4)-m_shared(cyc_numb-1,3)+m_shared(cyc_numb,2)+1) ones(1,m_shared(cyc_numb,3)-m_shared(cyc_numb,2)) tp_shared(cyc_numb,2)*ones(1,m_shared(cyc_numb,4)-m_shared(cyc_numb,3)+m_shared(cyc_numb+1,2))];

                        time_tot_exp(cyc_numb,1:m_shared(cyc_numb-1,4)-m_shared(cyc_numb-1,3)+m_shared(cyc_numb,4)+m_shared(cyc_numb+1,2)+1)=...
                            [time_tot( cyc_numb-1, m_shared(cyc_numb-1,3):m_shared(cyc_numb-1,4)) time_tot(cyc_numb,1:m_shared(cyc_numb,4)) time_tot(cyc_numb+1,1:m_shared(cyc_numb+1,2))];

                        qf1_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:)))=(Ar1_d(r)*Af)*cosd(time_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:)))/per*360);
                        qf2_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:)))=A2(r,j)*cosd(time_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:)))/per2*360-pl);
                        qdf1_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:)))=(Ar1_d(r)*Af)*(-2*pi/per*sind(time_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:)))/per*360));
                        qdf2_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:)))= A2(r,j)*(-2*pi/per2*sind(time_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:)))/per2*360-pl));
                    end

                    dissipation1m(r,cyc_numb) = -trapz(time_tot(cyc_numb,1:n_times(cyc_numb,3)),force_tot(cyc_numb,1:n_times(cyc_numb,3)).*qdf1_tot(cyc_numb,1:n_times(cyc_numb,3)));    % unit N.nm (nJ)
                    virial1m(r,cyc_numb)= (1/per)* trapz(time_tot(cyc_numb,1:n_times(cyc_numb,3)),force_tot(cyc_numb,1:n_times(cyc_numb,3)).*qf1_tot(cyc_numb,1:n_times(cyc_numb,3)));    % unit N.nm (nJ)
                    dissipation1m_ave(r,cyc_numb) = mean(dissipation1m(r,1:cyc_numb));    % unit N.nm (nJ)
                    virial1m_ave(r,cyc_numb)= mean(virial1m(r,1:cyc_numb));    % unit N.nm (nJ)
                    Ar1m(r,cyc_numb)=(((1/Q)./sqrt(((1/Q)+(dissipation1m_ave(r,cyc_numb)/pi/k/(Ar1_d(r)*Af)/(Ar1_d(r)*Af))).^2+(2*virial1m_ave(r,cyc_numb)/k/(Ar1_d(r)*Af)/(Ar1_d(r)*Af)).^2)));
                    Phase_lag1m(r,cyc_numb)=(atan( ((1/Q)+(dissipation1m_ave(r,cyc_numb)/pi/k/(Ar1_d(r)*Af)/(Ar1_d(r)*Af))) ./(-2*virial1m_ave(r,cyc_numb)/k/(Ar1_d(r)*Af)/(Ar1_d(r)*Af)))*180/pi);
                    if Phase_lag1m(r,cyc_numb)<0
                        Phase_lag1m(r,cyc_numb)=Phase_lag1m(r,cyc_numb)+180;
                    end

                    dissipation2m(r,cyc_numb) = -(per2/per)*trapz(time_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:))),force_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:))).*qdf2_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:)))...
                        .*coef_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:)))); % unit N.nm (nJ)
                    virial2m(r,cyc_numb)= (per2/per)*(1/per2)* trapz(time_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:))),force_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:))).*qf2_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:)))...
                        .*coef_tot_exp(cyc_numb,1:nnz(time_tot_exp(cyc_numb,:)))); % unit N.nm (nJ)  
                    dissipation2m_ave(r,cyc_numb) = mean(dissipation2m(r,1:cyc_numb));    % unit N.nm (nJ)
                    virial2m_ave(r,cyc_numb)= mean(virial2m(r,1:cyc_numb));    % unit N.nm (nJ)
                    Ar2m(r,cyc_numb)=(1/Q2)/sqrt(((1/Q2)+(dissipation2m_ave(r,cyc_numb)/pi/k2/A2(r,j)/A2(r,j))).^2+(-2*virial2m_ave(r,cyc_numb)/k2/A2(r,j)/A2(r,j)).^2);
                    Phase_lag2m(r,cyc_numb)= atan(((1/Q2)+(dissipation2m_ave(r,cyc_numb)/pi/k2/A2(r,j)/A2(r,j)))./(-2*virial2m_ave(r,cyc_numb)/k2/A2(r,j)/A2(r,j)))*180/pi;
                    if Phase_lag2m(r,cyc_numb)<0
                        Phase_lag2m(r,cyc_numb)=Phase_lag2m(r,cyc_numb)+180;
                    end  

                    cyc_numb=cyc_numb+1;                        
                end

                %% Conv. check

                if ( cyc_numb>2)
                    if or(and(Ar2m(r,cyc_numb-1)>Ar2m(r,cyc_numb),Ar2m(r,cyc_numb-1)>Ar2m(r,cyc_numb-2)),and(Ar2m(r,cyc_numb-1)<Ar2m(r,cyc_numb),Ar2m(r,cyc_numb-1)<Ar2m(r,cyc_numb-2)))
                        Ar2m_Ext(r,b)=Ar2m(r,cyc_numb-1);
                        if b>1
                            Ar2m_Ext_ave(r,b-1)=mean(Ar2m_Ext(r,b-1:b));
                            Ar2m_Ext_cyc(r,b-1)=cyc_numb-1;
                        end
                        b=b+1;
                    end
                    if or(and(Phase_lag2m(r,cyc_numb-1)>Phase_lag2m(r,cyc_numb),Phase_lag2m(r,cyc_numb-1)>Phase_lag2m(r,cyc_numb-2)),and(Phase_lag2m(r,cyc_numb-1)<Phase_lag2m(r,cyc_numb),Phase_lag2m(r,cyc_numb-1)<Phase_lag2m(r,cyc_numb-2)))
                        Phase_lag2m_Ext(r,bb)=Phase_lag2m(r,cyc_numb-1);
                        if bb>1
                            Phase_lag2m_Ext_ave(r,bb-1)=mean(Phase_lag2m_Ext(r,bb-1:bb));
                            Phase_lag2m_Ext_cyc(r,bb-1)=cyc_numb-1;
                        end
                        bb=bb+1;
                    end    
                end

                cyc_numb=cyc_numb+1;

                if (cyc_numb<=n_c)&&((nnz(Ar2m_Ext_ave(r,:))<5)||(nnz(Phase_lag2m_Ext_ave(r,:))<5))      % check for the convergence
                    condition3=1;

                elseif (cyc_numb>n_c) 
                    condition3=0;
                    converged_cycle=n_c-1;

                elseif (n_Phase_lag2m_Ext_ave_old~=nnz(Phase_lag2m_Ext_ave(r,:)))
                    if (nnz(Ar2m_Ext_ave(r,:))>conv_lim)&& (nnz(Phase_lag2m_Ext_ave(r,:))>conv_lim)
                        if any(abs(diff(Ar2m_Ext_ave(r,nnz(Ar2m_Ext_ave(r,:))-conv_lim:nnz(Ar2m_Ext_ave(r,:)))))>AR2_tol)+...
                                any(abs(diff(Phase_lag2m_Ext_ave(r,nnz(Phase_lag2m_Ext_ave(r,:))-conv_lim:nnz(Phase_lag2m_Ext_ave(r,:)))))>PL2_tol)==0
                            condition3=0;
                            for converged_cycle=min(max(Ar2m_Ext_cyc(r,:)),max(Phase_lag2m_Ext_cyc(r,:))):-1:1
                                if any(ismember(Ar2m_Ext_cyc(r,:),converged_cycle))*any(ismember(Phase_lag2m_Ext_cyc(r,:),converged_cycle))==0
                                    if converged_cycle>75
                                        break
                                    end
                                end
                            end

                        end                           
                    end

                    n_Phase_lag2m_Ext_ave_old=nnz(Phase_lag2m_Ext_ave(r,:));
                    n_Ar2m_Ext_ave_old=nnz(Ar2m_Ext_ave(r,:));

                end

%                     condition3=(cyc_numb<n_c-5);
%                     converged_cycle=n_c-5;
            end

            Ar1(r,i)= Ar1m(r,converged_cycle+1);
            Phase_lag1(r,i)= Phase_lag1m(r,converged_cycle+1);                  
            Ar2(r,j)= Ar2m(r,converged_cycle+1);
            Phase_lag2(r,j)= Phase_lag2m(r,converged_cycle+1);

            condition2= (abs(Ar2(r,j)-(A2(r,j)/Af2))>tol_Ar_2)&&(dA2>tol_Ar_2^2);

            table(converged_cycle, r, i, j, A2(r,j), Ar2(r,j)-(A2(r,j)/Af2), signchangeA2, Zguess(r,i),(Ar1(r,i)-Ar1_d(r)), Ar1_d(r)*Af, signchangeZ)

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
            E_ts_2_final(r)=dissipation2m_ave(r,converged_cycle)*6.242e9;
            V_ts_2_final(r)=virial2m_ave(r,converged_cycle)*6.242e9;
            E_ts_1_final(r)=dissipation1m_ave(r,converged_cycle)*6.242e9;
            V_ts_1_final(r)=virial1m_ave(r,converged_cycle)*6.242e9;
            errorfinal_1(r)=Ar1(r,i-1)-Ar1_d(r);
            errorfinal_2(r)=Ar2(r,j-1)-(A2(r,j-1)/Af2);
            it_i(r)=i-1;
            it_j(r)=j-1;
        end

        condition1=(abs(Ar1(r,i-1)-Ar1_d(r))>tol_Ar_1) && (dZ>tol_Ar_1^2);

        %%
        dissipation2m_1=zeros(1,n_c);
        virial2m_1=zeros(1,n_c);

        if ishandle(3)
            print('-f3','-dpng','-r1000',[address 'diss.jpg'])
            print('-f4','-dpng','-r1000',[address 'Virial.jpg'])
        end

    end 

 end   

%% plot area
S_M_P=0;

figure(2)
subplot(2,2,1)
plot(Z_final,Ar1_d,'LineWidth',2)
if S_M_P==1
    hold on
    plot(z_sm_1,Ar1_d,'o--r')
    legend('Bimodal','Single mode');
end
set(get(gca,'XLabel'),'String','Z (nm)','LineWidth',2)
set(get(gca,'YLabel'),'String','A_1/A_0_1','LineWidth',2);

subplot(2,2,3)
plot(Ar1_d,Phase_lag1_final,'LineWidth',2)
if S_M_P==1
    hold on
    plot(Ar1_d,Phase_lag_sm_1,'o--r')
    legend('Bimodal','Single mode');
end
set(get(gca,'XLabel'),'String','A_1/A_0_1','LineWidth',2)
set(get(gca,'YLabel'),'String','{\phi}_1 (degree)','LineWidth',2);
subplot(2,2,2)
plot(Ar1_d, A2_final/Af2,'LineWidth',2)
set(get(gca,'XLabel'),'String','A_1/A_0_1','LineWidth',2)
set(get(gca,'YLabel'),'String','A_2/A_0_2','LineWidth',2);

subplot(2,2,4)
plot(Ar1_d,Phase_lag2_final,'LineWidth',2)
set(get(gca,'XLabel'),'String','A_1/A_0_1','LineWidth',2)
set(get(gca,'YLabel'),'String','{\phi}_2 (degree)','LineWidth',2);

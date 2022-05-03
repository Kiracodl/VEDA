%these are the runs at 30 kHz and 1 N/m
clear
close all
rappture_reader_single( 'foo.xml');
%%


figureset(4,5,1 )

red_ndx = find( (Topog_x < 20) & (Topog_x > 10));
blue_ndx = setdiff( 1:length(Topog_x), red_ndx);


subplot(3,1,1);
plot(  Topog_x(red_ndx), Topog_y(red_ndx), 'r',  Topog_x(blue_ndx), Topog_y(blue_ndx), 'b.', 'LineWidth',2)
grid on
title( 'b) Measured Topography');
xlabel('X Position (nm)'); ylabel('Height (nm)');

subplot(3,1,2);
plot(Indent_x(red_ndx), Indent_y(red_ndx), 'r', Indent_x(blue_ndx), Indent_y(blue_ndx), 'b.','LineWidth', 2);
xlabel('X distance (nm)'); ylabel('Indentation (nm)');
title('c) Indentation');
grid on;

subplot(3,1,3);
plot(  LockinPhase_x(red_ndx), LockinPhase_y(red_ndx),'r', LockinPhase_x(blue_ndx), LockinPhase_y(blue_ndx),'b.', 'LineWidth', 2);
title('d) 1st Harmonic Phase');
xlabel('X Position (nm)'); 
ylabel('Phase (deg)');
ylim([65 67])
grid on
print -dpdf ../../doc/Manual/AMS_Ex2_H_P.pdf

%%
figure('PaperPositionMode', 'auto');

subplot(2,1,1);
plot( ErrorZ_x, ErrorZ_y, 'LineWidth', 2);
xlabel('X distance (nm)'); ylabel('Measurement error (nm)');

subplot(2,1,2);
plot( AmpDist_x, AmpDist_y, 'LineWidth', 2)
xlabel('X distance (nm)'); ylabel('Amplitude error (nm)');

print -dpdf ../../doc/Manual/AMS_Ex2_Merr_Aerr.pdf

%%
figure('PaperPositionMode', 'auto','Position', [0 500 900 300]);

subplot(1,3,1);
plot( MeanForce_x, MeanForce_y, 'LineWidth', 2);
xlabel('X distance (nm)'); ylabel('Mean interaction force (nN)');

subplot(1,3,2);
plot( FPeakAtt_x, FPeakAtt_y, FPeakRep_x, FPeakRep_y, 'LineWidth', 2);
legend('Attractive', 'Repulsive', 'Location', 'NorthWest');
xlabel('X distance (nm)'); ylabel('Peak interaction force (nN)');
ylim([-30 80])

subplot(1,3,3);
plot( Pts_x, Pts_y);
xlabel('X distance (nm)'); ylabel('Power dissipated (eV / drive cycle)');


print -dpdf ../../doc/Manual/AMS_Ex2_Mf_Pf_Wdot.pdf

%%
%old figure
%figure('PaperPositionMode', 'auto');
%plot(Indent_x, Indent_y, 'LineWidth', 2);
%xlabel('X distance (nm)'); ylabel('Indentation (nm)');
%print -dpdf ../../doc/Manual/AMS_Ex2_I.pdf


%%

clear
opt = 2;
if (opt == 1)
    E1 = 4e9; E2 = E1;
    eta1 =  400;
    eta2 =  4000;
else
    E1 = 0.1e9; E2 = E1;
    eta1 =  10;
    eta2 =  100;
end
    
nu = 0.3;

G1 = E1 / ( 2 * (1+nu));
G2 = E2 / ( 2 * (1+nu));

f = dlogspace( 1e3, 1e7, 600);
omega = f * 2 * pi;
tau1 = eta1 / (G1+G2);
tau2 = eta2 / (G1+G2);
Gstorage1 = (G1 * G2) / (G1+G2)  + (G1^2) / (G1+G2) * (tau1^2 * omega.^2) ./ (1+tau1^2 * omega.^2);
Gloss1    =                        (G1^2) / (G1+G2) * (tau1   * omega   ) ./ (1+tau1^2 * omega.^2);
Gstorage2 = (G1 * G2) / (G1+G2)  + (G1^2) / (G1+G2) * (tau2^2 * omega.^2) ./ (1+tau2^2 * omega.^2);
Gloss2    =                        (G1^2) / (G1+G2) * (tau2   * omega   ) ./ (1+tau2^2 * omega.^2);


ff1 =  30e3;
ff2 = 200e3;
figureset(6,4); 
subplot(3,1, 1:2)
semilogx( f, Gstorage1/1e6, 'r', f, Gloss1/1e6, 'r--', f, Gstorage2/1e6, 'b',f, Gloss2/1e6, 'b--', [ff1 ff1], [0 4e7]/1e6, 'k--', [ff2 ff2], [0 4e7]/1e6, 'k--', 'LineWidth', 2)
lh = legend('G_{storage}, \eta= 10', 'G_{loss}, \eta= 10', 'G_{storage}, \eta= 100', 'G_{loss}, \eta= 100', 'Drive Freq', 'Location', 'East');
set(lh, 'FontSize', 9);
xlabel('Freq (Hz)');
ylabel('Shear Modulus (MPa)');

subplot(3,1, 3)
semilogx( f, Gloss1 ./ Gstorage1, 'r', f, Gloss2 ./ Gstorage2, 'b', [ff1 ff1], [0 1], 'k--',[ff2 ff2], [0 1], 'k--', 'LineWidth', 2)
legend('\eta=10','\eta=100')
xlabel('Freq (Hz)');
ylabel('Loss Tangent')
ylim([0 0.5])
print -dpdf scanbasic_ex2_dynamic_modulus.pdf

%%

figureset(2.5, 1.5)
semilogx( f, Gloss1 ./ Gstorage1, 'r', f, Gloss2 ./ Gstorage2, 'b', [ff2 ff2], [0 1], 'k--', 'LineWidth', 2)
xlabel('Freq (Hz)');
ylabel('Loss Tangent')
ylim([0 0.4])
xlim([1e4 1e7])
grid on;
print -dpdf  scanbasic_ex2_tandel_small.pdf

%%

%for dow presentation
ff = 200e3;
figureset(5,4,1); 
subplot(3,1, 1:2)
semilogx( f, Gstorage1/1e6, 'r', f, Gstorage2/1e6, 'b',f, Gloss1/1e6, 'r--', f, Gloss2/1e6, 'b--', [ff ff], [0 4e7]/1e6, 'k--', 'LineWidth', 2)
legend('\eta=10','\eta=100', 'Location', 'East')
xlabel('Freq (Hz)');
ylabel('Shear Modulus (MPa)');

subplot(3,1, 3)
semilogx( f, Gloss1 ./ Gstorage1, 'r', f, Gloss2 ./ Gstorage2, 'b', [ff ff], [0 1], 'k--', 'LineWidth', 2)
legend('\eta=10','\eta=100')
xlabel('Freq (Hz)');
ylabel('Loss Tangent')
ylim([0 0.5])
print_png('scanbasic_ex2_dynamic_modulus_png.png', 100);


%%
%for dow training tutorial.  kelvin-voigt
f = dlogspace( 1e2, 1e8, 100);
omega = f * 2 * pi;

Gstorage = 100 * ones(size(f));
Gloss    = f / 5e5;
tandel   = Gloss ./ Gstorage;

ff = 30e3;
figureset(6,4); 
subplot(3,1, 1:2)
semilogx( f, Gstorage, 'r', f, Gloss, 'LineWidth', 2)
legend('G_{storage}', 'G_{loss}', 'Location', 'North')
xlabel('Freq (Hz)');
ylabel('Shear Modulus (MPa)');
grid on
set(gca, 'XTickLabel', []);

subplot(3,1, 3)
semilogx( f, tandel,  'LineWidth', 2)
ylabel('Loss Tangent')
ylim([0 2])
xlabel('Freq (Hz)');
grid on

%%
%these are the runs at 200 kHz and 40 N/m

clear
close all
rappture_reader_multi('foo', 1:2, '.xml');
%%


figureset(3.7,3,1 )

red_ndx = find( (Topog_x(:,1) < 20) & (Topog_x(:,1) > 10));
blue_ndx = setdiff( 1:length(Topog_x(:,1)), red_ndx);


subplot(2,1,1);
plot(  LockinPhase_x(red_ndx,1), LockinPhase_y(red_ndx,1),'r', LockinPhase_x(blue_ndx,1), LockinPhase_y(blue_ndx,1),'b.', 'LineWidth', 2);
grid on
title( 'Phase, 90% setpoint');
xlabel('X Position (nm)'); 
ylabel('Phase (deg)');
ylim([65 72]);

subplot(2,1,2);
plot(  LockinPhase_x(red_ndx,2), LockinPhase_y(red_ndx,2),'r', LockinPhase_x(blue_ndx,2), LockinPhase_y(blue_ndx,2),'b.', 'LineWidth', 2);
title('Phase, 30% setpoint');
xlabel('X Position (nm)'); 
ylabel('Phase (deg)');
ylim([24 28])
grid on

print -dpdf ../../doc/Manual/AMS_Ex2_higherfreq.pdf
print_png('../../doc/Manual/AMS_Ex2_higherfreq_png.png', 100)

%this just for dow presentation
figureset(3.7,3,1 )
subplot(2,1,1)
plot(  Topog_x(red_ndx,1), Topog_y(red_ndx,1),'r', Topog_x(blue_ndx,1), Topog_y(blue_ndx,1),'b.', 'LineWidth', 2);
xlabel('X Position (nm)'); 
title('Measured Topography');
ylabel('Height (nm)');
grid on

subplot(2,1,2)
plot(  Indent_x(red_ndx,1), Indent_y(red_ndx,1),'r', Indent_x(blue_ndx,1), Indent_y(blue_ndx,1),'b.', 'LineWidth', 2);
xlabel('X Position (nm)'); 
title('Indentation');
ylabel('Indentation (nm)');
grid on

print_png('../../doc/Manual/AMS_Ex2_higherfreq_topog.png', 100)

%%
%time histories.
figureset(5,4, 1);
subplot(2,2,1);
t = ((Fts1_x(:,1) - 0.6147) / (200e3)) / 1e-6;
plot( t,Fts1_y(:,1), 'r','LineWidth',2 );
xlim([0 1] / 200e3 / 1e-6 )
ylim([0 12])
xlabel('Time (\mu s)');
ylabel('Force (nN)');
grid on;
title('Setpoint 90%');

subplot(2,2,2)
[f,a] = genFFT(  Fts1_y(1:3750,1), [], 0, 750 *  200e3, [], @rectwin);
plot( f/1e6, a/max(a), 'r', 'LineWidth',2 );
xlim([0 2]);
xlabel('Frequency (MHz)')
ylabel('Force (normalized)')
ylim([0 1.1])
grid on
title('Setpoint 90%');

subplot(2,2,3);
t = (Fts1_x(:,2) - 0.592) / (200e3) / 1e-6;
plot(t ,Fts1_y(:,2),'r', 'LineWidth',2 );
xlim([0 1] / 200e3 / 1e-6 )
xlabel('Time (\mu s)');
ylim([0 12])
ylabel('Force (nN)');
grid on;
title('Setpoint 30%');

subplot(2,2,4)
[f,a] = genFFT(  Fts1_y(1:3750,2), [], 0, 750 *  200e3, [], @rectwin);
plot( f/1e6, a/max(a), 'r','LineWidth',2 );
xlim([0 2]);
xlabel('Frequency (MHz)')
ylim([0 1.1])
grid on
title('Setpoint 30%');

print_png(' ../../doc/Manual/AMS_Ex2_timehist.pdf', 100);

print_png(' ../../doc/Manual/AMS_Ex2_timehist_png.png', 100);

%%

%for dow presentation


figureset(4,3,1);
t = ((Fts1_x(:,1) - 0.6147) / (200e3)) / 1e-6;
p1 = @(x,y) plot(x,y,'m', 'LineWidth', 2);
p2 = @(x,y) plot(x,y,'k', 'LineWidth', 2);
ax = plotyy( t,Fts1_y(:,1),t, gap1_y(:,1 ) - mean( gap1_y(1:3750,1)),p1,p2);
title('Deflection and Force History at 90% setpoint');
xlim(ax(1), [0 2] / 200e3 / 1e-6 )
xlim(ax(2), [0 2] / 200e3 / 1e-6 )
xlabel('Time (\mu s)');
ylabel(ax(1), 'Force (nN)');
ylabel(ax(2), 'Deflection (nm)');

print_png('defl_and_force.png', 100);


figureset( 2.5, 1.6, 2);
[f,a] = genFFT(  Fts1_y(1:3750,1), [], 0, 750 *  200e3, [], @rectwin);
plot( f/1e6, a, 'm','LineWidth',2 );
xlim([0 2]);
ylabel('Force (nN)');
xlabel('Frequency (MHz)')
grid on
ylim([0 3])
title('Fourier Transform of Force');
print_png('force.png', 100);

figureset( 2.5, 1.6, 3);
[f,a] = genFFT(  gap1_y(1:3750,1), [], 0, 750 *  200e3, [], @rectwin);
plot( f/1e6, a, 'k','LineWidth',2 );
xlim([0 2]);
xlabel('Frequency (MHz)')
ylabel('Deflection (nm)');
ylim([0 30])
grid on
title('Fourier Transform of Deflection');
print_png('defl.png', 100);
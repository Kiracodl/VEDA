clear
close all
rappture_reader_single( 'AMS_Ex1.xml')

figure('PaperPositionMode', 'auto');
subplot(2,1,1);
plot(  Topog_x, Topog_y, Hsample_x, Hsample_y,'r--',  'LineWidth', 2);
legend('Sample Height', 'Measured Topography', 'Location', 'North');
xlabel('X Position (nm)'); ylabel('Height (nm)');

subplot(2,1,2);
plot(  Phase_x, Phase_y, LockinPhase_x, LockinPhase_y, 'LineWidth', 2);
legend('1st Harmonic Phase', 'Lockin 1st Harmonic Phase', 'Location', 'North');
xlabel('X Position (nm)'); ylabel('Phase (deg)');

print -dpdf ../../doc/Manual/AMS_Ex1_H_P.pdf

%%
figure('PaperPositionMode', 'auto');

subplot(2,1,1);
plot( ErrorZ_x, ErrorZ_y, 'LineWidth', 2);
xlabel('X distance (nm)'); ylabel('Measurement error (nm)');

subplot(2,1,2);
plot( AmpDist_x, AmpDist_y, 'LineWidth', 2)
xlabel('X distance (nm)'); ylabel('Amplitude error (nm)');

print -dpdf ../../doc/Manual/AMS_Ex1_Merr_Aerr.pdf

%%
figure('PaperPositionMode', 'auto');

subplot(1,2,1);
plot( MeanForce_x, MeanForce_y, 'LineWidth', 2);
xlabel('X distance (nm)'); ylabel('Mean interaction force (nN)');

subplot(1,2,2);
plot( FPeakAtt_x, FPeakAtt_y, FPeakRep_x, FPeakRep_y, 'LineWidth', 2);
legend('Attactive', 'Repulsive');
xlabel('X distance (nm)'); ylabel('Peak interaction force (nN)');

print -dpdf ../../doc/Manual/AMS_Ex1_Mf_Pf.pdf

%%
figure('PaperPositionMode', 'auto');
plot(Indent_x, Indent_y, 'LineWidth', 2);
xlabel('X distance (nm)'); ylabel('Indentation (nm)');
print -dpdf ../../doc/Manual/AMS_Ex1_I.pdf
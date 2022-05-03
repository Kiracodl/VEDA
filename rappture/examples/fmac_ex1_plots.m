rappture_reader_single('fm.xml')
%%

figureset(4,3);
plot( LockinPhase_x, LockinPhase_y, 'LineWidth', 2);
xlabel('Min gap (nm)');
ylabel('Phase (deg)');
grid on;
ylim([89.98 90.01]);

print -dpdf ../../doc/Manual/FMAC_Ex1_phase.pdf

%%

figure;
p = @(x,y) plot(x,y, 'LineWidth', 2);
ax = plotyy( DriveFreq_x, DriveFreq_y, DriveAmp_x, DriveAmp_y ./ DriveAmp_y(1), p);
ylabel(ax(1), 'Freq shift (kHz)');
ylabel(ax(2), 'Drive Amp');
xlabel('Min gap (nm)');
print -dpdf ../../doc/Manual/FMAC_Ex1_FMFeedback.pdf
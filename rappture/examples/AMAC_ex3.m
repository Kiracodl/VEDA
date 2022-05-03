clear
rappture_reader_multi('foo', 1:2, '.xml')

figureset(4,4);
subplot(2,1,1)
plot( Amp_x, Amp_y, 'LineWidth', 2)
text(0, 70, 'a)')
ylabel('Amplitude (nm)');
xlim([0 65])
ylim([0 65])
subplot(2,1,2)
plot( Amp_x, Phase_y, 'LineWidth', 2);
ylabel('Phase (deg)');
xlabel('Z distance (nm)');
text(0, 110, 'b)')
xlim([0 65])
legend('With Viscoelasticity', 'Without Viscoelasticity', 'Location', 'NorthWest')
print -dpdf ../../doc/Manual/DACex3_A_P.pdf

%%

figureset(4,5)
subplot(2,1,1);
plot(gap1_x, gap1_y, 'LineWidth', 2)
xlim([0.53 0.65])
ylim([-2 4]);
xlabel('Time (t/T)');
ylabel('Tip-Sample gap (nm)');
legend('With Viscoelasticity', 'Without Viscoelasticity', 'Location', 'North')
grid on

subplot(2,1,2)
plot(Fts1_x, Fts1_y, 'LineWidth', 2)
xlim([0.53 0.65])
ylim([0 18]);
xlabel('Time (t/T)');
ylabel('Tip-Sample Force (nN)');
grid on

print -dpdf ../../doc/Manual/DACex3_D_TS.pdf

%%

figureset(7, 3);
subplot(1,3,1);
plot( MeanForce_x, MeanForce_y, 'LineWidth', 2);
ylabel('Mean Interaction Force (nN)');
xlim([0 60])
title('a)')

subplot(1,3,2);
plot(FPeakRep_x, FPeakRep_y, 'LineWidth', 2);
ylabel('Peak Interaction Force (nN)');
xlim([0 60])
title('b)')
lh = legend('Viscoelasticity', 'Without', 'Location', 'SouthWest')
set(lh, 'FontSize', 9);


subplot(1,3,3);
plot( Pts_x, Pts_y, 'LineWidth', 2);
ylabel('Power Dissipation (eV/drive cycle)');
xlim([0 60])
title('c)')

print -dpdf ../../doc/Manual/DACex3_Mf_Pf_Wdot.pdf

%%
figureset(6, 3);
subplot(1,2,1);
plot( Indent_x, Indent_y, 'LineWidth', 2);
ylabel('Indentation (nm)');
xlim([0 60])
title('a)')

subplot(1,2,2);
plot(tcontact_x, tcontact_y, 'LineWidth', 2);
ylabel('Contact Time (us)');
xlim([0 60])
title('b)')
legend('With Viscoelasticity', 'Without Viscoelasticity', 'Location', 'NorthEast')

print -dpdf ../../doc/Manual/DACex3_I_C.pdf
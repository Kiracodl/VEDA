rappture_reader_multi('run', 1:2 , 'xml');
figureset(6,4);
subplot(1,2,1);
plot( LockinR_x, LockinR_y, 'LineWidth', 2);
legend('T_{lockin} = 500 \mus', 'T_{lockin} = 5 ms', 'Location', 'South');
ylabel('A_{1,1} (nm)')
xlabel('Z-distance (nm)')
text(1, 10.5, '(a)', 'FontSize', 16);
title('Amplitude @ 1st Drive Freq');


subplot(1,2,2);
plot(LockinAmp2_x, LockinAmp2_y, 'LineWidth', 2);
legend('T_{lockin} = 500 \mus', 'T_{lockin} = 5 ms', 'Location', 'South');
ylabel('A_{1,2} (nm)')
xlabel('Z-distance (nm)')
text(1, 1.05, '(b)', 'FontSize', 16);
title('Amplitude @ 2nd Drive Freq');

print -dpdf ../../doc/Manual/ADACex2_AT.pdf
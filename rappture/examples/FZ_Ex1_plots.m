rappture_reader_multi('FZ_Ex1_', 1:2, 'xml');

%%

figure('Position', [0 500 350 300], 'PaperPositionMode', 'auto')

plot( DeflZ_x(:,1), DeflZ_y(:,1), 'ro', DeflZ_x(:,2), DeflZ_y(:,2), 'b*','LineWidth', 2);
legend('Approach','Retract', 'Location', 'NorthEast');
%ylim([0 25]); xlim([0 25]);
xlabel('Z-distance (nm)');
ylabel('Observed Deflection (nm)');

print -dpdf ../doc/Manual/FZex1_deflZ.pdf
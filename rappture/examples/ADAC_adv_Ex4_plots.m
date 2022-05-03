rappture_reader_multi('output', 1:2, 'xml');

%%

figure('Position', [0 500 650 300])
subplot(1,2,1);
plot( Amp_x(:,1), Amp_y(:,1), 'r--', Amp_x(:,2), Amp_y(:,2), 'b','LineWidth', 2);
legend('E_{sample} = 1 GPa','E_{sample} = 5 GPa', 'Location', 'SouthEast');
ylim([0 25]); xlim([0 25]);
xlabel('Z-distance (nm)');
ylabel('A_{1,1} (nm)');
title('(a)');

subplot(1,2,2);
plot( Amp2_x(:,1), Amp2_y(:,1), 'r--', Amp2_x(:,2), Amp2_y(:,2), 'b','LineWidth',2);
legend('E_{sample} = 1 GPa','E_{sample} = 5 GPa', 'Location', 'North');
xlabel('A{1,1} / A{0,1} (nm)');
ylabel('A_{1,2} (nm)');
title('(c)');

print  -dpdf ~/app-adac-stable/doc/Manual/ADACex4_P.pdf

figure('Position', [0 800 650 300])
subplot(1,2,1);
plot( Phase_x(:,1), Phase_y(:,1), 'r--', Phase_x(:,2), Phase_y(:,2), 'b','LineWidth',2);
legend('E_{sample} = 1 GPa','E_{sample} = 5 GPa', 'Location', 'North');
xlabel('A{1,1} / A{0,1} (nm)');
ylabel('Phase_{1,1} (deg)');
title('(b)');

subplot( 1,2,2);
plot( Phase2_x(1:499,1), Phase2_y(1:499,1), 'r--', Phase2_x(1:499,2), Phase2_y(1:499,2), 'b','LineWidth',2);
legend('E_{sample} = 1 GPa','E_{sample} = 5 GPa', 'Location', 'North');
xlabel('A{1,1} / A{0,1} (nm)');
ylabel('Phase_{1,2} (deg)');
title('(d)');

print  -dpdf ~/app-adac-stable/doc/Manual/ADACex4_P.pdf

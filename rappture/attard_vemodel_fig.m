clear
close all

G0   = 10
Ginf =  2.5
tau  = 1;
omega = logspace( log10(1e-2), log10(1e2), 100)

t = linspace(0, 10, 100);

figure
subplot(1,3,1)
i = imread('viscoelastic_elements.png'); %, "PixelRegion", {[200 600], [300 700]});
image(i( 20:end ,210:end,:))
axis off

subplot(1,3,2)
psi = Ginf + (G0 - Ginf) * exp(-t / tau)
plot(t, psi, 'LineWidth', 2)
title('psi(t)');
xlabel('Time (s)');
ylim([0 10])
text( 8, 3, 'Ginf')
text( 1, 9, 'G0')
set(gca, 'XTick', []);
set(gca, 'YTick', []);


subplot(1,3,3)
Gst = Ginf + (G0-Ginf) * tau^2 * omega.^2 ./ ( 1 + tau^2 * omega.^2);
Gloss = (G0 - Ginf) * tau * omega ./ (1 + tau^2 * omega.^2);
semilogx(omega, Gst, 'LineWidth', 2); hold on;
semilogx(omega, Gloss, 'LineWidth', 2)
xlabel('Freq');
legend('Storage', 'Loss', 'Location', 'East')
%ylim([4 11])
text( 1e-2, 3, 'Ginf')
text( 20, 9.5, 'G0')
set(gca, 'XTick', []);
set(gca, 'YTick', []);

print('-dpng', '-r75', 'three_element_model.png')
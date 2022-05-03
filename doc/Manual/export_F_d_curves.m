
option = 3


%list = {'Chadwick_F-d',  'Hertz_F-d',  'LinAR_F_d', 'DMTDLVO_F-d', 'DMT_F-d', 'JKR_F-d',  'LinCon_F-d', 'Chadwick_Capillary_F-d', 'DMT_Capillary_F-d', 'Hertz_Capillary_F-d', 'LinearAR_Capillary_F-d', 'LinCon_Capillary_F-d'};
list = {'Hertz_Capillary_F-d'};
for i = 1:length(list)    
    s =cell2mat(list(i))
    openfig([s '.fig']);
    set( gcf, 'WindowStyle', 'normal');
    
    if (option == 1) 
        set(gcf, 'PaperSize', [8 6]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperPosition', [0 0 8 6]);
        print('-dpdf', [s '.pdf']);
    elseif (option == 2)
        set(gcf, 'PaperSize', [4 3]);
        print('-dpng', '-r60', [s '.png']);  %used -r70 initially.  images a little too big.  try 60.
    elseif (option == 3)        
        set(gcf, 'PaperPositionMode', 'manual');
        px = 4 * 0.4; py = 3 * 0.4;
        set(gcf, 'PaperSize', [px py]);
        set(gcf, 'PaperUnits', 'inches');
        set( gcf, 'PaperPosition', [0 0 px py]);
        %set( gcf, 'Position', 95 * [0 0 px py]);
        set( gca, 'XTick', [0]);
        set( gca, 'YTick', [0]);
        set( gca, 'XTickLabel', []);
        set( gca, 'YTickLabel', []);       
        h = get(gca, 'Title')
        set(h, 'FontSize', 10);
        set(h, 'FontName', 'Helvetica-Narrow')
        h = xlabel('Gap  d')
        set(h, 'FontSize', 10);
        ylabel('Force  F_{ts}', 'FontSize', 10);
        legend off
        
        hold on;
        plot( [-1 1], [0, 0], 'k--', 'LineWidth', 1)
        plot( [0 0], [-10000, 100], 'k--', 'LineWidth', 1)

        
        print_png([s, '.png'], 300, 'transparent')            
    end
    
    %close
end
    
    
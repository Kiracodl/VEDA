clear
fn = 'Hertz'
rappture_reader_single([fn '.xml'])
%%
option = 1;

if (option == 1)
    %pdf, for manual
    figureset(4,3);
else
    %gif, thumbnail for gui
    figureset( 1.6, 1);
end
    
plot( FtsGap_x, FtsGap_y, 'r', 'LineWidth', 2)

%Ron says the lack of labels in the VEDA manual confused him.
%so put the labels in
%set( gca, 'XTickLabel', []);
%set( gca, 'YTickLabel', []);
%set( gca, 'XTick', [0]);
%set( gca, 'YTick', [0]);

hold on;
plot( [-10 10], [0, 0], 'k--', 'LineWidth', 1)
plot( [0 0], [-10000, 100], 'k--', 'LineWidth', 1)

xlabel('Gap d (nm)')
ylabel('Force  F_{ts} (nN)')

if (strcmp( fn, 'osc'))
    if (option == 1)
        title('Oscillatory')    
    end
    xlim([-0.1 1.5])
    ylim([-60 100])
elseif (strcmp( fn, 'Fts_visco'))
    text(0.1, 1.8, {'Kelvin-','Voigt'})
    xlim([-1.1 1])
    ylim([-0.5 5])
elseif (strcmp( fn , 'morse'))
    if (option == 1)
        title('Morse')
    else
        text(0.6, 16, 'Morse', 'FontSize', 12)
    end
    
    xlim([-0.1 1.3])
    ylim([-3 20])
elseif (strcmp( fn , 'electrostatic'))
    if (option == 1)
        title('Electrostatic')    
    end
    xlim([-0.1 2])
    ylim([-3000 200])
elseif (strcmp( fn, 'Hertz') )
    xlim([-2 2])
    ylim([-0.1 1.5])
end

if  (option == 1)    
    print( '-dpdf', [fn, '_F-d.pdf'])
else    
    %export_fig([fn '_F-d.png'], '-r100')
   % print('-dpng',  [fn '_F-d.png'] );
    print_png([fn, '_F-d.png'], 300, 'transparent')    
end
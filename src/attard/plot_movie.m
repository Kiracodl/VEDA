data = importdata('surface_movie.txt');
%%

include_tip = 1;
want_all_three_plots = 0;

opt = 3


if (opt == 3)
    m = size(data,1);
    if (include_tip == 1)
        tip = data(:,1);
        d    = data((1:3:m)    , 2:end);
        pr   = data((2:3:(m+1)), 2:end);
        udot = data((3:3:(m+2)), 2:end );
    else
        d    = data((1:3:m)    , :);
        pr   = data((2:3:(m+1)), :);
        udot = data((3:3:(m+2)), :);
    end

else
    m = size(data,1);
    d = data(1:2:m, :);
    pr = data(2:2:(m+1), :);
    udot = 0;
end

R = 10; N =200;
rr = [quad_space(0,  R, 0.7 * N , 2), linspace(1.01*R, 5*R, 0.3 * N) ];   
%%


figure;
    j = 1;
   %for i = 1:2:size(d,1)
   for i = 320:1:400
        if (want_all_three_plots)
            subplot(opt,1,1)
        end
        plot( rr, d(i, :), 'k+', -rr, d(i,:),'k' )
        
        %ylim([min(min(d)), max(max(d)) ]);
        ylim([-0.1, 0.4 ]);

        %xlim([-30 30])
        xlim([-10 30]);
        
        if (include_tip)
            hold on;
            plot(0, tip(i), 'r+');
            hold off;
        end
        
        ylabel('Surface Displacement (nm)');
        title(num2str(i) );
        
        if (want_all_three_plots)
            subplot(opt,1,2)
            plot( rr, pr(i,:),'k+', -rr, pr(i,:), 'k' );
            ylim([min(min(pr)), max(max(pr)) ]);
            xlim([-30 30])
            ylabel('Pressure');
        
            if (opt == 3)
                subplot(3,1,3)
                plot( rr, udot(i,:) );
                ylim([min(min(udot)), max(max(udot)) ]);
                ylabel('udot');
            end
        end
         
%         
%        frames(j) = getframe(gcf); j=j+1;
        
        pause(0.2)
    end

%%

%interactive

figure;

i = 330;
while (i>0)    
        plot( rr, d(i, :), 'k+', -rr, d(i,:),'k' )

        hold on
        plot( rr, d(i-1, :), 'g')
        hold off

        %ylim([min(min(d)), max(max(d)) ]);
        ylim([-0.1, 0.4 ]);

        %xlim([-30 30])
        xlim([-10 30]);
        
        ylabel('Surface Displacement (nm)');
        title(num2str(i) );
                        
        ii = input('foo', 's');
        if ( ii == 's')
            i =0;
        elseif (ii == 'k')
            i = i+1;
        else
            i = i-1;
        end
    end


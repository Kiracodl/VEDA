%function prony_fit(filename, debugging, datasrc, data_in, N_prony)
%
% main purpose is usage inside rappture tool.  in which case will be called with filename as an argument and no other arguments
% can also be called directly from within matlab.  In that case, specify datasrc =1, and make data_in a matrix
% with frequency in col 1, e_storage in col 2, and e_loss in col 3.

function prony_fit(filename, debugging, datasrc, data_in, N_prony)

if (nargin < 3)
    debugging = 0;
    datasrc = 0;
end

if (nargin < 4)
  data_in = [];
end

if (datasrc == 0)
    %read from rappture
    [libHandle,err] = rpLib( filename );
    
    N_om = rpLibGetDouble( libHandle, 'input.phase(input).number(n_dma).current');
    N_prony = rpLibGetDouble( libHandle, 'input.phase(input).number(n_prony).current');
    
       
    f        = strread( rpLibGetString( libHandle, 'input.phase(input).string(freq_in).current'),   '%f', N_om, 'delimiter', ',');
    estorage = strread( rpLibGetString( libHandle, 'input.phase(input).string(G_storage).current'), '%f', N_om, 'delimiter', ',')';
    eloss    = strread( rpLibGetString( libHandle, 'input.phase(input).string(G_loss).current'),    '%f', N_om, 'delimiter', ',')';

    tandel = eloss ./ estorage;
    
    tauopt = 1;
    bound_option = 1;
elseif (datasrc == 1)
    %data passed as a matrix;
    
      f = data_in(:,1);
      estorage = data_in(:,2);
      eloss = data_in(:,3);
      
      N_om = length(f);
                 
      tandel = eloss ./ estorage;
    
      tauopt = 1;
    
  bound_option = 2;
elseif ((datasrc == 4) || (datasrc ==4.5))
    %this and all the following data sources are some built-in data for testing
    
    fxx = 10* [1,  5, 10,   100, 500, 1000, 2000, 1e4, 1e5, 2e5, 1e6];
    
    if (datasrc == 4.5)
        fx = (  10.^((log10(fxx)-1) * (3/6) + 2.5));
    else
        fx = fxx;
    end
    
    estoragex = 0.1 + [-0.01,  0,  0.005, 0.015, 0.03,   0.05, 0.07    0.085, 0.095, 0.1, 0.11];
    tandelx = 9*      [0.005, 0.007   0.009, 0.015, 0.022       0.03, 0.022    0.015, 0.009, 0.007, 0.005];
    
    f = dlogspace(min(fx), max(fx), 201);
    estorage =  interp1( log10(fx), estoragex, log10(f), 'linear');
    tandel =    interp1( log10(fx), tandelx,   log10(f), 'linear');
    
    eloss = tandel .* estorage;
    N_prony = 21;    
    
    tauopt =0;
    bound_option = 1;
 elseif (datasrc == 5)
    %for testing 3 elm fit
    f = [1 10 100 1000 10000];
    estorage = [1 1.1 1.5 1.9 2];
    eloss = [0.1 0.2 0.5 0.2 0.1];
    tandel = eloss ./ estorage;
    N_prony = 1;    
    
    tauopt =1;
    bound_option = 1;
elseif (datasrc == 6)
    %veda tool defaults
    f = [1 10 100 1000 10000];
    estorage = [1 1.1 1.5 1.9 2];
    eloss = [0.1 0.2 0.5 0.2 0.1];
    tandel = eloss ./ estorage;
    N_prony = 5;    
    
    tauopt =1;
    bound_option = 1;
end

wtopt = 2;
if (wtopt == 1)
    %this might favor one or the other depending on relative magnitudes.
    weight = max(eloss) * 2;
else
    %this weights eloss and estorage about evenly.
    weight = max(eloss)  / max(estorage) ;
end


if ( N_prony == 1)
    %hidden flag, fit a three element model.  i.e. fit tau instead of
    %having it be proscribed.

    x0 = [estorage(1), 0.1, 1/(min(f)+max(f))];
    
    els = @(X, f)  els_t(X(1:2), 2*pi*f, X(3));
    est = @(X, f)  est_t(X(1:2), 2*pi*f, X(3));    
    
    fun = @(E, f) [est(E, f),  els(E, f) / weight];
    
    ub = [ 100, 100, 1/ (0.01*min(f))];
    lb = [ 0, 0,     1/ (100*max(f))];

    
     %try to work equally well with matlab as with octabe
    if (exist('octave_config_info'))
      %octave
      opt = optimset('ubound', ub', 'lbound', lb');
      [X] = nonlin_curvefit (fun, x0', f, [estorage, eloss / weight], opt)       
    else  
      %matlab
      X = lsqcurvefit(fun,x0,f,[estorage, eloss / weight], lb, ub);
    end 
    
    

else
    %this really is a real prony series
   if ((tauopt == 1) || (N_prony < 3))
        gamma = 5;
    
     tau_j = dlogspace( 1/(min(f) / gamma), 1/(max(f)*gamma), N_prony) /(2*pi) ;  

%     tau_j = dlogspace( 1/(min(f) / gamma), 1/(max(f)*gamma), N_prony)  ;  
    else
        %make sure to put one relaxation at the peak.  doesn't actually seem to work very well
        if ( floor(N_prony/2) == ceil(N_prony/2))
            %even number
            midndx = floor(N_prony / 2);
        else
            %odd
            midndx = floor(N_prony / 2)+1;
        end
        
        
        [foo, peakndx] = max( tandel );
                
        tau_j(midndx) = peakndx;
        
        tau_j(1:midndx) = dlogspace( 1/min(f),     1/f(peakndx), midndx);
        tau_j(midndx:N_prony) = dlogspace( 1/f(peakndx), 1/(max(f)*10), N_prony - midndx + 1);
    
   end



    x0 = [estorage(1), 1/N_prony * ones(1,N_prony)];
    
    if (bound_option == 1)
        %bounded optimization.
        %8/7/2012. earlier I had convinced myself that the 2pi was not needed, but I think that
        %was a mistake. this has to be in radians
        els = @(E, f)  els_t(E, 2*pi*f, tau_j);
        est = @(E, f)  est_t(E, 2*pi*f, tau_j);
    elseif (bound_option == 2)
        %no lower bound on modulli, so ignore sign.
        els = @(E, f)  els_t(abs(E), 2*pi*f, tau_j);
        est = @(E, f)  est_t(abs(E), 2*pi*f, tau_j);
    end
    

    fun = @(E, f) [est(E, f),  els(E, f) / weight];


    if (bound_option == 1)
        %8/24/2014  ub(1) is actually lower than lb(1).  not sure how that managed to actually work this whole time
        %can only test it in octave right now, don't have matlab access.  therefore leave as-is for now, and just do
        %the unbounded option in octave for now.  
        ub = [ 10 * ones(1,1+N_prony)];  
        lb = [ 0.5 * estorage(1),     zeros(1,N_prony)];
    else
        ub = [];
        lb = [];
    end

    %try to work equally well with matlab as with octabe
    if (exist('octave_config_info'))
      %octave
      opt = optimset('ubound', ub', 'lbound', lb', 'MaxFunEvals', 10000);
      [X] = nonlin_curvefit (fun, x0', f, [estorage, eloss / weight], opt)       
    else  
      %matlab
      opt = optimset('MaxFunEvals', 10000);
      X = lsqcurvefit(fun,x0,f,[estorage, eloss / weight], lb, ub, opt);
    end 
    

%8/24/2014.  these parameters were not being used, so don't calculate them.
    %for j = 1:length(X)
    %    eee(j,:) = est( [zeros(1,j-1), X(j), zeros(1,N_prony - j+1)], f);
    %    eel(j,:) = els( [zeros(1,j-1), X(j), zeros(1,N_prony - j+1)], f);
    %end
    
    
end
%add some more resolution
fup = logspace( log10(min(f)), log10(max(f)), max( 100, length(f)));

tnd = @(E, f)  els(E, f) ./ est(E, f);

%%
if (~ debugging)
    SetupGenericPlot(libHandle, 'estorage', 'G storage raw', 'Modulii', 'Freq', '(Hz)', 'Modulus', 'GPa');
    SetupGenericPlot(libHandle, 'eloss',    'G loss raw',    'Modulii', 'Freq', '(Hz)', 'Modulus', 'GPa') ;
    SetupGenericPlot(libHandle, 'tandel',   'tan del raw',   'Loss Tangent', 'Freq', '(Hz)', 'tan del', '(nd)')        ;

    SetupGenericPlot(libHandle, 'estorage_fit', 'G storage fit', 'Modulii', 'Freq', '(Hz)', 'Modulus', 'GPa')       ;
    SetupGenericPlot(libHandle, 'eloss_fit',    'G loss fit',    'Modulii', 'Freq', '(Hz)', 'Modulus', 'GPa')        ;
    SetupGenericPlot(libHandle, 'tandel_fit',   'tan del fit',   'Loss Tangent', 'Freq', '(Hz)', 'tan del', '(nd)')   ;
    
    for i = 1:length(f)
        s = sprintf('%f %f\n', f(i), estorage(i) );
        rpLibPutString( libHandle, 'output.curve(estorage).component.xy', s, 1);
        
        s = sprintf('%f %f\n', f(i), eloss(i) );
        rpLibPutString( libHandle, 'output.curve(eloss).component.xy', s, 1);
               
        s = sprintf('%f %f\n', f(i), tandel(i) );
        rpLibPutString( libHandle, 'output.curve(tandel).component.xy', s, 1);      
    end
    
    for i = 1:length(fup)
        s = sprintf('%f %f\n', fup(i), est(X, fup(i)) );
        rpLibPutString( libHandle, 'output.curve(estorage_fit).component.xy', s, 1);
        s = sprintf('%f %f\n', fup(i), els(X, fup(i) ));
        rpLibPutString( libHandle, 'output.curve(eloss_fit).component.xy', s, 1);
        s = sprintf('%f %f\n', fup(i), tnd(X, fup(i) ));
        rpLibPutString( libHandle, 'output.curve(tandel_fit).component.xy', s, 1)    ;
    end
end

%%

if (debugging)
    figure; 
    subplot(3,1,1);
    semilogx( f, estorage, fup, est(X, fup), 'LineWidth', 2);
    title('estorage')
    xlabel('Frequency (Hz)'); ylabel('GPa')
    legend('DMA data', 'Prony fit', 'Location', 'SouthEast')

    subplot(3,1,2);
    semilogx( f, eloss, fup,  els(X, fup), 'LineWidth', 2 );
    title('eloss')
    xlabel('Frequency (Hz)'); ylabel('GPa')
    legend('DMA data', 'Prony fit','Location', 'SouthEast')

    subplot(3,1,3);
    semilogx( f, tandel, fup, tnd(X, fup), 'LineWidth', 2);
    title('tandel');
    xlabel('Frequency (Hz)');
    legend('DMA data', 'Prony fit','Location', 'SouthEast')        
end
%%
%rappture string outputs here.
if (N_prony ==1)
    rslt = ['G_infinity = ' num2str(X(1)) '(GPa)' char(10) 'Modulus = ' num2str(X(2))  'GPa' char(10) 'Relaxation time = ' num2str(X(3)) '(s)'];
else
    X(1)
    s =''; s2 ='';
    for j = 2:length(X)
        s = [s sprintf('%3.2e, ', X(j))];
        s2 = [s2 sprintf('%3.2e, ', tau_j(j-1))];
    end

    rslt = ['G_infinity = ' num2str(X(1)) '(GPa)' char(10) 'Modulii = ', s, '(GPa)' char(10) 'Relaxation Times = ', s2, '(s)'];
    
    if (debugging)
        display(rslt);
    end
end


if (~ debugging)
  rpLibPutString( libHandle, 'output.string(series).current', rslt, 1)
  rpLibPutString( libHandle, 'output.string(series).about.label', 'Output Data', 1)
  rpLibPutString( libHandle, 'output.string(series).size', '80x100', 1)

  rpLibResult(libHandle);

    quit
end

end

function SetupGenericPlot(libHandle, name, label, group, xlabel, xunits, ylabel, yunits)        
  rpLibPutString(libHandle, ['output.curve(' name ').about.label'], label, 0) 
  rpLibPutString(libHandle, ['output.curve(' name ').about.group'], group, 0) 
  rpLibPutString(libHandle, ['output.curve(' name ').xaxis.label'], xlabel, 0) 
  rpLibPutString(libHandle, ['output.curve(' name ').xaxis.units'], xunits, 0) 
  rpLibPutString(libHandle, ['output.curve(' name ').xaxis.scale'], 'log', 0) 
  rpLibPutString(libHandle, ['output.curve(' name ').yaxis.label'], ylabel,0)
  rpLibPutString(libHandle, ['output.curve(' name ').yaxis.units'], yunits, 0) 
end

function rslt = est_t(E, omega, tau_j) 
 %try to work equally well with matlab as with octabe
if (exist('octave_config_info'))
    tau_j = tau_j';%the transpose is necessary in octave.  matlab doesn't need it.    
end
for j = 1:length(omega)    
    rslt(j) = E(1) + sum(E(2:end) * omega(j)^2  ./ (1./tau_j.^2 + omega(j)^2)); 
end
rslt = rslt';
end

function rslt = els_t(E, omega, tau_j)        
if (exist('octave_config_info'))
    tau_j = tau_j';%the transpose is necessary in octave.  matlab doesn't need it.    
end
for j = 1:length(omega)    
    rslt(j) = sum((E(2:end) * omega(j) ./ tau_j) ./ (1./tau_j.^2 + omega(j)^2));
end
rslt = rslt';
end

%matlab's linspace function and logspace function have differing definitions of the input values.
%this is a version of logspace that has same inputs as linspace.
function out = dlogspace( x1, x2, n)
out = logspace( log10(x1), log10(x2), n);
end

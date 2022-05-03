%
%we have had a preisach module in veda for some time, but we have not had a good way to 
%determine the weights.  This function generates the weights for a proposed hysteretic
%hydration force model.  It can be adapted to generate a preisach model for any type of force.
%
% the idea is that we first explicitly compute the approach and retract curves (force versus distance)
%for a series of different approach distances. This is called the "first order reversal curves".
%the weights are then solved such that the preisach model will reproduce these curves exactly.
%second order reversal curves (e.g. bimodal) then fall out of the model naturally.
%
%in this case we assume that the approach curves and the retract curves are both exponential functions,
%but the decay length is different on the approach and the retract.  There are five parameters:
%the approach decay length (l1), the retract decay length (l2), and approach scaling (ph), the number of terms (nside)
% and the maximum tip-sample gap considered (max_d) (any gaps above that are considered to have zero force)
%
%input is from the command line, and output goes to a file named preisach.txt
%
%nside=70 takes about 1 second to run.  nside = 80 is about 2, nside=100 is about 3.5

function preisach( l1, l2, ph ,nside, max_d, do_plots )
'starting now'
if (nargin < 6)
    do_plots = [];
end

if (isdeployed)
    l1 = str2double( l1)
    l2 = str2double( l2);
    ph = str2double(ph)
    nside = str2double(nside)
    max_d = str2double(max_d)
    do_plots = []
end

%nside = 120;

%this is the vector of the switch points
d = linspace(max_d, 0, nside);
%this is the vector of the trajectory we use to setup the coeff matrix
%one point on either side
dd = [d (d(end)-1) ] + d(end-1)/2;

%l1 = 0.25;
%l2 = 1/16;

%the forward curve is always the same 
fwd = ph * exp(- d / l1);

rev = zeros(nside, nside);

%generate reversal curves
for i = 1:nside
    rev(i, i:nside) = NaN;
    rev(i, 1:i) = fwd(i) * exp(- (d(1:i)-d(i)) /  l2);    
end

%figure; plot(d, fwd, '+-', d, rev', '+-', 'LineWidth', 2)
%xlabel('gap');
%ylabel('Force');


% we can determine the weight function uniquely from the first order reversal curves
%just need to set up the appropriate matrix and solve.

%setup the alpha and beta vectors
opt = 2;
if (opt == 1)
    %this generates the states as one big long vector
    k = 1;
    for i = 1:nside
        for j = i:nside
            alpha(k) = d(j);
            beta(k)  = d(i);
            k = k + 1;
        end
    end    
    states = zeros(size(alpha));
else
    %this generates the states as a matrix
    alpha = zeros(1, nside);
    beta = zeros(1, nside);
    for j = 1:nside
        alpha(j) = d(j);
        beta(j)  = d(j);
    end
    states = zeros(nside, nside);
end

k = 1;

coeff = zeros( nside*nside, nside*nside);

%setup the states & coeffs for the forward part
b = zeros( nside, 1);
for j = 2:nside+1
    states = update_states(dd(j), dd(j-1), states, alpha, beta);
    
    if (opt == 1)
        coeff(k,:) = states;
    else
        coeff(k,:) = reshape(states, nside*nside, 1);
    end
    
    b(k) = fwd(j-1);
    
    k = k + 1;
end

%fixme. opt=2 does not come up with a good weight matrix?
for i = 2:(nside+1)
    disp(i);  
    if (opt == 1)
        states = zeros(size(alpha));
    else
        states = zeros(nside, nside);
    end
    
    %setup the states for the forward part
    for j = 2:i
        states = update_states(dd(j), dd(j-1), states, alpha, beta);
    end
    
    %setup the states for the reverse part, with one row in
    %the coeff matrix for each step
    for j = (i-1):-1:2        
        states = update_states(dd(j), dd(j+1), states, alpha, beta);
        
        if (opt == 1)
            coeff(k,:) = states;
        else
            coeff(k,:) = reshape(states, nside*nside, 1);
        end
                
        b(k) = rev(i-1, j-1);
        
        k = k + 1;
    end
end



%fill in to a square matrix.  much faster
for i = 1:(nside^2)
    if ( nnz(coeff( :, i)) == 0)
        coeff(k,i) = 1;
        b(k) = 0;
        k=k+1;
    end
end


weights = coeff \ b;


if (opt == 2)
    weights = reshape(weights, nside, nside);
end

fid = fopen('preisach.txt', 'w+');

if (opt == 1)

    fprintf( fid, '%i\n', length(weights) );
    for i = 1:length(weights)
        fprintf(fid, '%e %e %e\n', alpha(i) * 1e-9, beta(i)* 1e-9, weights(i)* 1e-9);
    end
else
    fprintf( fid, '%i\n', length(weights) );
    for i = 1:length(weights)
        fprintf(fid, '%e %e\n', alpha(i) * 1e-9, beta(i)* 1e-9);
    end
    
    for i = 1:length(weights)
        for j = 1:i
            fprintf(fid, '%e\n', weights(i,j) * 1e-9);
        end
    end
            
end
fclose(fid);

%%
% %this section runs an example to check the weights
if (~ isempty(do_plots))
    clear dt f
    
    t = linspace(0, 1, 2*nside);
    
    dt =  1.2-1.19*sin( pi * t)    ;
    
    f = zeros(1, length(dt) );
    
    if (opt == 1)
        states = zeros(size(alpha));
        f(1) = sum(sum(weights' .*  states));
    else
        states = zeros(nside, nside);
        f(1) = sum(sum(weights .*  states));
    end        
    
    for i = 2:length(dt)
        
        states = update_states(dt(i), dt(i-1), states, alpha, beta);
        if (opt == 1)
            f(i) = sum(sum(weights' .*  states));
        else
            f(i) = sum(sum(weights .*  states));
        end
    end
    
    figure
    plot( dt, f, 'LineWidth', 2);
    xlabel('gap')
    ylabel('force')
    
    
    figure; plot( t, dt, t, f, 'LineWidth', 2)
    legend('gap', 'force');
    xlabel('time')
end

if (isdeployed)
    quit
end

function states = update_states(d, old_d, old_states, alpha, beta)


states = old_states;

if (size( states, 1) > 1)
    %states is a 2-d matrix, alpha is a 1-d vectors    
    
    if ( d < old_d )
        %advancing                
       ndx = find( d < alpha);
        for i = ndx
           states(i, 1:i) = ones(1, i);
        end
        
%         i = find( d < alpha, 1, 'last');
%         if (i == 1)
%             states(1,1) = 1;
%         else
%             states(1:i, 1:i) = tril( ones(i) );
%         end
%          
    else
        %retracting
        ndx = find( d > beta);
        states(:,ndx) = zeros(size(states,1),  length(ndx));
    end
else
    %states is a 1-d vector

    if ( d < old_d )
        %gap decreasing
        ndx = find( d < alpha);
        states(ndx) = ones(length(ndx),1);
    else
        %gap increasing
        ndx = find( d > beta);
        states(ndx) = zeros(length(ndx),1);
    end
end
%function rappture_reader_multi(prefix, filendx, suffix)
%this file reads multiple rappture runs into matlab at one time.
%they need to have a common prefix and suffix and be numbered in order
%for example run1.xml, run2.xml, run3.xml, etc.
%then call rappture_reader_multi('run', 1:3, 'xml')
function rappture_reader_multi(prefix, filendx, suffix, curve_suffix, ignore_impact)

if (nargin < 5)
    ignore_impact = 0;
end

if (nargin < 4)
    curve_suffix = '';
end

if (nargin < 3)
    suffix = 'xml';
end


system(['gzip -d ' prefix '*.gz']);
for ndx = filendx
    filename = [prefix num2str(ndx) suffix]

    % for very large exponents, such as 1.234E-200, fortran likes to drop
    % the E to save one character, so it gets 1.2345-200.  matlab doesn't
    % understand this, and I don't know how to make fortran stop this, so we have this hack
   %windows machines will not typically have sed installed. If you get errors,
   %you can comment this line out, (but if you have such large exponents in your
   %file you'll have to edit them manually)

    system(['sed -e ''s/\(\.[0-9]*\)[0-9]\([-+]\)\([0-9]*\)/\1E\2\3/g'' ' filename ' > tmp.out']);
    
    xdoc = xmlread('tmp.out');
    
    numCurves = xdoc.getElementsByTagName('curve').getLength;
    for i = 0:numCurves-1
        curve = xdoc.getElementsByTagName('curve').item(i);
        curveName = char(curve.getAttribute('id'));
               
        if ( (ignore_impact==1) && ~isempty(strfind( curveName, '_Impact')))
            'ignored impact'
        else        
            clear curvetxt
            comp = curve.getElementsByTagName('component');
            if (comp.getLength() > 0)                                
                curvetxt = char(comp.item(0).getElementsByTagName('xy').item(0).getTextContent);
                [x,y] = strread( curvetxt, '%f %f');
                
                s = [ curveName '_x(1:length(x),ndx) = x;']; eval(s);
                s = [ curveName '_y(1:length(y),ndx) = y;']; eval(s);
        
                %this way is more elegant, but it forces all curves to have same length, which causes us to
                %run out of memory when there are time history or poincare curves that are very long.
                %xdata(i+1,1:length(x), ndx) = x;
                %ydata(i+1,1:length(x), ndx) = y;
        
                %pad end with NaNs. solves problems when some vectors are longer than others,
                %xdata( i+1, length(x)+1, ndx) = NaN;
                %ydata( i+1, length(x)+1, ndx) = NaN;
                s = [ curveName '_x(length(x)+1,ndx) = NaN;']; eval(s);
                s = [ curveName '_y(length(x)+1,ndx) = NaN;']; eval(s);
            end
        end
    end
    
    numStrings = xdoc.getElementsByTagName('string').getLength;
    for i = 0:numStrings-1
        curve = xdoc.getElementsByTagName('string').item(i);
        nod = curve.getElementsByTagName('current').item(0);
        if (~ isempty(nod))
            c{i+1,ndx} = char(nod.getTextContent);
        end
    end
    
    numNums = xdoc.getElementsByTagName('number').getLength;
    for i = 0:numNums-1
        curve = xdoc.getElementsByTagName('number').item(i);
        nod = curve.getElementsByTagName('current').item(0);
        if (~ isempty(nod))
            n{i+1,ndx} = char(nod.getTextContent);
        end
    end
end

for ndx = filendx
    ss = who('*_x');
    for i = 1:size(ss,1)
        s = [ 'nn = find( isnan( ' ss{i} ' ( :, ndx)), 1); ']; 
        try 
            eval(s);
        catch
            continue;
        end
 
        if ( ~ isempty(nn))
            prefix = ss{i}(1:end-2);
            
            s =  [ prefix, '_x(nn:end, ndx) = NaN;']; eval(s);
            s =  [ prefix, '_y(nn:end, ndx) = NaN;']; eval(s);
        end
%        if (~ isempty( nn) )
%            xdata(i+1, nn:end, ndx) = NaN;
%            ydata(i+1, nn:end, ndx) = NaN;
%        end
    end
end

% for i = 0:numCurves-1
%     curve = xdoc.getElementsByTagName('curve').item(i);
%     curveName = char(curve.getAttribute('id'));
%     assignin('base', [curveName '_x' curve_suffix], squeeze(xdata(i+1,:,:)));
%     assignin('base', [curveName '_y' curve_suffix], squeeze(ydata(i+1,:,:)));       
% end

ss = who('*_x');
for i = 1:size(ss,1)
    curveName = ss{i}(1:end-2);
    assignin( 'base', [curveName '_x' curve_suffix], eval([curveName '_x']));
    assignin( 'base', [curveName '_y' curve_suffix], eval([curveName '_y']));
end

for i = 0:numStrings-1
     curve = xdoc.getElementsByTagName('string').item(i);
     curveName = char(curve.getAttribute('id'));
     assignin('base', [curveName curve_suffix], squeeze(c(i+1,:)));
end

for i = 0:numNums-1
     curve = xdoc.getElementsByTagName('number').item(i);
     curveName = char(curve.getAttribute('id'));
     assignin('base', [curveName curve_suffix], squeeze(n(i+1,:)));
end
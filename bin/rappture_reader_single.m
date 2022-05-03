%function rappture_reader_single( filename, suffix) 
%this file reads a single rappture file into matlab.  
%suffix is optional.
% example rappture_reader_single('run123456.xml')

function rappture_reader_single( filename, suffix) 

if (nargin < 2) 
    suffix = [];
end

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
        if (~ isempty( curve.getElementsByTagName('component').item(0) ))
            curvetxt = char(curve.getElementsByTagName('component').item(0).getElementsByTagName('xy').item(0).getTextContent);
            [x,y] = strread( curvetxt, '%f %f');

            assignin( 'caller',  [curveName suffix '_x' ], x);
            assignin( 'caller',  [curveName suffix '_y'], y);
        end
    end

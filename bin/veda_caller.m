%if want_offset is true, them param_vals(1) is a Z-offset to match up the curves with
%
% example:
%results =         veda_caller( [1e-5],    curve_x_data, {'TAG_C'} ,   'Amp_x', {'Phase_y'}, 'driver1.xml' ,    false,        [],              [],          ydata);
%ydata is just for plotting
function results = veda_caller( param_vals, curve_x_data, param_names, x_name, result_names, template_filename, want_offset, param_scalars, result_scalars, ydata)

if (nargin < 8)
    param_scalars = ones(1, length(param_vals)+1);
end

if (isempty(param_scalars))
    param_scalars = ones(1, length(param_vals)+1);
end



 % substitute parameters values into a text file and then run veda
 s = [];

 
 fid=fopen('runs.txt', 'w+');
 
 if (want_offset)
     disp(param_vals .* param_scalars(2:end))
     poff = 1;
 else
     disp(param_vals .* param_scalars)
     poff = 0;
 end
 
 for i = 1:length(param_vals)     
     fprintf(fid, '%s = %16.15e , ', param_names{i}, param_vals(i+poff) * param_scalars(i+poff) );
 end
 
 fclose(fid);
 
 system(['perl ~/app-adac-dev/bin/param-explore.pl runs.txt ' template_filename ' veda_wrapper.pl ']);
 

 % read results
rappture_reader_single('1.out');
 
%figure(1);
set(0, 'CurrentFigure', 1);

 % collect results into a vector, interpolating each onto the curve_x_data
 results = [];
 N = length( result_names);
 
 if (want_offset)
         xx = curve_x_data + param_vals(1);
 else
         xx = curve_x_data;
 end
 
 for i = 1:N
     %extrapolating usually bites us.  just make the simulation Z distance longer
     s = ['rr = interp1(' x_name ', (' result_names{i} ')*result_scalars(i), xx ,''linear'' );' ];
     eval(s);
     
     results = [results rr];
               
     subplot(N,1,i)
     plot( xx, rr, xx, ydata( (i-1) * length(rr) + (1:length(rr)) ))
 end

     
 
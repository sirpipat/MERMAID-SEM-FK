function writetimeseries(t, x, fname)
% WRITETIMESERIES(t, x, fname)
%
% Write time series data x = x(t) to a file
%
% INPUT:
% t         time
% x         data
% fname     file name
%
% EXAMPLE:
% % real time series data written to 'timeseries_demo_real.txt'
% % the line format is '%g %g\n'
% writetimeseries('demo1');
%
% % complex time series data written to 'timeseries_demo_complex.txt'
% % the line format is '%g (%g,%g)\n'
% writetimeseries('demo2');
%
% SEE ALSO:
% READTIMESERIES
%
% Last modified by sirawich-at-princeton.edu, 02/17/2025

% demos set up
if (ischar(t) || isstring(t)) && strcmp(t, 'demo1')
    t = (0:10)';
    x = 5 - t;
    fname = 'timeseries_demo_real.txt';
elseif (ischar(t) || isstring(t)) && strcmp(t, 'demo2')
    t = (0:10)';
    x = sqrt(5 - t);
    fname = 'timeseries_demo_complex.txt';
end

% make sure that t and x are column vectors
if size(t, 1) == 1
    t = t';
end
if size(x, 1) == 1
    x = x';
end

% default file name that the time series data are written to
DEFAULT_FNAME = 'timeseries.txt';
defval('fname', DEFAULT_FNAME)

fid = fopen(fname, 'w');
if strcmp(fname, DEFAULT_FNAME)
    fprintf('timeseries is being written to "timeseries.txt"\n');
end

% write the time series data
if isreal(x)
    fprintf(fid, '%g %g\n', [t x]');
else
    fprintf(fid, '%g (%g,%g)\n', [t real(x) imag(x)]');
end

fclose(fid);
end
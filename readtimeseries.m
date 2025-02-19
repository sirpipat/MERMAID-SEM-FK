function [t, x] = readtimeseries(fname, iscomplex)
% [t, x] = READTIMESERIES(fname, iscomplex)
%
% Reads a time series data set from a file as x = x(t).
%
% INPUT:
% fname         name of the file
% iscomplex     whether the data are complex numbers [default: false]
%       - if false, line format is '%g %g\n'
%       - if true,  line format is '%g (%g,%g)\n'
%
% OUTPUT:
% t             time
% x             data
%
% EXAMPLE:
% % read real time series data
% writetimeseries('demo1');
% [t, x] = readtimeseriesdata('timeseries_demo_real.txt', false);
%
% % read complex time seriesd data
% writetimeseries('demo2');
% [t, x] = readtimeseriesdata('timeseries_demo_complex.txt', true);
%
% SEE ALSO:
% WRITETIMESERIES
%
% Last modified by sirawich-at-princeton.edu, 02/17/2025

defval('iscomplex', false)

fid = fopen(fname, 'r');
if iscomplex
    data = fscanf(fid, '%g (%g,%g)\n');
    data = reshape(data, [3 length(data)/3])';
    t = data(:,1);
    x = data(:,2) + 1i * data(:,3);
else
    data = fscanf(fid, '%g %g\n');
    data = reshape(data, [2 length(data)/2])';
    t = data(:,1);
    x = data(:,2);
end
fclose(fid);
end
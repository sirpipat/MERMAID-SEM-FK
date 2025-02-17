function [t, vx, vy, vz, hdr] = readplotFKfile(fname)
% [t, vx, vy, vz] = READPLOTFKFILE(fname)
%
% Reads a plot FK file and formats to easy to use data. Output values are
% time and either velocity or traction depending on the file.
%
% INPUT:
% fname         name to a plot FK file
%
% OUTPUT:
% t             time
% vx            value in x direction
% vy            value in y direction
% vz            value in z direction
% hdr           file header
%
% Last modified by sirawich-at-princeton.edu, 02/17/2025

% open file
fid = fopen(fname, 'r');
% read header
hdr = '';
for ii = 1:4
    hdr = [hdr fgets(fid)];
end
% read data
data = fscanf(fid, '%g');
data = reshape(data, [4 length(data)/4])';

% format to output
t = data(:,1);
vx = data(:,2);
vy = data(:,3);
vz = data(:,4);
end

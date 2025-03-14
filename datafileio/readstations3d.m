function [n, name, network, x, y, z] = readstations3d(fname)
% [n, name, network, x, y, z] = READSTATIONS3D(fname)
% 
% Reads a STATION file for SPECFEM3D simulation.
%
% INPUT
% fname         full filename of a STATION file
%
% OUTPUT
% n             the number of stations
% name          station names
% network       station network names
% x             x-coordinates
% z             z-coordinates
%
% SEE ALSO:
% WRITESTATIONS3D, READ_STATION
%
% Last modified by Sirawich Pipatprathanporn, 03/14/2025

% read the station file as a table
opts = detectImportOptions(fname, 'FileType', 'text');
T = readtable(fname, opts);

% if the table is empty, assume the table has only one entry
if isempty(T)
    fid = fopen(fname, 'r');
    line = fgetl(fid);
    words = split(line);
    T = struct('Var1', {words(1)}, 'Var2', {words(2)}, ...
        'Var3', str2double(words{3}), ...
        'Var4', str2double(words{4}), ...
        'Var5', str2double(words{5}));
end

name = T.Var1;
network = T.Var2;
x = T.Var3;
y = T.Var4;
z = T.Var5;
n = size(T, 1);
end
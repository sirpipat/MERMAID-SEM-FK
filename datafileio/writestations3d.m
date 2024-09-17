function writestations3d(stations, fname)
% WRITESTATIONS3D(stations, fname)
%
% Writes a STATION file of SPECFEM3D_Cartesian.
% 
% INPUT:
% stations          struct containing following fields
%       name            station name
%       network         network name
%       lat             latitude  or y-coordinate
%       lon             longitude or x-coordinate
%       elev            elevation
%       z               burial or z-coorndiate or depth 
% fname             filename of a STATION file
%
% Last modified by sirawich-at-princeton.edu, 09/17/2024

if isempty(fname)
    % standard output aka console output
    fid = 1;
else
    fid = fopen(fname, 'w');
end

for ii = 1:length(stations.name)
    fprintf(fid, '%5s %2s %11.4f %11.4f %11.4f %11.4f\n', ...
        stations.name{ii}, stations.network{ii}, stations.lat(ii), ...
        stations.lon(ii), stations.elev(ii), stations.z(ii));
end

% close the file
if fid >= 3
    fclose(fid);
end
end
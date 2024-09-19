function [itfs, layers] = loadinterfacefiles3d(fname)
% [itfs, layers] = LOADINTERFACEFILES3D(fname)
%
% Reads interfaces from an interface file and other accompanied files for a
% SPECFEM3D_Cartesian simulation. Please make sure that other accompanied
% files are in the same directory as the interface file.
%
% INPUT:
% fname         name of the interface file
%
% OUTPUT:
% itfs          interfaces, an array of struct with following fields
%       SUPPRESS_UTM_PROJECTION     whether to suppress UTM projection
%       NXI                         number of elements in x-direction
%       NETA                        number of elements in y-direction
%       LON_MIN                     minimum longitude (or x value)
%       LAT_MIN                     minimum latitude  (or y value)
%       SPACING_XI                  spacing in x-direction
%       SPACING_ETA                 spacing in y-direction
%       FILE                        elevation file name
%       Z                           elevation grid at (X,Y) or (LON,LAT)
% layers        number of vertical spectral elements for each layer
%
% SEE ALSO:
% MAKEINTERFACES3D, WRITEINTERFACEFILES3D
%
% Last modified by sirawich-at-princeton.edu, 09/19/2024

ddir = strcat(fileparts(fname), filesep);

%% open the file
fid = fopen(fname, 'r');
line = strip(fgetl(fid));
% skip comments / headers
while isempty(line) || strcmp(line(1), '#')
    line = strip(fgetl(fid));
end

%% read number of interfaces
numinterfaces = sscanf(line, '%d', 1);

itfs = cell(numinterfaces, 1);
layers = nan(numinterfaces, 1);

%% read the interfaces
for ii = 1:numinterfaces
    line = strip(fgetl(fid));
    % skip comments / headers
    while isempty(line) || strcmp(line(1), '#')
        line = strip(fgetl(fid));
    end
    [boolword, ~, ~, ni] = sscanf(line, '%s', 1);
    if strcmp(boolword, '.true.')
        itfs{ii}.SUPPRESS_UTM_PROJECTION = true;
    elseif strcmp(boolword, '.false.')
        itfs{ii}.SUPPRESS_UTM_PROJECTION = false;
    else
        warning(['SUPPRESS_UTM_PROJECTION flag could not be ' ...
            'interpreted, assumed to be true then.'])
        itfs{ii}.SUPPRESS_UTM_PROJECTION = true;
    end
    nums = sscanf(line(ni:end), '%f');
    itfs{ii}.NXI = nums(1);
    itfs{ii}.NETA = nums(2);
    itfs{ii}.LON_MIN = nums(3);
    itfs{ii}.LAT_MIN = nums(4);
    itfs{ii}.SPACING_XI = nums(5);
    itfs{ii}.SPACING_ETA = nums(6);
    
    line = strip(fgetl(fid));
    itfs{ii}.FILE = strcat(ddir, sscanf(line, '%s', 1));
    % try to read the elevation from the file
    try
        fid_ii = fopen(itfs{ii}.FILE, 'r');
        z = fscanf(fid_ii, '%f');
        itfs{ii}.Z = reshape(z, itfs{ii}.NXI, itfs{ii}.NETA)';
        fclose(fid_ii);
    catch ME
        fprintf('Encounter problems while reading the file %s\n', ...
            itfs{ii}.FILE);
        getReport(ME)
        continue
    end
end

%% read the number of spectral elements in the vertical direction
for ii = 1:numinterfaces
    line = strip(fgetl(fid));
    % skip comments / headers
    while isempty(line) || strcmp(line(1), '#')
        line = strip(fgetl(fid));
    end
    layers(ii) = sscanf(line, '%d', 1);
end

%% close the file
fclose(fid);
end
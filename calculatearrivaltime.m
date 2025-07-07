function t0 = calculatearrivaltime(ddir)
% t0 = CALCULATEARRIVALTIME(ddir)
%
% Calculates the arrival time at the ocean bottom station in MERMAID 
% FK-SPECFEM3D runs given the origin time and origin wavefront.
%
% INPUT:
% ddir          directory to a FK-SPECFEM3D run
%
% OUTPUT:
% t0            rrival time at the ocean bottom station in FK-SPECFEM3D run
%
% Last modified by sirwich-at-princeton.edu, 07/07/2025

% determine the location of wave front at origin time
fkmodel = loadfkmodel(fullfile(ddir, 'DATA', 'FKMODEL'));
t0 = fkmodel.origin_time;
x0 = fkmodel.origin_wavefront(2);
z0 = fkmodel.origin_wavefront(3);
theta = fkmodel.theta * pi/180; % in radians

% determine the location of the ocean bottom right below the float
[~, ~, ~, ~, y, z] = readstations3d(fullfile(ddir, 'DATA', 'STATIONS'));
x1 = y(1);
z1 = z(1);

if fkmodel.nlayers < 1
    error('No layers found in the model.');
elseif fkmodel.nlayers == 1
    vp = fkmodel.layers{1}.vp;
    t0 = t0 + ((x1 - x0) * sin(theta) + (z1 - z0) * cos(theta)) / vp;
elseif fkmodel.nlayers == 2
    vp = fkmodel.layers{2}.vp;
    t0 = t0 + ((x1 - x0) * sin(theta) + (z1 - z0) * cos(theta)) / vp;
else
    % Calculate arrival time for multiple layers
    vp = fkmodel.layers{fkmodel.nlayers}.vp;
    for ii = fkmodel.nlayers:-1:3
        theta = asin(fkmodel.layers{ii}.vp * sin(theta) / vp);
        vp = fkmodel.layers{ii}.vp;
        t0 = t0 + (fkmodel.layers{ii}.ztop - z0) * cos(theta) / vp;
        x0 = x0 + (fkmodel.layers{ii}.ztop - z0) * tan(theta);
        z0 = fkmodel.layers{ii}.ztop;
    end
    % Add the travel time for the last layer
    theta = asin(fkmodel.layers{2}.vp * sin(theta) / vp);
    vp = fkmodel.layers{2}.vp;
    t0 = t0 + ((x1 - x0) * sin(theta) + (z1 - z0) * cos(theta)) / vp;
end
end
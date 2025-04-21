function seis = fluidsolidseismogram(z_interface, z_station, rho_f, vp_f, rho_s, vp_s, vs_s, theta, s)
% seis = FLUIDSOLIDSEISMOTRAM(z_interface, z_station, rho_f, vp_f, ...
%     rho_s, vp_s, vs_s, theta, stf)
%
% Generate seismogram of the waves arriving at
% a stations inside a fluid layer on top of a solid half-space with the
% P-wave entering the system from the bottom.
%
% INPUT:
% z_interface   depth of the interface (given the top of fluid layer is 0)
% z_station     depth of the station
% rho_f         density of the fluid
% vp_f          P-wave velocity in the fluid
% rho_s         density of the solid
% vp_s          P-wave velocity in the solid
% vs_s          S-wave velocity in the solid
% theta         incidence angle of the P-wave entering from the bottom
%               (in radians)
% s             time-series [t; d] of an incoming P-wave where d is the
%               displacement on the direction of wave propagation
%
% OUTPUT:
% seis          time-series [t; x; z] of the incoming wave and reflected
%               waves in both x and z directions
%
% Last modified by sirawich-at-princeton.edu: 04/01/2025

t = s(:,1);
x = s(:,2) * sin(theta);
z = s(:,2) * cos(theta);

% determine max time
t_max = s(end,1) - s(1,1);

% compute arrival times and relative amplitudes
[ta, xa, za, pha] = fluidsolidarrivals(z_interface, z_station, rho_f, ...
    vp_f, rho_s, vp_s, vs_s, theta, t_max);

% construct the response functions
xw = placeweight(ta, xa, t - t(1));
zw = placeweight(ta, za, t - t(1));

x = indeks(conv(x, xw), 1:length(x));
z = indeks(conv(z, zw), 1:length(z));

seis = [t, x, z];
end
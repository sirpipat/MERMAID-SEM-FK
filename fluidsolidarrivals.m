function [t, x, z] = fluidsolidarrivals(z_interface, z_station, rho_f, vp_f, rho_s, vp_s, vs_s, theta, t_max)
% [t, x, z] = FLUIDSOLIDARRIVALS(z_interface, z_station, rho_f, vp_f, ...
%     rho_s, vp_s, vs_s, theta, t_max)
%
% Compute arrival times and displacement amplitude of the waves arriving at
% a stations inside a fluid layer on top of a solid half-space with the
% P-wave entering the system from the bottom. The first P-wave arrival is
% zero and its displacement amplitudes are one.
%
% INPUT:
% z_interface   depth of the interface (given the top of fluid layer is 0)
% z_station     depth of the station
% rho_f         density of the fluid
% vp_f          P-wave velocity in the fluid
% rho_s         density of the solid
% vp_s          P-wave velocity in the solid
% vs_s          S-wave velocity in the solid
% theta         incident angle of the P-wave entering from the bottom
%               (in radians)
% t_max         maximum time for consideration
%
% OUTPUT:
% t             arrival time
% x             displacement amplitude in the horizontal direction
% z             displacement amplitude in the vertical direction
%
% Last modified by sirawich-at-princeton.edu 03/05/2025

%% demos
if ischar(z_interface) && strcmp(z_interface, 'demo')
    return
end

%% angles of incident for each wave
theta_s = asin(vs_s / vp_s * sin(theta));
theta_f = asin(vp_f / vp_s * sin(theta));

%% displacement amplitude coefficients
% solid to fluid
[~, D_S2F] = fluidsolidcoefficients2(rho_f, vp_f, rho_s, vp_s, vs_s, theta);
D_S2F_RP   = D_S2F(1);   % reflected P
D_S2F_RSV  = D_S2F(2);   % reflected SV
D_S2F_T    = D_S2F(3);   % transimitted acoustic

% free-surface reflection in an acoustic media
D_FREE_R   = -1;

% fluid to solid
[~, D_F2S] = fluidsolidcoefficients(rho_f, vp_f, rho_s, vp_s, vs_s, theta_f);
D_F2S_R    = D_F2S(1);   % reflected acoustic
D_F2S_TP   = D_F2S(2);   % transmitted P
D_F2S_TSV  = D_F2S(3);   % transmitted SV

%% sign to convert down-going wave to displacement
SGN_X      =  1;
SGN_Z      = -1;

%% determine which layer the station is in
if (z_interface > z_station)
    layer = "fluid";
else
    layer = "solid";
end

t = 0;
x = 1;
z = 1;
if strcmp(layer, "solid")
    % travel time of first down-going P-wave
    tPP = 2 * (z_station - z_interface) / (vp_s * cos(theta));
    
    % travel time of first down-going SV-wave
    tPSV = (z_station - z_interface) * (1 / (vp_s * cos(theta)) + ...
        1 / (vs_s * cos(theta_s)));
    
    % round-trip travel time of acoustic wave in the fluid layer
    tFF = 2 * z_interface / (vp_f * cos(theta_f));
    
    % number of rounds
    nP = floor((t_max - tPP) / tFF);
    nSV = floor((t_max - tPSV) / tFF);
    
    % down-going P-wave
    % [P; P^bathP; P(A^A)P]     % A denote acoustic phase
    t = [t; tPP + (0:nP)' * tFF];
    x = [x; [D_S2F_RP; D_S2F_T * (D_FREE_R .^ (1:nP)') .* ...
        (D_F2S_R .^ (0:nP-1)') * D_F2S_TP] * SGN_X];
    z = [z; [D_S2F_RP; D_S2F_T * (D_FREE_R .^ (1:nP)') .* ...
        (D_F2S_R .^ (0:nP-1)') * D_F2S_TP] * SGN_Z];
    
    % down-going SV-wave
    % [P; P^bathS; P(A^A)S]     % A denote acoustic phase
    t = [t; tPSV + (0:nSV)' * tFF];
    x = [x; [D_S2F_RSV; D_S2F_T * (D_FREE_R .^ (1:nP)') .* ...
        (D_F2S_R .^ (0:nP-1)') * D_F2S_TSV] * SGN_X * ...
        cos(theta_s) / sin(theta)];
    z = [z; [D_S2F_RSV; D_S2F_T * (D_FREE_R .^ (1:nP)') .* ...
        (D_F2S_R .^ (0:nP-1)') * D_F2S_TSV] * SGN_Z * ...
        sin(theta_s) / cos(theta)];
else
    % travel time of the first down-going acoustic wave
    tAA = 2 * z_station / (vp_f * cos(theta_f));
    
    % round-trip travel time of acoustic wave in the fluid layer
    tFF = 2 * z_interface / (vp_f * cos(theta_f));
    
    % number of rounds
    n_up = floor(t_max / tFF);
    n_down = floor((t_max - tAA) / tFF);
    
    % up-going acoustic wave arrivals
    t = [t; (1:n_up)' * tFF];
    x = [x; (D_FREE_R * D_F2S_R) .^ (1:n_up)'];
    z = [z; (D_FREE_R * D_F2S_R) .^ (1:n_up)'];
    
    % down-going acoustic wave arrivals
    t = [t; tAA + (0:n_down)' * tFF];
    x = [x; D_FREE_R * (D_F2S_R * D_FREE_R) .^ (0:n_down)' * SGN_X];
    z = [z; D_FREE_R * (D_F2S_R * D_FREE_R) .^ (0:n_down)' * SGN_Z];
end

% sort arrival by time
[t, it] = sort(t);
x = x(it);
z = z(it);

% vertical incident does not have any horizontal displacement
if abs(theta) == 0
    x = 0 * t;
end
end
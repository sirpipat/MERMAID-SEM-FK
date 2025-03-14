function y = placeweight(loc, w, x)
% y = PLACEWEIGHT(loc, w x)
%
% Places weight on grid points with equally spacing. When the weight is
% not on any grid point, it will be distributed to two adjacent grid points
% with the ratio being the inverse distance from the weight.
%
% INPUT:
% loc           locations of the weights
% w             weights
% x             grid points with equally spacing         
%
% OUTPUT:
% y             weights at grid points
%
% EXAMPLE:
% % basic function call
% locs = sqrt(5) + (0:4)' * 2*pi;
% weights = [1 -0.8 0.64 -0.512 0.4096]';
% x = (0:0.1:40)';
% y = placeweight(locs, weights, x);
%
% % Useful for convolution of delta functions with off-grid locations
% t = (-10:0.1:40)';
% s = exp(-t.^2 / 0.25^2);
% u = indeks(conv(s, y), 1:length(s));
%
% Last modified by sirawich-at-princeton.edu: 03/14/2025

% place each weight to the nearest grid point(s)
y = zeros(size(x));
dx = (x(end) - x(1)) / (length(x) - 1);
for ii = 1:length(loc)
    % Find the distance scaled by sampling time from a stem point to each 
    % grid point
    dist = abs(x - loc(ii)) / dx;
    % Weight is maximum of 1 when the stem point is right at a grid point
    % and decrease to zero when the scaled distance is greater than 1
    weight = max(1 - dist, 0);
    y = y + weight * w(ii);
end
end
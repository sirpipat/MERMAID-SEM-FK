function [A, B, C, ZZ] = detrend2d(zz, xx, yy)
% [A, B, C, ZZ] = detrend2d(zz, xx, yy)
%
% Detrends the surface zz given a mesh [xx, yy]. It first use least-squared
% method to fits the plane equation zz = A + B * xx + C * yy. Then it 
% computes the residue after the plane is removed:
%
% ZZ = zz - (A + B * xx + C * yy)
%
% INPUTS:
% zz            2D grid
% xx            meshgrid for x-value [Default: (1:size(zz,2))']
% yy            meshgrid for y-value [Default: (1:size(zz,1))']
%
% OUTPUTS:
% A, B, C       best-fitted parameters of a plane z = A + Bx + Cy
% ZZ            residue once linear trend is removed
%
% Example
% x = (-1000:25:1000)';
% y = (-2000:25:2000)';
% [xx, yy] = meshgrid(x, y);
% zz = 1200 + 0.4 * xx - 0.25 * yy + 50 * randn(size(xx));
% 
% % Quick detrend
% [~, ~, ~, zz2] = detrend2d(zz);
%
% % Estimate slope
% [A, B, C] = detrend2d(zz, xx, yy);
%
% Last modified by sirawich-at-princeton.edu, 06/12/2025

defval('xx', (1:size(zz,2))')
defval('yy', (1:size(zz,1))')

if size(xx, 2) == 1 && size(yy, 2) == 1
    [xx, yy] = meshgrid(xx, yy);
end

% Least-squared fitting to a plane
G = [ones(numel(xx),1) xx(:) yy(:)];
m = (G' * G) \ (G' * zz(:));

A = m(1);
B = m(2);
C = m(3);

ZZ = zz - (A + B * xx + C * yy);
end
function [DD, p, zz_gaussian] = degaussian2d(zz, xx, yy)
% [zz_residue] = DEGAUSSIAN2D(zz, xx, yy)
%
% Removes the best fit Gaussian 2D from 2D data and returns the remaining
% data.
%
% INPUT:
% zz            2D data
% xx            meshgrid x-coordinates [optional, default: (1:size(zz, 2))]
% yy            meshgrid y-coordinates [optional, default: (1:size(zz, 1))]
%
% OUTPUT:
% DD            2D data after the best fit Gaussian 2D is removed
% p             best-fitted parameters for Gaussian 2D curve
%               z = p(1) + p(2) * exp(-((x-p(3))^2 / (2*p(5)^2) + (y-p(4))^2 / (2*p(6)^2)))
% zz_gaussian   the Gaussian 2D hill that is removed
% 
% Last modified by sirawich-at-princeton.edu, 06/25/2025

defval('xx', (1:size(zz, 2)))
defval('yy', (1:size(zz, 1)))

if size(xx, 1) == 1
    [xx, yy] = meshgrid(xx, yy);
end

xy = [xx(:) yy(:)];
z = zz(:);

% Initial parameter guess for the Gaussian fit
p0 = [min(z), range(z), mean(xx(:)), mean(yy(:)), std(xx(:)), std(yy(:))];

% Perform the fitting using nonlinear least squares
options = optimset('Display', 'off');
p = lsqcurvefit(@gaussian2d, p0, xy, z, [], [], options);

zz_gaussian = reshape(gaussian2d(p, xy), size(yy));
DD = zz - zz_gaussian;
end

function z = gaussian2d(p, xy)
z = p(1) + p(2) * exp(-((xy(:,1) - p(3)).^2 / (2 * p(5)^2) + ...
    (xy(:,2) - p(4)).^2 / (2 * p(6)^2)));
end
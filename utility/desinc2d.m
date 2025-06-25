function [DD, p, zz_sinc] = desinc2d(zz, xx, yy)
% [DD, p, zz_sinc2d] = DESINC2D(zz, xx, yy)
%
% Removes the best fit sinc 2D from 2D data and returns the remaining
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
%               z = p(1) + p(2) * sinc(pi * ((x-p(3))^2 / (2*p(5)^2) + (y-p(4))^2 / (2*p(6)^2)))
% zz_sinc       the sinc 2D hill that is removed
% 
% Last modified by sirawich-at-princeton.edu, 06/25/2025
%
% EXAMPLE:
% % quick function call to remove sinc function
% DD = desinc2d(zz);
%
% % remove sinc function and estimate sinc parameters for z(x,y)
% [DD, p, zz_sinc] = desinc2d(zz, xx, yy);
%
% SEE ALSO:
% DEGAUSSIAN2D
%
% Last modified by sirawich-at-princeton.edu, 06/25/2025

defval('xx', (1:size(zz, 2)))
defval('yy', (1:size(zz, 1)))

if size(xx, 1) == 1
    [xx, yy] = meshgrid(xx, yy);
end

% flatten the data (z) and observation (x,y)
xy = [xx(:) yy(:)];
z = zz(:);

% Initial parameter guess for the Gaussian fit
p0 = [min(z), range(z), indeks(xy(z == max(z), 1), 1), ...
    indeks(xy(z == max(z), 2), 1), std(xx(:)), std(yy(:))];

% Perform the fitting using nonlinear least squares
options = optimset('Display', 'off');
p = lsqcurvefit(@sinc2d, p0, xy, z, [], [], options);

zz_sinc = reshape(sinc2d(p, xy), size(yy));
DD = zz - zz_sinc;
end

function z = sinc2d(p, xy)
z = p(1) + p(2) * sinc((xy(:,1) - p(3)).^2 / p(5)^2 + ...
    (xy(:,2) - p(4)).^2 / p(6)^2);
end
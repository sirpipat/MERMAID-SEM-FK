function makegenericoceanbottom(boxsize, elemsize, x0, y0, z0, rot)
% MAKEGENERICOCEANBOTTOM(boxsize, elemsize, x0, y0, z0, rot)
%
% Makes generic ocean bottom topography files for SPECFEM3D runs
%
% ---- ISOTROPIC ----
% 1. Gaussian hill
% 2. Gaussian basin
% 3. Gaussian crator -- hill with hole on top, (A * exp(-r^2) * cosh(r))
% 4. Hill that bends the ground -- sinc(r)
% --- ANISOTROPIC ---
% 5. Elliptic hill -- exp(-(x/a)^2 - (y/b)^2)
% 6. ridge -- 1D Gaussian where the ridge is along y-axis when rot=0
% 7. trench -- upside-down ridge
% 8. rolling ridges -- 1D sine wave where the ridges is along y-axis
% 9. Slope -- inclined plane
%
% INPUT:
% boxsize           width of the square box
% elemsize          element size, a square element
% x0, y0            center location of the feature
% z0                height of such the feature
% rot               rotation angle (in degrees) clockwise from north
%                   (equivalent to azimuth angle) for anisotropic feature
% 
% SEE ALSO:
% MAKESTOCKFILES3D, COOKSPECFEM3DRUN
%
% Last modified by sirawich-at-princeton.edu, 06/04/2025

defval('x0', 0)
defval('y0', 0)
defval('z0', 0)
defval('rot', 60)

rot = rot * pi/180;

stockdir = fullfile(getenv('REMOTE3D'), 'stockfiles');
system(sprintf('mkdir -p %s', stockdir));

topodir = fullfile(stockdir, 'topofiles2');
system(sprintf('mkdir -p %s', topodir));

% mesh
n = round(boxsize / elemsize) + 1;
x = linspace(-1, 1, n)';
y = x;
[xx, yy] = meshgrid(x, y);
rr = sqrt((xx-x0/(boxsize/2)).^2 + (yy-y0/(boxsize/2)).^2);

% rotated clockwise
xx_rot = cos(rot) * (xx - x0/(boxsize/2)) - sin(rot) * (yy - y0/(boxsize/2));
yy_rot = sin(rot) * (xx - x0/(boxsize/2)) + cos(rot) * (yy - y0/(boxsize/2));

% tapering window
w = shanning(n, 0.15);
wn = w' .* w;

% circular tapering window
rr0 = sqrt(xx.^2 + yy.^2);
w1 = [w(floor(n/2):end); 0];
x1 = [x(floor(n/2):end); 3];
wr = interp1(x1, w1, rr0);

% hill (Gaussian)
% circular
zz = z0 * exp(-(5*rr).^2 / 2);

figure(1)
clf
surface(x, y, wr.*zz);
colorbar
set(gca, 'DataAspectRatio', [1 1 z0]);
colormap(kelicol)

fname = sprintf('topo_GENERIC_HILL_BOX-%d_ELEM-%d_X-%d_Y-%d_Z-%d_ROT-%03d.dat', ...
    boxsize, elemsize, x0, y0, z0, mod(round(rot*180/pi), 360));
fid = fopen(fullfile(topodir, fname), 'w');
fprintf(fid, '%g\n', reshape((wr .* zz)', n * n, 1));
fclose(fid);

% basin
zz = -z0 * exp(-(5*rr).^2 / 2);

figure(1)
clf
surface(x, y, wr.*zz);
colorbar
set(gca, 'DataAspectRatio', [1 1 z0]);
colormap(kelicol)

fname = sprintf('topo_GENERIC_BASIN_BOX-%d_ELEM-%d_X-%d_Y-%d_Z-%d_ROT-%03d.dat', ...
    boxsize, elemsize, x0, y0, z0, mod(round(rot*180/pi), 360));
fid = fopen(fullfile(topodir, fname), 'w');
fprintf(fid, '%g\n', reshape((wr .* zz)', n * n, 1));
fclose(fid);

% crator inside a hill
% equivalent to multiplying Gaussian with hypoberlic cosine
zz = z0 * exp(-(5*rr-1.5).^2/2) + z0 * exp(-(5*rr+1.5).^2/2);

figure(1)
clf
surface(x, y, wr.*zz);
colorbar
set(gca, 'DataAspectRatio', [1 1 z0]);
colormap(kelicol)

fname = sprintf('topo_GENERIC_CRATOR_BOX-%d_ELEM-%d_X-%d_Y-%d_Z-%d_ROT-%03d.dat', ...
    boxsize, elemsize, x0, y0, z0, mod(round(rot*180/pi), 360));
fid = fopen(fullfile(topodir, fname), 'w');
fprintf(fid, '%g\n', reshape((wr .* zz)', n * n, 1));
fclose(fid);

% hill sinking into the ground (flexure)
zz = z0 * sinc(pi * rr);

figure(1)
clf
surface(x, y, wr.*zz);
colorbar
set(gca, 'DataAspectRatio', [1 1 z0]);
colormap(kelicol)

fname = sprintf('topo_GENERIC_FLEXURE_BOX-%d_ELEM-%d_X-%d_Y-%d_Z-%d_ROT-%03d.dat', ...
    boxsize, elemsize, x0, y0, z0, mod(round(rot*180/pi), 360));
fid = fopen(fullfile(topodir, fname), 'w');
fprintf(fid, '%g\n', reshape((wr .* zz)', n * n, 1));
fclose(fid);

% elliptic hill
zz = z0 * exp(-((5 * xx_rot).^2 + (10 * yy_rot).^2) / 2);

figure(1)
clf
surface(x, y, wn.*zz);
colorbar
set(gca, 'DataAspectRatio', [1 1 z0]);
colormap(kelicol)

fname = sprintf('topo_GENERIC_ELLIPTIC_BOX-%d_ELEM-%d_X-%d_Y-%d_Z-%d_ROT-%03d.dat', ...
    boxsize, elemsize, x0, y0, z0, mod(round(rot*180/pi), 360));
fid = fopen(fullfile(topodir, fname), 'w');
fprintf(fid, '%g\n', reshape((wn .* zz)', n * n, 1));
fclose(fid);

% ridge
zz = z0 * exp(-(5 * xx_rot).^2 / 2);

figure(1)
clf
surface(x, y, wn.*zz);
colorbar
set(gca, 'DataAspectRatio', [1 1 z0]);
colormap(kelicol)

fname = sprintf('topo_GENERIC_RIDGE_BOX-%d_ELEM-%d_X-%d_Y-%d_Z-%d_ROT-%03d.dat', ...
    boxsize, elemsize, x0, y0, z0, mod(round(rot*180/pi), 360));
fid = fopen(fullfile(topodir, fname), 'w');
fprintf(fid, '%g\n', reshape((wn .* zz)', n * n, 1));
fclose(fid);

% trench
zz = -z0 * exp(-(5 * xx_rot).^2 / 2);

figure(1)
clf
surface(x, y, wn.*zz);
colorbar
set(gca, 'DataAspectRatio', [1 1 z0]);
colormap(kelicol)

fname = sprintf('topo_GENERIC_RIDGE_BOX-%d_ELEM-%d_X-%d_Y-%d_Z-%d_ROT-%03d.dat', ...
    boxsize, elemsize, x0, y0, z0, mod(round(rot*180/pi), 360));
fid = fopen(fullfile(topodir, fname), 'w');
fprintf(fid, '%g\n', reshape((wn .* zz)', n * n, 1));
fclose(fid);

% rolling ridges
zz = z0 * sin(2*pi* xx_rot * 1.2);

figure(1)
clf
surface(x, y, wn.*zz);
colorbar
set(gca, 'DataAspectRatio', [1 1 z0]);
colormap(kelicol)

fname = sprintf('topo_GENERIC_ROLLING_BOX-%d_ELEM-%d_X-%d_Y-%d_Z-%d_ROT-%03d.dat', ...
    boxsize, elemsize, x0, y0, z0, mod(round(rot*180/pi), 360));
fid = fopen(fullfile(topodir, fname), 'w');
fprintf(fid, '%g\n', reshape((wn .* zz)', n * n, 1));
fclose(fid);

% slope
zz = z0 * xx_rot;
zz = zz - mean(zz, 'all');

figure(1)
clf
surface(x, y, wn.*zz);
colorbar
set(gca, 'DataAspectRatio', [1 1 z0]);
colormap(kelicol)

fname = sprintf('topo_GENERIC_SLOPE_BOX-%d_ELEM-%d_X-%d_Y-%d_Z-%d_ROT-%03d.dat', ...
    boxsize, elemsize, x0, y0, z0, mod(round(rot*180/pi), 360));
fid = fopen(fullfile(topodir, fname), 'w');
fprintf(fid, '%g\n', reshape((wn .* zz)', n * n, 1));
fclose(fid);
end
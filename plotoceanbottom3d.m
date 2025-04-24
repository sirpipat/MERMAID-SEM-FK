function plotoceanbottom3d(ddir3d, unit, plotdir)
% PLOTOCEANBOTTOM3D(ddir3d, unit, plotdir)
%
% Plots the ocean bottom topography from a SPECFEM3D run.
%
% INPUT:
% ddir3d        directory for a SPECFEM3D run
% unit          unit of lengths in the plots (either 'm' [default] or 'km')
% plotdir       directory of the printed figure [Default: getenv('EPS')]
%
% Last modified by sirawich-at-princetonedu: 04/24/2025

defval('unit', 'm')
defval('plotdir', getenv('EPS'))

if strcmp(unit, 'km')
    divisor = 1000;
else
    divisor = 1;
    unit = 'm';
end

% read interface file
fname = fullfile(ddir3d, 'DATA', 'meshfem3D_files', 'interfaces.dat');
[itfs, ~] = loadinterfacefiles3d(fname);

% determine the grid coordinate
X = ((0:itfs{1}.NXI-1) * itfs{1}.SPACING_XI + itfs{1}.LON_MIN) / divisor;
Y = ((0:itfs{1}.NETA-1) * itfs{1}.SPACING_ETA + itfs{1}.LAT_MIN) / divisor;
Z = itfs{1}.Z / divisor;

% determine plot title
runname = removepath(ddir3d);
words = split(runname, '_');
if strcmpi(words{1}, 'RANDOM')
    s2 = str2double(cindeks(split(words{2}, '-'), 2));
    nu = str2double(cindeks(split(words{3}, '-'), 2));
    rho = str2double(cindeks(split(words{4}, '-'), 2));
    bottom = -5000;
    titlestring = sprintf('\\mu = %d | \\sigma^2 = %d | \\nu = %.1f | \\rho = %d', ...
        bottom / divisor, s2 / divisor^2, nu, rho / divisor);
    clim = bottom / divisor + [-4 4] * sqrt(s2);
else
    bottom = -str2double(cindeks(split(words{2}, '-'), 2));
    titlestring = sprintf('bottom = %d', bottom);
    clim = (bottom + [-200 200]) / divisor;
end

% slice for 2D cross-section profile along x=0
Z2D = Z(:, ceil(itfs{1}.NXI/2));

% read the bathymetry from 2D directory
ddir2d = fullfile(getenv('REMOTE2D'), sprintf('%s2D', ...
    cindeks(split(ddir3d, filesep), 'end-1')), runname);
[itfs2d, ~] = loadinterfacefile(fullfile(ddir2d, 'DATA', ...
    'interfaces_from3d.dat'));

figure(1)
clf
set(gcf, 'Units', 'inches', 'Position', [18 1 6 9])
ax_bath3d = subplot('Position', [0.15 0.29 0.80 0.69]);
imagesc(X, Y, Z)
grid on
hold on
plot([0 0], [Y(1) Y(end)], 'LineWidth', 1)
xlabel(sprintf('easting (%s)', unit))
ylabel(sprintf('northing (%s)', unit))
cb = colorbar;
colormap(gca, kelicol)
cb.Label.String = sprintf('elevation (%s)', unit);
cb.Location = 'southoutside';
title(titlestring)
set(ax_bath3d, 'TickDir', 'out', 'FontSize', 12, 'Box', 'on',  ...
    'DataAspectRatio', [1 1 1], ...
    'XLim', [-10000 10000] / divisor, ...
    'YLim', [-10000 10000] / divisor, ...
    'CLim', clim, ...
    'XTick', [-10000 -5000 0 5000 10000] / divisor, ...
    'XTickLabel', [-10000 -5000 0 5000 10000] / divisor, ...
    'YTick', [-10000 -5000 0 5000 10000] / divisor, ...
    'YTickLabel', [-10000 -5000 0 5000 10000] / divisor)
set(cb, 'FontSize', 12, 'Limit', clim, 'TickDirection', 'out')

ax_x2d = subplot('Position', [0.14 0.06 0.81 0.18]);
%plot(Ykm, Z2D, 'LineWidth', 1)
grid on
hold on
plot((itfs2d{2}.pts(:,1) - 10000) / divisor, ...
    (itfs2d{2}.pts(:,2) - 9600) / divisor, 'LineWidth', 1)
xlabel(sprintf('northing (%s)', unit))
ylabel(sprintf('elevation (%s)', unit))
xlim([-10000 10000] / divisor)
title('cross section along x = 0')
set(ax_x2d, 'TickDir', 'out', 'FontSize', 12, 'Box', 'on', ...
    'XTick', [-10000 -5000 0 5000 10000] / divisor, ...
    'XTickLabel', [-10000 -5000 0 5000 10000] / divisor)

set(ax_x2d, 'TickDir', 'out', 'FontSize', 12)

set(gcf, 'Renderer', 'painters')

% save figure
system(sprintf('mkdir -p %s', plotdir));
eps_filename = fullfile(plotdir, strcat(mfilename, '_', runname, '.epsc'));
pdf_filename = fullfile(plotdir, strcat(mfilename, '_', runname, '.pdf'));
print(eps_filename, '-depsc');
system(sprintf('epstopdf %s %s', eps_filename, pdf_filename));
system(sprintf('rm %s', eps_filename));
end
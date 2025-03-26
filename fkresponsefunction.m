function [t_rf, rf, d] = fkresponsefunction(ddir, id, reg, plt)
% [t_rf, rf, d] = FKRESPONSEFUNCTION(ddir, id, reg, plt)
%
% Computes response function due to flat layers over half-space. It
% deconvolve the injected wave (subjected to the layers) and the 
% source-time function used for injection (not subjected to the layers).
%
% INPUT:
% ddir          directory of a FK-SPECFEM3D simulation
% id            station ID ['Default': 1]
% reg           regularization method: 'water' [default] or 'damp'
% plt           whether to plot or not [Default: false]
%
% OUTPUT
% t_rf          time for response function
% rf            response funtion
% d             regularization value for spectral division
%
% SEE ALSO:
% SPECTRALDIVISION
%
% Last modified by sirawich-at-princeton.edu, 03/21/2025

defval('method', 'water')
defval('plt', false)

% read the injected wave
[t, ~, ~, vz, hdr] = readplotFKfile(strcat(ddir, ...
    sprintf('OUTPUT_FILES/plot_FK_Veloc.%d.dat', id)));
hdrwords = split(hdr);

if str2double(hdrwords{20}) < -4128
    z = cumsum(vz) * (t(2) - t(1));
else
    z = vz;
end

% read the source-time-function file
fkmodel = loadfkmodel(strcat(ddir, 'DATA/FKMODEL'));
if fkmodel.stf_type == 4
    % Source-time function is defined by a STF file
    [t0, z0] = readtimeseries(strcat(ddir, fkmodel.stf_file), false);
    
    % interpolate stf to the injected wave's time
    z1 = shannon(t0, z0, t);
else
    % Source-time funciton is Gaussian with a characteristic frequency
    z1 = fkmodel.fmax / sqrt(pi) * exp(-(t * fkmodel.fmax) .^ 2);
end


% deconvolve
[rf, ~, ~, d] = spectraldivision(detrend(z, 0), detrend(z1, 0), ...
    shanning(length(z), 0.1), reg, []);

% time for response function
t_rf = t - t(1);

if plt
    figure
    clf
    set(gcf, 'Units', 'inches', 'Position', [0 1 8 4])
    
    ax1 = subplot('Position', [0.09 0.70 0.88 0.22]);
    plot(t-t(1), z, 'LineWidth', 1, 'Color', 'b');
    grid on
    xlim([t_rf(1) t_rf(end)])
    nolabels(ax1, 1)
    ylabel('disp. (m)')
    set(ax1, 'TickDir', 'out', 'FontSize', 12)
    title(sprintf('(x,y,z) = (%s,%s,%s)', hdrwords{18}, hdrwords{19}, hdrwords{20}))
    
    ax2 = subplot('Position', [0.09 0.42 0.88 0.22]);
    plot(t-t(1), z1, 'LineWidth',1 ,'Color', [0 0.4 0]);
    grid on
    xlim([t_rf(1) t_rf(end)])
    nolabels(ax2, 1)
    ylabel('disp. (m)')
    set(ax2, 'TickDir', 'out', 'FontSize', 12)
    
    ax3 = subplot('Position', [0.09 0.14 0.88 0.22]);
    plot(t_rf, rf, 'LineWidth',1 ,'Color', 'k');
    grid on
    xlim([t_rf(1) t_rf(end)])
    xlabel('time (s)')
    ylabel('resp.')
    set(ax3, 'TickDir', 'out', 'FontSize', 12)
    
    set(gcf, 'Renderer', 'painters')
    figdisp(sprintf('%s_%02d', mfilename, id), [], [], 2, [], 'epstopdf');
end
end
function benchmarkflatocean(f_gaussian2d, f_gaussian3d, fc, npoles, filter_zobs, reg, tmax)
% BENCHMARKFLATOCEAN(f_gaussian2d, f_gaussian3d, fc, npoles, filter_zobs, reg, tmax)
%
% Makes plots showing the benchmarking FK-SPECFEM3D to SPECFEM2D on a flat
% bathymetry and Gaussian source-time function by reading the outputs from
% FK-SPECFEM3D and SPECFEM2D simulation runs.
%
% The first two plots are like Figure 8 from Pipatprathanporn & Simons
% (2024) paper for each simulation run. The other two plots are the
% comparison of the response functions in time and frequency domains.
%
% INPUT:
% f_gaussian2d      either 1 (default), 2, 4, or 10
% f_gaussian3d      either 1, 2 (default), 4, or 10
% fc                corner frequency for low-pass filtering [default: 2]
% npoles            number of poles for low-pass filtering  [default: 4]
% filter_zobs       whether to apply the filter to the first arrival of
%                   vertical displacement at ocean bottom [default: true]
% reg               regularization method for deconvolution
%                   either 'damp' (default) or 'water'
% tmax              maximum time to consider for analysis and visualization
%                   [default: 30]
%
% SEE ALSO:
% SPECTRALDIVISION
%
% Last modified by sirawich-at-princeton.edu, 10/31/2024

% define more PARAMETERS here
defval('f_gaussian2d', 1);
defval('f_gaussian3d', 2);
defval('fc', 1)
defval('npoles', 4);
defval('filter_zobs', true);
defval('reg', 'damp');
defval('tmax', 30);
ii_cutoff = 800;    % cutoff index
correct_demean = true;

% SPECFEM2D
if f_gaussian2d == 1
    outputdir2d = '/Users/sirawich/research/remote_specfem2d/flat_10936816_P0009_Gaussian1Hz/';
elseif f_gaussian2d == 2
    outputdir2d = '/Users/sirawich/research/remote_specfem2d/flat_10936816_P0009_Gaussian2Hz/';
elseif f_gaussian2d == 4
    outputdir2d = '/Users/sirawich/research/remote_specfem2d/flat_10936816_P0009_Gaussian4Hz/';
elseif f_gaussian2d == 10
    outputdir2d = '/Users/sirawich/research/remote_specfem2d/flat_10936816_P0009_Gaussian10Hz/';
else
    error('f_gaussian2d must be 1, 2, 4, or 10')
end
[tobs2d, zobs2d] = getarrivaltemplate(outputdir2d);
[tmh2d, pmh2d] = read_seismogram(fullfile(outputdir2d, 'OUTPUT_FILES', 'AA.S0001.PRE.semp'));
t = tobs2d - tobs2d(1);
dt = tobs2d(2) - tobs2d(1);
fs = 1 / dt;

% SPECFEM3D
if f_gaussian3d == 1
    outputdir3d = '/Users/sirawich/research/remote_specfem3d/FK-FLAT_10936816_P0009_1Hz/';
elseif f_gaussian3d == 2
    outputdir3d = '/Users/sirawich/research/remote_specfem3d/FK-FLAT_10936816_P0009_2Hz/';
elseif f_gaussian3d == 4
    outputdir3d = '/Users/sirawich/research/remote_specfem3d/FK-FLAT_10936816_P0009_4Hz/';
elseif f_gaussian3d == 10
    outputdir3d = '/Users/sirawich/research/remote_specfem3d/FK-FLAT_10936816_P0009_10Hz/';
else
    error('f_gaussian3d must be 1, 2, 4, or 10')
end
[tobs3d, zobs3d] = read_seismogram(fullfile(outputdir3d, 'OUTPUT_FILES', 'AA.OBS01.HXZ.semd'));
[tmh3d, pmh3d] = read_seismogram(fullfile(outputdir3d, 'OUTPUT_FILES', 'MH.P0009.HXP.semp'));
t3d = tobs3d - tobs3d(1);

% resample to the same sampling rate
zobs3d = shannon(t3d, zobs3d, t);
pmh3d = shannon(t3d, pmh3d, t);

% lowpass filtering
if fc > 0
    % applying a filter to zobs (template, input) will perturb the nice
    % Gaussian source-time function. 
    if filter_zobs
        % keeps the mean and apply after lowpass filter since it demeans the
        zobs3d_mean = mean(zobs3d);
        zobs2d_mean = mean(zobs2d);
        zobs3dlp = lowpass(zobs3d, fs, fc, npoles, 2, 'butter', 'constant') + correct_demean * zobs3d_mean;
        zobs2d = lowpass(zobs2d, fs, fc, npoles, 2, 'butter', 'constant') + correct_demean * zobs2d_mean;
    else
        zobs3dlp = zobs3d;
    end
    % keeps the mean and apply after lowpass filter since it demeans the
    pmh3d_mean = mean(pmh3d);
    pmh2d_mean = mean(pmh2d);
    pmh3dlp = lowpass(pmh3d, fs, fc, npoles, 2, 'butter', 'constant') + correct_demean * pmh3d_mean;
    pmh2d = lowpass(pmh2d, fs, fc, npoles, 2, 'butter', 'constant') + correct_demean * pmh2d_mean;
else
    zobs3dlp = zobs3d;
    pmh3dlp = pmh3d;
end

% removes the trailing signals
zobstemp = zobs3dlp(1:ii_cutoff) .* shanning(ii_cutoff, 0.1);
zobstemp = [zobstemp; zeros(length(zobs3dlp)-ii_cutoff, 1)];

% applies tiny random noise to the seismograms at the ocean bottom to avoid
% error when calling pchave vi aspecdensplot.
zobs2d = zobs2d + randn(size(zobs2d)) * 1e-24 * std(zobs2d);
zobstemp = zobstemp + randn(size(zobstemp)) * 1e-24 * std(zobstemp);

% computes response functions
r2d = spectraldivision(pmh2d, zobs2d, [], reg, []);
r3d = spectraldivision(pmh3dlp, zobstemp, [], reg, []);

% Makes a similar plot to Fig. 8 of Pipatprathanporn & Simons (2024)
figure(10)
clf
set(gcf, 'Unit', 'inches', 'Position', [14 6.1 9 4]);
subplot('Position', [0.06 0.16 0.88 0.76])
plot(t, pmh2d / max(abs(pmh2d)), 'LineWidth', 1, 'Color', [0 0.5 1]);
hold on
grid on
plot(t, zobs2d / max(abs(zobs2d)) - 2, 'LineWidth', 1, 'Color', [0 0.7 0]);
plot(t, r2d / max(abs(r2d)) - 4, 'LineWidth', 1, 'Color', 'k');
xlabel('time (s)')
xlim([0 tmax])
legend('MERMAID Pressure', 'OBS Z-disp 1st arrival', ...
    sprintf('Response function (%s)', reg), ...
    'Position', [0.6794 0.3535 0.1928 0.1371]);
set(gca, 'TickDir', 'out', 'FontSize', 12)
title(sprintf('SPECFEM2D: STF-Gaussian %d Hz, lp -c %d -n %d -p 2', ...
    f_gaussian2d, fc, npoles))
nolabels(gca, 2)
set(gcf, 'Renderer', 'painters')
figdisp(sprintf('%s_SPECFEM2D_TEST_Gaussian_%dHz_Response', mfilename, ...
    f_gaussian2d), [], [], 2, [], 'epstopdf')

figure(11)
clf
set(gcf, 'Unit', 'inches', 'Position', [14 1.1 9 4]);
subplot('Position', [0.06 0.16 0.88 0.76])
plot(t, pmh3dlp / max(abs(pmh3dlp)), 'LineWidth', 1, 'Color', [0 0.5 1]);
hold on
grid on
plot(t, zobstemp / max(abs(zobstemp)) - 2, 'LineWidth', 1, 'Color', [0 0.7 0]);
plot(t, r3d / max(abs(r3d)) - 4, 'LineWidth', 1, 'Color', [.75 0 0]);
xlabel('time (s)')
xlim([0 tmax])
legend('MERMAID Pressure', 'OBS Z-disp 1st arrival', ...
    sprintf('Response function (%s)', reg), ...
    'Position', [0.6794 0.3535 0.1928 0.1371]);
set(gca, 'TickDir', 'out', 'FontSize', 12)
title(sprintf('SPECFEM3D: STF-Gaussian %d Hz, lp -c %d -n %d -p 2', ...
    f_gaussian3d, fc, npoles))
nolabels(gca, 2)
set(gcf, 'Renderer', 'painters')
figdisp(sprintf('%s_FK-SPECFEM3D_TEST_Gaussian_%dHz_Response', ...
    mfilename, f_gaussian3d), [], [], 2, [], 'epstopdf');


% Make a PSD plot comparing the PSD of two response functions
figure(13)
clf
set(gcf, 'Unit', 'inches', 'Position', [14 11.1 9 5]);
subplot('Position', [0.10 0.13 0.80 0.68])
p2d = specdensplot(r2d(t <= tmax), 500, fs, 500, 70, 10, 's');
hold on
grid on
p3d = specdensplot(r3d(t <= tmax), 500, fs, 500, 70, 10, 's');
for ii = 1:length(p2d)
    p2d(ii).Color = [0 0 0];
    p3d(ii).Color = [.75 0 0];
    if ii == 1
        p2d(ii).LineWidth = 1;
        p3d(ii).LineWidth = 1;
    elseif ii == 4
        p2d(ii).MarkerFaceColor = p2d(ii).Color;
        p3d(ii).MarkerFaceColor = p3d(ii).Color;
        p2d(ii).Marker = 'o';
        p3d(ii).Marker = 'o';
        p2d(ii).MarkerSize = 4;
        p3d(ii).MarkerSize = 4;
    end
end
set(gca, 'TickDir', 'out', 'FontSize', 12, 'XTick', [.2 1 10 49.99], ...
    'XTickLabel', {'0.2', '1', '10', '50'})
legend([p2d(1) p3d(1)], sprintf('SPECFEM2D (%d Hz)', f_gaussian2d), ...
    sprintf('SPECFEM3D (%d Hz)', f_gaussian3d))
ax = gca;
ax2 = doubleaxes(ax);
inverseaxis(ax2.XAxis, 'period (s)')
ax2.YAxis.Label.String = ax.YAxis.Label.String;
[~, y_title] = norm2trueposition(gca, 0, 1.16);
title(ax, 'Response function PSD', 'Position', [3.1623 y_title 0])
set(gcf, 'Renderer', 'painters')
figdisp(sprintf('%s_SPECFEM3D_Compare_Responses_PSD_2D_%dHz_3D_%dHz', ...
    mfilename, f_gaussian2d, f_gaussian3d), [], [], 2, [], 'epstopdf')

% Make plot response function comparison
scale = 1.2 * max([max(abs(r2d)) max(abs(r3d)) max(abs(r3d - r2d))]);
figure(14)
clf
set(gcf, 'Unit', 'inches', 'Position', [5 1.1 9 6]);
subplot('Position', [0.06 0.72 0.88 0.22])
plot(t, r3d, 'LineWidth', .5, 'Color', [.9 0.6 0.6])
hold on
plot(t, r2d, 'LineWidth', 1, 'Color', 'k')
grid on
xlim([0 tmax])
ylim([-1 1] * scale)
title(sprintf('Response function: SPECFEM2D %d Hz', f_gaussian2d))
set(gca, 'TickDir', 'out', 'FontSize', 12)
nolabels(gca, 1)

subplot('Position', [0.06 0.42 0.88 0.22])
plot(t, r2d, 'LineWidth', .5, 'Color', [.6 .6 .6])
hold on
plot(t, r3d, 'LineWidth', 1, 'Color', [.75 0 0])
grid on
xlim([0 tmax])
ylim([-1 1] * scale)
title(sprintf('Response function: SPECFEM3D %d Hz', f_gaussian3d))
set(gca, 'TickDir', 'out', 'FontSize', 12)
nolabels(gca, 1)


subplot('Position', [0.06 0.12 0.88 0.22])
plot(t, r3d - r2d, 'LineWidth', 1, 'Color', 'b')
grid on
hold on
xlim([0 tmax])
ylim([-1 1] * scale)
hline(gca, [-1; 1] * max(abs(r3d(t <= tmax) - r2d(t <= tmax))), ...
    'LineStyle', '-', 'LineWidth', 0.25, 'Color', [0.8 0.4 0]);
[x_text_norm, y_text_norm] = true2normposition(gca, 2, ...
    max(abs(r3d(t <= tmax) - r2d(t <= tmax))));
[x_text, y_text] = norm2trueposition(gca, x_text_norm, y_text_norm + 0.12);
e_text = floor(log10(max(abs(r3d(t <= tmax) - r2d(t <= tmax)))));
b_text = max(abs(r3d(t <= tmax) - r2d(t <= tmax))) / 10^e_text;
t_text = sprintf('%.2f \\times 10^{%d}', b_text, e_text);
text(x_text, y_text, t_text, 'FontSize', 12, 'Color', [0.8 0.4 0])
xlabel('time (s)')
title(sprintf('Response function difference: SPECFEM3D (%d Hz) - SPECFEM2D (%d Hz)', ...
    f_gaussian2d, f_gaussian3d))
set(gca, 'TickDir', 'out', 'FontSize', 12)
set(gcf, 'Renderer', 'painters')
figdisp(sprintf('%s_SPECFEM3D_Compare_Responses_TIME_2D_%dHz_3D_%dHz', ...
    mfilename, f_gaussian2d, f_gaussian3d), [], [], 2, [], 'epstopdf')

% Make a PSD plot comparing the PSD of two response functions
figure(15)
clf
set(gcf, 'Unit', 'inches', 'Position', [5 11.1 9 5]);
subplot('Position', [0.10 0.13 0.80 0.68])
p2d = specdensplot(zobs2d(1:1500) / std(zobs2d(1:1500)), 500, fs, 500, 70, 10, 's');
hold on
grid on
p3d = specdensplot(zobstemp(1:1500) / std(zobstemp(1:1500)), 500, fs, 500, 70, 10, 's');
for ii = 1:length(p2d)
    p2d(ii).Color = [0 0 0];
    p3d(ii).Color = [.75 0 0];
    if ii == 1
        p2d(ii).LineWidth = 1;
        p3d(ii).LineWidth = 1;
    elseif ii == 4
        p2d(ii).MarkerFaceColor = p2d(ii).Color;
        p3d(ii).MarkerFaceColor = p3d(ii).Color;
        p2d(ii).Marker = 'o';
        p3d(ii).Marker = 'o';
        p2d(ii).MarkerSize = 4;
        p3d(ii).MarkerSize = 4;
    end
end
set(gca, 'TickDir', 'out', 'FontSize', 12, 'XTick', [.2 1 10 49.99], ...
    'XTickLabel', {'0.2', '1', '10', '50'})
legend([p2d(1) p3d(1)], sprintf('SPECFEM2D (%d Hz)', f_gaussian2d), ...
    sprintf('SPECFEM3D (%d Hz)', f_gaussian3d))
ax = gca;
ax2 = doubleaxes(ax);
inverseaxis(ax2.XAxis, 'period (s)')
ax2.YAxis.Label.String = ax.YAxis.Label.String;
[~, y_title] = norm2trueposition(ax, 0.50, 1.16);
title(ax, 'OBS Z-disp function PSD', 'Position', [3.1623 y_title 0])
set(gcf, 'Renderer', 'painters')
figdisp(sprintf('%s_SPECFEM3D_Compare_zdisp_PSD_2D_%dHz_3D_%dHz', ...
    mfilename, f_gaussian2d, f_gaussian3d), [], [], 2, [], 'epstopdf')
end
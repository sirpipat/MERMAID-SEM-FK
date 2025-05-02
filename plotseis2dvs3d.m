function [tshift1, tshift2, cc1, cc2, ss1, ss2] = plotseis2dvs3d(ddir2d, ddir3d, useflag, plt)
% [tshift1, tshift2, cc1, cc2, ss1, ss2] = PLOTSEIS2DVS3D(ddir2d, ddir3d, useflag, plt)
%
% Plots the following:
% 1) SPECFEM2D seismograms vs scaled, time shifted SPECFEM3D seismograms
% 2) Ray-theory seismograms (flat-ocean) vs scaled, time shifted SPECFEM3D 
%    seismograms
%
% 3)-4) Correlation coefficients vs timeshifts for 1) and 2)
%
% For each plot, the seismograms are split to six windows centering at the
% midpoint between the upgoing arrival and the next downgoing arrival.
%
% INPUT:
% ddir2d        directory for a SPECFEM2D run OR a cell of directories
% ddir3d        directory for a SPECFEM3D run OR a cell of directories
% useflag       either 'obs' or 'hydrophone' [default]
% plt           whether to plot or not [Default: true]
%
% OUTPUT:
% tshift1       timeshift for SPECFEM3D to align SPECFEM2D seimosgram
% tshift2       timeshift for SPECFEM3D to align Ray-theory seimosgram
% cc1           max correlation coefficient of SPECFEM3D and SPECFEM2D seimosgrams
% cc2           max correlation coefficient of SPECFEM3D and Ray-theory seimosgrams
% ss1           scaling for SPECFEM3D to align SPECFEM2D seimosgram
% ss2           scaling for SPECFEM3D to align Ray-theory seimosgram
%
% SEE ALSO:
% PLOTLOCAMP, CCSCALE
%
% Last modified by sirawich-at-princeton.edu, 05/02/2025

defval('useflag', 'hydrophone')
defval('plt', true)

%% loop over directories if they are cell arrays
if iscellstr(ddir2d) || iscellstr(ddir3d)
    % input validation
    if ~iscellstr(ddir2d) || ~iscellstr(ddir3d)
        error('Directories for 2D and 3D mush be both strings or both cell arrays');
    end
    n1 = length(ddir2d);
    n2 = length(ddir3d);
    if n1~= n2
        error('Lengths of the directory lists mush be equal.')
    end
    
    % determine the hash file name form input parameters
    defval('sname', sprintf('%s_%s.mat', mfilename, ...
        hash([n1 n2 double(cell2mat(ddir2d)) double(cell2mat(ddir3d)) ...
        useflag], 'SHA-1')))
    
    pname = fullfile(getenv('IFILES'), 'HASHES', sname);
    
    % Loop correlations if the hash file does not exist or the plots are
    % requested.
    if plt || ~exist(pname, 'file')
        tshift1 = zeros(n1,6);
        tshift2 = zeros(n1,6);
        cc1 = zeros(n1,6);
        cc2 = zeros(n1,6);
        ss1 = zeros(n1,6);
        ss2 = zeros(n1,6);

        for ii = 1:n1
            [tshift1(ii,:), tshift2(ii,:), cc1(ii,:), cc2(ii,:), ...
                ss1(ii,:), ss2(ii,:)] = plotseis2dvs3d(ddir2d{ii}, ...
                ddir3d{ii}, useflag, plt);
        end
        
        % save
        fprintf('save the output to a file to %s ...\n', pname);
        save(pname, 'tshift1', 'tshift2', 'cc1', 'cc2', 'ss1', 'ss2');
    % Otherwise, load the existing hash file to speed up the process.
    else
        % load
        fprintf('found the save in a file in %s\n', pname);
        load(pname, 'tshift1', 'tshift2', 'cc1', 'cc2', 'ss1', 'ss2');
    end
    
    % plot statistics
    figure(4)
    clf
    set(gcf, 'Unit', 'inches', 'Position', [12 1 6 8]);
    
    tshift1_grid = zeros(41, 6);
    tshift2_grid = zeros(41, 6);
    cc1_grid = zeros(41, 6);
    cc2_grid = zeros(41, 6);
    ss1_grid = zeros(61, 6);
    ss2_grid = zeros(61, 6);
    
    for ii = 1:n1
        for jj = 1:6
            ii_tshift1 = round(tshift1(ii, jj) * 20) + 21;
            ii_tshift1 = min(max(ii_tshift1, 1), 41);
            tshift1_grid(ii_tshift1, jj) = tshift1_grid(ii_tshift1, jj) + 1;
            
            ii_tshift2 = round(tshift2(ii, jj) * 20) + 21;
            ii_tshift2 = min(max(ii_tshift2, 1), 41);
            tshift2_grid(ii_tshift2, jj) = tshift2_grid(ii_tshift2, jj) + 1;
            
            ii_cc1 = round(cc1(ii, jj) * 20) + 21;
            cc1_grid(ii_cc1, jj) = cc1_grid(ii_cc1, jj) + 1;
            
            ii_cc2 = round(cc1(ii, jj) * 20) + 21;
            cc2_grid(ii_cc2, jj) = cc2_grid(ii_cc2, jj) + 1;
            
            ii_ss1 = round(log2(ss1(ii, jj)) * 10) + 31;
            ii_ss1 = min(max(ii_ss1, 1), 61);
            ss1_grid(ii_ss1, jj) = ss1_grid(ii_ss1, jj) + 1;
            
            ii_ss2 = round(log2(ss2(ii, jj)) * 10) + 31;
            ii_ss2 = min(max(ii_ss2, 1), 61);
            ss2_grid(ii_ss2, jj) = ss2_grid(ii_ss2, jj) + 1;
        end
    end
    
    maxval = max([max(tshift1_grid, [], 'all') ...
        max(tshift2_grid, [], 'all') ...
        max(cc1_grid, [], 'all') ...
        max(cc2_grid, [], 'all') ...
        max(ss1_grid, [], 'all') ...
        max(ss2_grid, [], 'all')]);
    
    subplot('Position', [0.13 0.70 0.35 0.23])
    imagesc(1:6, (-1:0.05:1), tshift1_grid);
    axis xy
    colormap((1:-0.005:0)' * [1 1 1]);
    grid on
    hold on
    errorbar(1:6, median(tshift1, 1), 1.00 * std(tshift1, 1), ...
        1.00 * std(tshift1, 1), '-o', 'LineWidth', 1, 'Color', rgbcolor('2'))
    ylabel('timeshift (s)')
    set(gca, 'TickDir', 'out', 'FontSize', 12, 'CLim', [0 maxval], ...
        'XTick', 1:6)
    title('SPECFEM2D vs SPECFEM3D')
    
    subplot('Position', [0.62 0.70 0.35 0.23])
    imagesc(1:6, (-1:0.05:1), tshift2_grid);
    axis xy
    colormap((1:-0.005:0)' * [1 1 1]);
    grid on
    hold on
    errorbar(1:6, median(tshift2, 1), 1.00 * std(tshift2, 1), ...
        1.00 * std(tshift2, 1), '-o', 'LineWidth', 1, 'Color', rgbcolor('2'))
    ylabel('timeshift (s)')
    set(gca, 'TickDir', 'out', 'FontSize', 12, 'CLim', [0 maxval], ...
        'XTick', 1:6)
    title('FLAT vs SPECFEM3D')
    
    subplot('Position', [0.13 0.42 0.35 0.23])
    imagesc(1:6, (-1:0.05:1), cc1_grid);
    axis xy
    colormap((1:-0.005:0)' * [1 1 1]);
    grid on
    hold on
    errorbar(1:6, median(cc1, 1), 1.00 * std(cc1, 1), ...
        1.00 * std(cc1, 1), '-o', 'LineWidth', 1, 'Color', rgbcolor('2'))
    ylabel('correlation coefficient')
    ylim([0.475-0.5 1.025]);
    set(gca, 'TickDir', 'out', 'FontSize', 12, 'CLim', [0 maxval], ...
        'XTick', 1:6)
    
    subplot('Position', [0.62 0.42 0.35 0.23])
    imagesc(1:6, (-1:0.05:1), cc2_grid);
    axis xy
    colormap((1:-0.005:0)' * [1 1 1]);
    grid on
    hold on
    errorbar(1:6, median(cc2, 1), 1.00 * std(cc2, 1), ...
        1.00 * std(cc2, 1), '-o', 'LineWidth', 1, 'Color', rgbcolor('2'))
    ylabel('correlation coefficient')
    ylim([0.475-0.5 1.025]);
    set(gca, 'TickDir', 'out', 'FontSize', 12, 'CLim', [0 maxval], ...
        'XTick', 1:6)
    
    subplot('Position', [0.13 0.14 0.35 0.23])
    imagesc(1:6, linspace(-3, 3, 61), ss1_grid);
    axis xy
    colormap((1:-0.005:0)' * [1 1 1]);
    grid on
    hold on
    errorbar(1:6, log2(median(ss1, 1)), ...
        -log2(1 - 1.00 * std(ss1, 1) ./ median(ss1, 1)), ...
        log2(1 + 1.00 * std(ss1, 1)./ median(ss1, 1)), '-o', ...
        'LineWidth', 1, 'Color', rgbcolor('2'))
    xlabel('cycle')
    ylabel('scaling')
    set(gca, 'TickDir', 'out', 'FontSize', 12, 'CLim', [0 maxval], ...
        'XTick', 1:6, 'YTick', -3:3, 'YTickLabel', 2.^(-3:3))
    
    subplot('Position', [0.62 0.14 0.35 0.23])
    imagesc(1:6, linspace(-3, 3, 61), ss2_grid);
    axis xy
    colormap((1:-0.005:0)' * [1 1 1]);
    grid on
    hold on
    errorbar(1:6, log2(median(ss2, 1)), ...
        -log2(1 - 1.00 * std(ss2, 1) ./ median(ss2, 1)), ...
        log2(1 + 1.00 * std(ss2, 1)./ median(ss2, 1)), '-o', ...
        'LineWidth', 1, 'Color', rgbcolor('2'))
    xlabel('cycle')
    ylabel('scaling')
    set(gca, 'TickDir', 'out', 'FontSize', 12, 'CLim', [0 maxval], ...
        'XTick', 1:6, 'YTick', -3:3, 'YTickLabel', 2.^(-3:3))
    
    % add colorbar
    ax = subplot('Position', [0.12 0.09 0.85 0.03]);
    im = imagesc(1:6, (-1:0.01:1), cc1_grid);
    colormap((1:-0.005:0)' * [1 1 1]);
    cb = colorbar(ax, 'southoutside');
    set(im, 'Visible', 'off')
    set(cb.Label, 'String', 'counts', 'FontSize', 12)
    set(ax, 'FontSize', 12, 'Color', 'none', 'CLim', [0 maxval])
    set(ax.XAxis, 'Visible', 'off')
    set(ax.YAxis, 'Visible', 'off')
    
    % add title
    ax = subplot('Position', [0.12 0.96 0.85 0.01]);
    set(ax, 'FontSize', 12, 'Color', 'none')
    set(ax.XAxis, 'Visible', 'off')
    set(ax.YAxis, 'Visible', 'off')
    
    words = split(removepath(ddir3d{1}), '_');
    s2 = str2double(cindeks(split(words{2}, '-'), 2));
    nu = str2double(cindeks(split(words{3}, '-'), 2));
    rho = str2double(cindeks(split(words{4}, '-'), 2));
    titlestring = sprintf('\\mu = %d | \\sigma^2 = %d | \\nu = %.1f | \\rho = %d', ...
        -5000, s2, nu, rho);
    title(titlestring)
    return
end

%% Read seismograms from SPECFEM2D run
try
    [tobs2d, zobs2d] = read_seismogram(fullfile(ddir2d, 'OUTPUT_FILES', ...
        'AB.S0001.BXZ.semd'));
catch ME
    par_file = fullfile(ddir2d, 'DATA', 'Par_file');
    source_file = fullfile(ddir2d, 'DATA', 'SOURCE');
    tobs2d = specfem2dtime(par_file, 2, source_file);
    seis_file = fullfile(ddir2d, 'OUTPUT_FILES', 'Uz_file_single_d.bin');
    try
        data = freadseismograms(seis_file, par_file);
    catch ME
        seis_file = fullfile(ddir2d, 'OUTPUT_FILES', ...
            'Uz_file_double_d.bin');
        data = freadseismograms(seis_file, par_file);
    end
    zobs2d = data(:,2);
end
try
    [tmh2d, pmh2d] = read_seismogram(fullfile(ddir2d, 'OUTPUT_FILES', ...
        'AA.S0001.PRE.semp'));
catch ME
    par_file = fullfile(ddir2d, 'DATA', 'Par_file');
    source_file = fullfile(ddir2d, 'DATA', 'SOURCE');
    tmh2d = specfem2dtime(par_file, 2, source_file);
    seis_file = fullfile(ddir2d, 'OUTPUT_FILES', 'Up_file_single_p.bin');
    try
        data = freadseismograms(seis_file, par_file);
    catch ME
        seis_file = fullfile(ddir2d, 'OUTPUT_FILES', ...
            'Up_file_double_p.bin');
        data = freadseismograms(seis_file, par_file);
    end
    pmh2d = data(:,1);
end

% remove high-frequency artifact
zobs2d = lowpass(zobs2d, (length(tobs2d) - 1) / ...
    (tobs2d(end) - tobs2d(1)), 2, 2, 2, 'butter', 'linear');
pmh2d = lowpass(pmh2d, (length(tmh2d) - 1) / (tmh2d(end) - tmh2d(1)), ...
    2, 2, 2, 'butter', 'linear');

%% Read seismograms from SPECFEM3D run
[tobs3d, zobs3d] = read_seismogram(fullfile(ddir3d, 'OUTPUT_FILES', ...
    'AA.OBS01.HXZ.semd'));
[tmh3d, pmh3d] = read_seismogram(fullfile(ddir3d, 'OUTPUT_FILES', ...
    'MH.P0009.HXP.semp'));

% remove high-frequency artifact
zobs3d = lowpass(zobs3d, (length(tobs3d) - 1) / ...
    (tobs3d(end) - tobs3d(1)), 2, 2, 2, 'butter', 'linear');
pmh3d = lowpass(pmh3d, (length(tmh3d) - 1) / (tmh3d(end) - tmh3d(1)), ...
    2, 2, 2, 'butter', 'linear');

% Get the geometry
fkmodel = loadfkmodel(fullfile(ddir3d, 'DATA', 'FKMODEL'));
rho_f = fkmodel.layers{1}.rho;
vp_f = fkmodel.layers{1}.vp;
rho_s = fkmodel.layers{2}.rho;
vp_s = fkmodel.layers{2}.vp;
vs_s = fkmodel.layers{2}.vs;
z_interface = -fkmodel.layers{2}.ztop;
theta = fkmodel.theta;

%% Get stations locations
[n, name, network, sx, sy, sz] = readstations3d(fullfile(ddir3d, ...
    'DATA', 'STATIONS'));

% station elevation
z_station_mh = -indeks(sz(strcmp(network, 'MH')), 1);
z_station_obs = -indeks(sz(strcmp(network, 'AA')), 1);

%% choose which seismograms to use
if strcmpi(useflag, 'hydrophone')
    t2d = tmh2d;
    seis2d = pmh2d;
    t3d = tmh3d;
    seis3d = pmh3d;
    z_station = z_station_mh;
else
    t2d = tobs2d;
    seis2d = zobs2d;
    t3d = tobs3d;
    seis3d = zobs3d;
    z_station = z_station_obs;
end

% correct time for SPECFEM2D
% the source is at (200, 720), but it should be at (0,0)
% TODO: Correct this by setting SPECFEM2D geometry such that first source
% is at (-14000, -10000) like in SPECFFEM3D.
sources = loadsource(fullfile(ddir2d, 'DATA', 'SOURCE'));
xs = sources{1}.xs + 4000;
zs = sources{1}.zs + 400;
t0_2d = (xs * sin(theta * pi/180) + zs * cos(theta * pi/180)) / vp_s;
t2d = t2d + t0_2d;
    
%% calculate the ray theory arrivals
[tarr, xarr, zarr, parr] = fluidsolidarrivals(z_interface, z_station, ...
    rho_f, vp_f, rho_s, vp_s, vs_s, theta * pi/180, t3d(end));

% compute the time of the first arrival
meshparams3d = loadmeshparfile3d(fullfile(ddir3d, 'DATA', 'meshfem3D_files', 'Mesh_par_file'));
z0 = -meshparams3d.DEPTH_BLOCK_KM * 1000;
x0 = meshparams3d.LONGITUDE_MIN;
if strcmpi(useflag, 'hydrophone')
    x_station = indeks(sx(strcmp(network, 'MH')), 1);
    t0 = ((x_station - x0) * sin(theta * pi/180) + ... 
        (-z_interface - z0) * cos(theta * pi/180)) / vp_s + ...
        (-z_station_mh - (-z_interface)) * ...
        cos(asin(vp_f/vp_s * sin(theta * pi/180))) / vp_f;
else
    x_station = indeks(sx(strcmp(network, 'AA')), 1);
    t0 = ((x_station - x0) * sin(theta * pi/180) + ... 
        (-z_station_obs - z0) * cos(theta * pi/180)) / vp_s;
end
tarr = tarr + t0;

% calculate the ray theory waveforms
if fkmodel.stf_type == 4
    [t, x] = read_seismogram(fullfile(ddir3d, fkmodel.stf_file));
    x = shannon(t, x, t3d);
else
    x = fkmodel.fmax / sqrt(pi) * exp(-(t3d * fkmodel.fmax).^2);
end
if strcmpi(useflag, 'hydrophone')
    zarr(parr < 0) = -zarr(parr < 0);
end
w = placeweight(tarr, zarr, t3d - t3d(1));
if strcmpi(useflag, 'hydrophone')
    seis = indeks(diff(conv(x, w)) / (t3d(2) - t3d(1)), 1:length(x));
else
    seis = indeks(conv(x, w), 1:length(x));
end

%% rescaling amplitudes
% TODO: find a way to normalize the amplitudes in SPECFEM2D, SPECFEM3D-FK,
% and ray-theory prediction
% Number to be figured out (i.e. magic numbers)
% 18.7681 -- conversion from ray-theory to FK-injection in SPECFEM3D
% 2.2158e+23 / 18.7681 -- conversion from SPECFEM2D to SPECFEM3D
%
% These numbers are determined from amplitude of first arrival in the
% vertical displacement seismograms at the flat ocean bottom. The depth is
% 5000m, incident angles of 0, 10, and 40 degrees, the frequencies of 1, 2,
% and 4 Hz.
%
% [2025/04/09] The amplitude from FK-injection is doubled when the FK
% sampling frequency reduced from 20 to 10 Hz.

% SPECFEM2D
seis2d = seis2d * fkmodel.fmax^2 * (2.2158e+23 / 18.7681);

% SPECFEM3D
seis3d = seis3d / 18.7681; %19.99;

% Ray Theory
if strcmpi(useflag, 'hydrophone')
    % multiply the cosine for z-component
    seis = seis * cos(theta * pi/180);
    
    % now take care of displacement amplitude coefficient of P-to-acoustic
    % wave
    [~, D, ~] = fluidsolidcoefficients2(rho_f, vp_f, rho_s, vp_s, vs_s, ...
        theta * pi/180);
    seis = seis * D(3);
    
    % now convert from du_z / dt to P = -\lambda (\Nabla \cdot u)
    seis = seis * rho_f * vp_f / cos(asin(sin(theta*pi/180) * vp_f/vp_s));
else
    % multiply the cosine for z-component
    seis = seis * cos(theta * pi/180);
end

%% determine the arrivals in the seismograms
[pks3d, locs3d] = findpeakstopbottom(seis3d, t3d, 'MinPeakHeight', ...
    0.1 * max(abs(seis3d)), 'MinPeakProminence', 0.05 * max(abs(seis3d)));
[pks2d, locs2d] = findpeakstopbottom(seis2d, t2d, 'MinPeakHeight',...
    0.1 * max(abs(seis2d)), 'MinPeakProminence', 0.05 * max(abs(seis2d)));
[pks, locs] = findpeakstopbottom(seis, t3d, 'MinPeakHeight', ...
    0.1 * max(abs(seis)), 'MinPeakProminence', 0.05 * max(abs(seis)));
locs_full = locs;

%% remove unmatched peaks
% if strcmpi(use, 'hydrophone')
%     pks3d  = pks3d(1:25);
%     locs3d = locs3d(1:25);
%     pks2d  = pks2d([1:8 10:26]);
%     locs2d = locs2d([1:8 10:26]);
%     locs   = locs(1:25);
%     pks    = pks(1:25);
% else
%     pks2d = pks2d([1:2 4:6]');
%     locs2d = locs2d([1:2 4:6]');
% end

% locations of peaks relative to the first peak's location
locs_rel = [0; locs(2:end) - locs(1:end-1)];
locs2d_rel = [0; locs2d(2:end) - locs2d(1:end-1)];
locs3d_rel = [0; locs3d(2:end) - locs3d(1:end-1)];

% remove any peaks that do not quite match the prediction
ii = 1;
while ii <= length(locs_rel)
    while ii <= length(locs2d_rel) && abs(locs2d_rel(ii) - locs_rel(ii)) > 0.1
        if ii + 1 <= length(locs2d_rel)
            locs2d_rel(ii+1) = locs2d_rel(ii+1) + locs2d_rel(ii);
        end
        locs2d_rel(ii) = [];
        locs2d(ii) = [];
        pks2d(ii) = [];
    end
    while ii <= length(locs3d_rel) && abs(locs3d_rel(ii) - locs_rel(ii)) > 0.5
        if ii + 1 <= length(locs3d_rel)
            locs3d_rel(ii+1) = locs3d_rel(ii+1) + locs3d_rel(ii);
        end
        locs3d_rel(ii) = [];
        locs3d(ii) = [];
        pks3d(ii) = [];
    end
    % update the condition
    ii = ii + 1;
end

%% output containers
tshift1 = zeros(1,6);
tshift2 = zeros(1,6);
cc1 = zeros(1,6);
cc2 = zeros(1,6);
ss1 = zeros(1,6);
ss2 = zeros(1,6);

%% plot CC
% SPECFEM2D vs SPECFEM3D
if plt
    fig2 = figure(2);
    clf
    set(gcf, 'Units', 'inches', 'Position', [0 1 6 8])
    % list of axes
    ax = subplot(6, 2, 1);
    fig2axs = repmat(ax, 13, 1);
    for ii = 2:12
        fig2axs(ii) = subplot(6, 2, ii);
    end
    fig2axs(13) = subplot('Position', [0.09 0.95 0.82 0.01]);

    % FLAT v SPECFEM3D
    fig3 = figure(3);
    clf
    set(gcf, 'Units', 'inches', 'Position', [6 1 6 8])
    % list of axes
    ax = subplot(6, 2, 1);
    fig3axs = repmat(ax, 13, 1);
    for ii = 2:12
        fig3axs(ii) = subplot(6, 2, ii);
    end
    fig3axs(13) = subplot('Position', [0.09 0.95 0.82 0.01]);
end

% reverberation period
if strcmp(useflag, 'hydrophone')
    width = tarr(3) - tarr(1);
else
    width = tarr(2) - tarr(1);
end

% central position of the first+second arrivals
if strcmp(useflag, 'hydrophone')
    pos = mean(locs_full(1:4));
else
    pos = locs_full(1);
end

% interpolate to the common time
ttime = (0:0.01:50)';
fs = 100;
sseis = shannon(t3d, seis, ttime);
sseis2d = shannon(t2d, seis2d, ttime);
sseis3d = shannon(t3d, seis3d, ttime);

seisplot_index = [1 3 5 7 9 11];
ccplot_index = [2 4 6 8 10 12];

dt0 = datetime('now');

% plot seismograms and cc plots
for ii = 1:6
    tlim = pos + ([-0.5 0.5] + (ii-1)) * width;
    wh = and(ttime >= tlim(1), ttime < tlim(2));
    ttime_cut = ttime(wh);
    sseis_cut = sseis(wh);
    sseis2d_cut = sseis2d(wh);
    sseis3d_cut = sseis3d(wh);
    
    % SPECFEM2D vs SPECFEM3D
    [tshift1(ii), cc1(ii), lag, CC, ss1(ii)] = ccscale(sseis2d_cut, sseis3d_cut, dt0, ...
        dt0, fs, seconds(1), 'soft', false);
    
    if plt
        ax = fig2axs(seisplot_index(ii));
        plot(ax, ttime_cut, sseis2d_cut, 'LineWidth', 1)
        hold(ax);
        grid(ax);
        plot(ax, ttime_cut + tshift1(ii), sseis3d_cut * ss1(ii), 'LineWidth', 1)
        xlim(ax, tlim)
        if ii == 1
            plim = [-1.1 1.1] * max(get(ax, 'YLim'));
            set(ax, 'YLim', plim)
        else
            set(ax, 'YLim', [-1.1 1.1] * max(get(ax, 'YLim')))
        end
        if ii == 6
            xlabel(ax, 'time (s)')
            ylim(ax, plim .* [4 1]);
            legend(ax, 'SPECFEM2D', 'SPECFEM3D', 'FontSize', 9, 'Location', 'south')
        end
        ylabel(ax, 'pressure (Pa)')
        title(ax, sprintf('Cycle %d', ii));

        ax = fig2axs(ccplot_index(ii));
        plot(ax, lag, CC, 'LineWidth', 1)
        grid(ax);
        ylim(ax, [-1 1])
        if ii == 6
            xlabel(ax, 'lag (s)')
        end
        ylabel(ax, 'cc')
        title(ax, sprintf('\\tau = %.2f s, cc = %.2f, scale = %.2f', ...
            tshift1(ii), cc1(ii), ss1(ii)));
    end
    
    % FLAT vs SPECFEM3D
    [tshift2(ii), cc2(ii), lag, CC, ss2(ii)] = ccscale(sseis_cut, sseis3d_cut, dt0, ...
        dt0, fs, seconds(1), 'soft', false);
    
    if plt
        ax = fig3axs(seisplot_index(ii));
        plot(ax, ttime_cut, sseis_cut, 'LineWidth', 1)
        grid(ax);
        hold(ax);
        plot(ax, ttime_cut + tshift2(ii), sseis3d_cut * ss2(ii), 'LineWidth', 1)
        xlim(ax, tlim)
        if ii == 1
            plim = [-1.1 1.1] * max(get(ax, 'YLim'));
            set(ax, 'YLim', plim)
        else
            set(ax, 'YLim', [-1.1 1.1] * max(get(ax, 'YLim')))
        end
        if ii == 6
            xlabel(ax, 'time (s)')
            ylim(ax, plim .* [4 1]);
            legend(ax, 'Flat ocean', 'SPECFEM3D', 'FontSize', 9, 'Location', 'south')
        end
        ylabel(ax, 'pressure (Pa)')
        title(ax, sprintf('Cycle %d', ii));

        ax = fig3axs(ccplot_index(ii));
        plot(ax, lag, CC, 'LineWidth', 1)
        grid(ax);
        ylim(ax, [-1 1])
        if ii == 6
            xlabel(ax, 'lag (s)')
        end
        ylabel(ax, 'cc')
        title(ax, sprintf('\\tau = %.2f s, cc = %.2f, scale = %.2f', ...
            tshift2(ii), cc2(ii), ss2(ii)));
    end
end

if ~plt
    return
end

% resize the figure
for ii = 1:12
    if mod(ii, 2) == 0
        fig2axs(ii).Position = fig2axs(ii).Position + [0 0 0.06 0];
        fig3axs(ii).Position = fig3axs(ii).Position + [0 0 0.06 0];
    else
        fig2axs(ii).Position = fig2axs(ii).Position + [-0.03 0 0.03 0];
        fig3axs(ii).Position = fig3axs(ii).Position + [-0.03 0 0.03 0];
    end
    fig2axs(ii).Position(4) = 0.0856;
    fig3axs(ii).Position(4) = 0.0856;
end

% title
words = split(removepath(ddir3d), '_');
s2 = str2double(cindeks(split(words{2}, '-'), 2));
nu = str2double(cindeks(split(words{3}, '-'), 2));
rho = str2double(cindeks(split(words{4}, '-'), 2));
titlestring = sprintf('\\mu = %d | \\sigma^2 = %d | \\nu = %.1f | \\rho = %d', ...
    -5000, s2, nu, rho);
title(fig2axs(13), titlestring, 'FontSize', 12);
set(fig2axs(13), 'Color', 'none')
set(fig2axs(13).('XAxis'), 'Visible', 'off')
set(fig2axs(13).('YAxis'), 'Visible', 'off')
set(fig2, 'Renderer', 'painters')
figure(fig2)
savename = strcat(mfilename, '_', removepath(ddir3d), '_', '2Dvs3D', '_', useflag);
savename = replace(savename, '.', 'p');
figdisp(savename, [], [], 2, [], 'epstopdf')

title(fig3axs(13), titlestring, 'FontSize', 12);
set(fig3axs(13), 'Color', 'none')
set(fig3axs(13).('XAxis'), 'Visible', 'off')
set(fig3axs(13).('YAxis'), 'Visible', 'off')
set(fig3, 'Renderer', 'painters')
figure(fig3)
savename = strcat(mfilename, '_', removepath(ddir3d), '_', 'FLATvs3D', '_', useflag);
savename = replace(savename, '.', 'p');
figdisp(savename, [], [], 2, [], 'epstopdf')
end

function [Ypk,Xpk,Wpk,Ppk] = findpeakstopbottom(Yin,varargin)
% top peaks
[Ypk_top,Xpk_top,Wpk_top,Ppk_top] = findpeaks(Yin,varargin{:});

% bottom peaks
[Ypk_bot,Xpk_bot,Wpk_bot,Ppk_bot] = findpeaks(-Yin,varargin{:});

% merge
Ypk = [Ypk_top; -Ypk_bot];
Xpk = [Xpk_top; Xpk_bot];
Wpk = [Wpk_top; Wpk_bot];
Ppk = [Ppk_top; -Ppk_bot];

% sort by X-value
[Xpk, iXpk] = sort(Xpk);
Ypk = Ypk(iXpk);
Wpk = Wpk(iXpk);
Ppk = Ppk(iXpk);
end
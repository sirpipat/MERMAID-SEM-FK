function plotlocamp(ddir2d, ddir3d, useflag, flipflag)
% PLOTLOCAMP(ddir2d, ddir3d, useflag, flipflag)
%
% Plots seismograms from SPECFEM2D, SPECFEM3D, and ray-theory prediction
% (assuming flat ocean bottom) and then compares the locations of the peaks
% in the seismogram. This function assumes that SPECFEM2D and SPECFEM3D 
% directories and files are created by SPECFEM2D_INPUT_SETUP* or
% SPECFEM3D_INPUT_SETUP*, respectively.
%
% INPUT:
% ddir2d        directory for a SPECFEM2D run
% ddir3d        directory for a SPECFEM3D run
% useflag       either 'obs' or 'hydrophone' [default]
% flipflag      whether to flip the seismograms from SPECFEM2D run or not
%               [Default: false]
%
% SEE ALSO:
% SPECFEM2D_INPUT_SETUP, SPECFEM3D_INPUT_SETUP
%
% Last modified by sirawich-at-princeton.edu, 05/02/2025

defval('useflag', 'hydrophone')
defval('flipflag', false)

% Demo calls from the past. Now obsoleted. I keep here for example calls
% TODO: move to documentation
% if strcmp(ddir2d, 'demo1')
%     ddir2d = fullfile(getenv('REMOTE2D'), ...
%         'flat_10936816_P0009_Gaussian4Hz_new/');
%     ddir3d = fullfile(getenv('REMOTE3D'), ...
%         'FK-FLAT_10936816_P0009_Gaussian4Hz_ORIGIN_TIME/');
%     useflag = 'hydrophone';
%     flipflag = false;
%     plotlocamp(ddir2d, ddir3d, useflag, flipflag);
%     figdisp(strcat(mfilename, '_BOTTOM_4128_THETA_10-4_', useflag), [], [], 2, [], 'epstopdf')
%     return
% elseif strcmp(ddir2d, 'demo2')
%     ddir2d = fullfile(getenv('REMOTE2D'), 'FLAT_20250320', ...
%         'FLAT_BOTTOM-4000-THETA-00/');
%     ddir3d = fullfile(getenv('REMOTE3D'), 'ORIGIN_TIME_TEST_20250312', ...
%         'FK_FLAT_BAZ-225_BOTTOM-4000_THETA-00_STF-1_ORIGIN_TIME--2_NSTEPS-20000/');
%     useflag = 'hydrophone';
%     flipflag = true;
%     plotlocamp(ddir2d, ddir3d, useflag, flipflag);
%     figdisp(strcat(mfilename, '_BOTTOM_4000_THETA_00_', useflag), [], [], 2, [], 'epstopdf')
%     return
% elseif strcmp(ddir2d, 'demo3')
%     ddir2d = fullfile(getenv('REMOTE2D'), 'FLAT_20250320', ...
%         'FLAT_BOTTOM-4000-THETA-40/');
%     ddir3d = fullfile(getenv('REMOTE3D'), 'ORIGIN_TIME_TEST_20250401', ...
%         'FK_FLAT_BAZ-180_BOTTOM-4000_THETA-40_STF-2_ORIGIN_TIME-00_NSTEPS-20000/');
%     useflag = 'hydrophone';
%     flipflag = true;
%     plotlocamp(ddir2d, ddir3d, useflag, flipflag);
%     figdisp(strcat(mfilename, '_BOTTOM_4000_THETA_40_', useflag), [], [], 2, [], 'epstopdf')
%     return
% elseif strcmp(ddir2d, 'demo4')
%     ddir2d = fullfile(getenv('REMOTE2D'), 'FLAT_20250320', ...
%         'FLAT_BOTTOM-4000-THETA-10/');
%     ddir3d = fullfile(getenv('REMOTE3D'), 'ORIGIN_TIME_TEST_20250330', ...
%         'FK_FLAT_BAZ-180_BOTTOM-4000_THETA-10_STF-2_ORIGIN_TIME-00_NSTEPS-20000/');
%     useflag = 'hydrophone';
%     flipflag = true;
%     plotlocamp(ddir2d, ddir3d, useflag, flipflag);
%     figdisp(strcat(mfilename, '_BOTTOM_6000_THETA_10_', useflag), [], [], 2, [], 'epstopdf')
%     return
% elseif strcmp(ddir2d, 'demo5')
%     ddir2d = fullfile(getenv('REMOTE2D'), 'FLAT_20250320', ...
%         'FLAT_BOTTOM-4000-THETA-10/');
%     ddir3d = fullfile(getenv('REMOTE3D'), 'ORIGIN_TIME_TEST_20250312_SMALL', ...
%         'FK_FLAT_BAZ-225_BOTTOM-4000_THETA-10_STF-1_ORIGIN_TIME--6_NSTEPS-20000/');
%     useflag = 'obs';
%     flipflag = true;
%     plotlocamp(ddir2d, ddir3d, useflag, flipflag);
%     figdisp(strcat(mfilename, '_BOTTOM_4000_THETA_10_', useflag), [], [], 2, [], 'epstopdf')
%     return
% elseif strcmp(ddir2d, 'demo6')
%     ddir2d = fullfile(getenv('REMOTE2D'), 'FLAT_20250320', ...
%         'FLAT_BOTTOM-6000-THETA-00/');
%     ddir3d = fullfile(getenv('REMOTE3D'), 'ORIGIN_TIME_TEST_20250312_SMALL', ...
%         'FK_FLAT_BAZ-225_BOTTOM-6000_THETA-00_STF-1_ORIGIN_TIME--6_NSTEPS-20000/');
%     useflag = 'hydrophone';
%     flipflag = true;
%     plotlocamp(ddir2d, ddir3d, useflag, flipflag);
%     figdisp(strcat(mfilename, '_BOTTOM_6000_THETA_00_', useflag), [], [], 2, [], 'epstopdf')
%     return
% elseif strcmp(ddir2d, 'demo7')
%     ddir2d = fullfile(getenv('REMOTE2D'), 'FLAT_20250320', ...
%         'FLAT_BOTTOM-4000-THETA-00/');
%     ddir3d = fullfile(getenv('REMOTE3D'), 'ORIGIN_TIME_TEST_20250330', ...
%         'FK_FLAT_BAZ-180_BOTTOM-4000_THETA-00_STF-2_ORIGIN_TIME-00_NSTEPS-20000/');
%     useflag = 'hydrophone';
%     flipflag = true;
%     plotlocamp(ddir2d, ddir3d, useflag, flipflag);
%     figdisp(strcat(mfilename, '_BOTTOM_4000_THETA_00_', useflag), [], [], 2, [], 'epstopdf')
%     return
% end

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

% flip SPECFEM2D seismograms when STF for FKSPECFEM3D is upside-down
% Gaussian
if flipflag
    zobs2d = -zobs2d;
    pmh2d = -pmh2d;
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
    while ii <= length(locs2d_rel) && abs(locs2d_rel(ii) - locs_rel(ii)) > 0.5
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

%% plot
figure(1)
clf
set(gcf, 'Units', 'inches', 'Position', [0 1 6 8])

ax1 = subplot('Position', [0.09 0.82 0.84 0.12]);
plot(t2d, seis2d, 'LineWidth', 1, 'Color', [0.65 0.65 0.65])
%plot(t2d + (10-locs2d(1)), seis2d, 'LineWidth', 1)
grid on
hold on
scatter(locs2d, pks2d)
%scatter(locs2d + (10-locs2d(1)), pks2d)
xlim([0 50])
set(gca, 'TickDir', 'out', 'FontSize', 12)
nolabels(gca, 1)
title('SPECFEM2D')

ax2 = subplot('Position', [0.09 0.64 0.84 0.12]);
plot(t3d, seis3d, 'LineWidth', 1, 'Color', 'k')
%plot(t3d + (10-locs3d(1)), seis3d, 'LineWidth', 1)
grid on
hold on
scatter(locs3d, pks3d)
%scatter(locs3d + (10-locs3d(1)), pks3d)
xlim([0 50])
set(gca, 'TickDir', 'out', 'FontSize', 12)
nolabels(gca, 1)
title('SPECFEM3D')

ax3 = subplot('Position', [0.09 0.46 0.84 0.12]);
plot(t3d, seis, 'LineWidth', 1, 'Color', [0.65 0.65 0.65])
%plot(t3d + (10-locs(1)), seis, 'LineWidth', 1)
grid on
hold on
scatter(locs, pks)
%scatter(locs + (10-locs(1)), pks)
xlabel('time (s)')
xlim([0 50])
set(gca, 'TickDir', 'out', 'FontSize', 12)
title('Flat ocean bottom')

% set the ylimit of the 3 plots to be the same
ylim1 = get(ax1, 'YLim');
ylim2 = get(ax2, 'YLim');
ylim3 = get(ax3, 'YLim');

ylim_all = [min([ylim1(1) ylim2(1) ylim3(1)]) ...
    max([ylim1(2) ylim2(2) ylim3(2)])];
% ylim_all = [-1 1] * max(abs(ylim_all));

set(ax1, 'YLim', ylim_all)
set(ax2, 'YLim', ylim_all)
set(ax3, 'YLim', ylim_all)

%% rescale and shift peaks
%locs2d = locs2d - locs2d(1);
%locs3d = locs3d - locs3d(1);
pks2d = pks2d / pks2d(1);
pks3d = pks3d / pks3d(1);
%locs = locs - locs(1);
pks = pks / pks(1);

% find the common peaks between simulations and predictions
keeplength = min(min(length(locs3d), length(locs2d)), length(locs));
locs = locs(1:keeplength);
pks = pks(1:keeplength);
locs2d = locs2d(1:keeplength);
pks2d = pks2d(1:keeplength);
locs3d = locs3d(1:keeplength);
pks3d = pks3d(1:keeplength);

subplot('Position', [0.09 0.07 0.36 0.25])
scatter(locs, locs2d, 'filled')
hold on
scatter(locs, locs3d, 'filled')
grid on
XLIM = get(gca, 'XLim');
YLIM = get(gca, 'YLim');
set(gca, 'XLim', [min(XLIM(1), YLIM(1)) max(XLIM(2), YLIM(2))], ...
    'YLim', [min(XLIM(1), YLIM(1)) max(XLIM(2), YLIM(2))])
hold on
refline(gca, 1, 0);
legend('SPECFEM2D', 'SPECFEM3D', 'Location', 'best')
set(gca, 'TickDir', 'out', 'FontSize', 12, 'Box', 'on')
ylabel('SPECFEM (s)')
xlabel('Flat ocean bottom (s)')
title('Arrival times')


subplot('Position', [0.59 0.07 0.36 0.25])
scatter(locs, locs2d - locs, 'filled')
hold on
scatter(locs, locs3d - locs, 'filled')
grid on
legend('SPECFEM2D', 'SPECFEM3D', 'Location', 'best')
set(gca, 'TickDir', 'out', 'FontSize', 12, 'Box', 'on')
ylabel('arrival time delay (s)')
xlabel('Flat ocean bottom (s)')
title('Time delays')

set(gcf, 'Renderer', 'painters')
savename = strcat(mfilename, '_', removepath(ddir3d), '_', useflag);
savename = replace(savename, '.', 'p');
figdisp(savename, [], [], 2, [], 'epstopdf')
end
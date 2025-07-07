function stfmaker(obs_struct, synmasterdir, specmasterdir)
% STFMAKER(obs_struct, synmasterdir, specmasterdir)
%
% Prepares stf_file for FK-SPECFEM3D runs from Instaseis seimogram. This
% will override the existing stf_file.txt in the directories.
%
% INPUT:
% obs_struct        a struct containing
%   - snr               signal-to-noise ratio
%   - fcorners          [lower upper] corner frequencies
%   - CCmaxs            maximum correlation coefficients for
%                       [flat bath] cases
%   - metadata          SAC Headers associated to the obsfile
%   - presiduals        InstaSeis arrival - TauP prediction for first P
%                       arrival
% synmasterdir      the master directory to the synthetic files sorted into
%                   IRIS event ID folders
% specmasterdir     the master directory to the input folders from
%                   FK-SPECFEM3D run for fluid-solid setting. Must have the
%                   following subdirectories
%   - LAYERED_OC_MERMAID_HIGHCC_**     ** is 01 to 15
%   - LAYERED_OC_MERMAID_LOWCC_**      ** is 01 to 15
% 
% Last modified by sirawich-at-princeton.edu, 07/07/2025

deval('specmasterdir', ...
    fullfile(getenv('REMOTE3D'), '20250629_MERMAID_INSTASEIS'))

% sort the observations by decreasing correlation coefficient
[~, ic] = sort(obs_struct.CCmaxs(:,2), 'descend');
% loop through the runs
for ii = 1:15
    % identify which Instaseis syntheteic seismogram to read
    stationid = indeks(obs_struct.metadata.KSTNM{ic(ii)}, '2:5');
    eventid = obs_struct.metadata.USER7(ic(ii));
    synfile = cindeks(ls2cell(sprintf('%s%d/*_%s_*.sac', synmasterdir, ...
        eventid, stationid), 1), 1);
    
    % read the Instaseis synthetic vertical displacement seismogram
    [seis_s, hdr_s, ~, ~, tims_s] = readsac(synfile);
    tims_s = tims_s';

    % adjust the time based on the ray-theory arrival time of first phase
    t = tims_s - hdr_s.T0;

    % compute the tapering window for FK-SPECFEM3D injection
    seis_s = seis_s .* stf_window(t);

    % write the source-time function from 200 s before to 200 s after the
    % ray-theory arrival time
    wh = and(t >= -200, t<= 200);
    ddir = fullfile(specmasterdir, ...
        sprintf('LAYERED_OC_MERMAID_HIGHCC_%02d', ii));
    writetimeseries(t(wh), seis_s(wh), fullfile(ddir, 'stf_file.txt'));
end

% sort the observations by increasing correlation coefficient
[~, ic] = sort(obs_struct.CCmaxs(:,2), 'ascend');
for ii = 1:15
    % identify which Instaseis syntheteic seismogram to read
    stationid = indeks(obs_struct.metadata.KSTNM{ic(ii)}, '2:5');
    eventid = obs_struct.metadata.USER7(ic(ii));
    synfile = cindeks(ls2cell(sprintf('%s%d/*_%s_*.sac', synmasterdir, ...
        eventid, stationid), 1), 1);
    
    % read the Instaseis synthetic vertical displacement seismogram
    [seis_s, hdr_s, ~, ~, tims_s] = readsac(synfile);
    tims_s = tims_s';

    % adjust the time based on the ray-theory arrival time of first phase
    t = tims_s - hdr_s.T0 - hdr_s.B;

    % compute the tapering window for FK-SPECFEM3D injection
    seis_s = seis_s .* stf_window(t);

    % write the source-time function from 200 s before to 200 s after the
    % ray-theory arrival time
    wh = and(t >= -200, t<= 200);
    ddir = fullfile(specmasterdir, ...
        sprintf('LAYERED_OC_MERMAID_LOWCC_%02d', ii));
    writetimeseries(t(wh), seis_s(wh), fullfile(ddir, 'stf_file.txt'));
end
end
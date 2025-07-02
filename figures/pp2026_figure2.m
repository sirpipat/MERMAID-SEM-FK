function pp2026_figure2(obs_struct, obsmasterdir, synmasterdir, presmasterdir)
% PP2026_FIGURE2(obs_struct, obsmasterdir, synmasterdir, presmasterdir)
%
% Makes figure 2A and 2B for Pipatprathanporn+2026 paper.
%
% Last modified by sirawich-at-princeton.edu, 07/02/2025

keywords = {'ascend', 'descend'};
keyfolders = {'LOWCC', 'HIGHCC'};
ii_cases = 11:15;
for kk = 1:2
    [~, ic] = sort(obs_struct.CCmaxs(:,2), keywords{kk});
    
    figure;
    set(gcf, 'Units', 'inches', 'Position', [0 1 12 7])
    
    xposshift = [0 0.2 0.4 0.6 0.8];
    
    for ii = 1:5
        %% Top row: bird-eye-view bathymetry plot
        [ll, tt, zz, llons, llats] = ...
            bathymetryprofile2d([20000 20000], [81 81], ...
            [obs_struct.metadata.STLO(ic(ii_cases(ii))) ...
            obs_struct.metadata.STLA(ic(ii_cases(ii)))], ...
            obs_struct.metadata.BAZ(ic(ii_cases(ii)))+180);
    
        subplot('Position', [0.030+xposshift(ii) 0.65 0.16 0.25])
        imagesc([-10 10], [-10 10], rot90(zz'));
        colormap(kelicol);
        clim(mean(zz(:) + std(zz(:)) * [-3 3]))
        % cb = colorbar('Location', 'southoutside');
        % cb.Label.String = 'elevation (m)';
        % cb.TickDirection = 'out';
        grid on;
        box on;
        hold on;
        plot([-10 10], [0 0], 'LineWidth', 2, 'Color', 'k')
        axis image;
        axis xy;
        xticks(-10:5:10)
        yticks(-10:5:10)
        if ii > 1
            nolabels(gca, 3);  % only keep y-axis label for leftmost plot
        else
            nolabels(gca, 1);
            yticklabels(-10:5:10)
            ylabel('transverse position (km)')
        end
        set(gca, 'TickDir', 'out', 'FontSize', 9, ...
            'Position', [0.033+xposshift(ii) 0.68 0.16 0.25])
        title(sprintf('%d-%s', obs_struct.metadata.USER7(ic(ii_cases(ii))), ...
            obs_struct.metadata.KSTNM{ic(ii_cases(ii))}))
    
        %% Middle row: bathymetry profile along great-circle path
        subplot('Position', [0.038+xposshift(ii) 0.48 0.15 0.18])
        hold on
        for jj = 1:size(zz,1)
            hold on
            if jj ~= 41
                plot(ll(1,:)/1000, zz(jj,:)/1000, 'LineWidth', 0.25, 'Color', ...
                    [0.7 0.7 0.7]);
            end
        end
        plot(ll(1,:)/1000, zz(41,:)/1000, 'LineWidth', 1, 'Color', 'k');
        grid on
        box on
        xticks(-10:5:10)
        xticklabels(-10:5:10)
        xlabel('radial position (km)')
        ylabel('elevation (km)')
        
        %% Bottom row: seismograms
        subplot('Position', [0.038+xposshift(ii) 0.07 0.15 0.26])
    
        % identify observed pressure file
        eventid = obs_struct.metadata.USER7(ic(ii_cases(ii)));
        stationid = indeks(obs_struct.metadata.KSTNM{ic(ii_cases(ii))}, '2:5');
        try
            obsfile = cindeks(ls2cell(sprintf('%s%d/*.%s_*.sac', ...
                obsmasterdir, eventid, stationid), 1), 1);
        catch ME
            if strcmp(ME.message, 'This directory or file does not exist')
                obsfile = cindeks(ls2cell(sprintf('%s%d/*.%s_*.sac', ...
                    obsmasterdir, eventid, stationid(end-1:end)), 1), 1);
            end
        end
    
        % read observed pressure file
        [seis_o, hdr_o] = readsac(obsfile);
        [dt_ref_o, ~, ~, fs_o, ~, dts_o] = gethdrinfo(hdr_o);
        t_relative = seconds(dts_o - dt_ref_o) - hdr_o.T0;
        
        % remove instrument response
        seis_o = counts2pa(seis_o, fs_o, [0.01 0.02 5 10], [], [], false);

        % corner frequency
        fc = obs_struct.fcorners(ic(ii_cases(ii)), :);

        % badpass filtering
        seis_o2 = bandpass(detrend(seis_o), fs_o, fc(1), fc(2), 4, 2, ...
            'butter', 'linear');
        
        % identify synthetic pressure file
        stationid = indeks(obs_struct.metadata.KSTNM{ic(ii_cases(ii))}, '2:5');
        presfile = cindeks(ls2cell(sprintf('%s/%d/*_%s_0_*.sac', ...
            presmasterdir, eventid, stationid), 1), 1);
        
        % read synthetic pressure file
        [seis_p, hdr_p] = readsac(presfile);
        [dt_ref_p, ~, ~, fs_p, ~, dts_p] = gethdrinfo(hdr_p);

        % bandpass filtering
        seis_p2 = bandpass(detrend(seis_p), fs_p, fc(1), fc(2), 4, 2, ...
            'butter', 'linear');

        % normalize
        t_shift = obs_struct.t_shifts(ic(ii_cases(ii)), 2);
        wh_s = and(t_relative + t_shift >= -10, t_relative + t_shift <= 30);
        seis_p2 = seis_p2 / max(abs(seis_p2(wh_s)));
    
        % normalize
        wh = and(t_relative >= -10, t_relative <= 30);
        plot(t_relative, seis_o2 / max(abs(seis_o2(wh))), 'LineWidth', 1, 'Color', 'k')
        hold on
        plot(t_relative + t_shift, seis_p2 + 2.5, 'LineWidth', 0.75, 'Color', 'r')
        xlim([-10 30])
        grid on
        hold on
        set(gca, 'FontSize', 9, 'TickDir', 'out')

        ddir = fullfile(getenv('REMOTE3D'), '20250629_MERMAID_INSTASEIS', ...
            sprintf('LAYERED_OC_MERMAID_%s_%02d', keyfolders{kk}, ...
            ii_cases(ii)));
        % read modeled acoustic pressure seismograms from SPECFEM3D
        [t,x] = read_seismogram(fullfile(ddir, 'OUTPUT_FILES/MH.P0009.CXP.semp'));
        % read the locations of the stations
        [~, ~, ~, ~, ~, z] = readstations3d(fullfile(ddir, 'DATA', 'STATIONS'));
        % Calculate the time of the simulation relative to estimated
        % P-wave arrival time from Simon et al. (2022) by
        % first calculate the time in the SPECFEM3D simulation when
        % P-wave arrive the ocean bottom below the float for the first
        % time
        t0 = calculatearrivaltime(ddir);
        tTauP = indeks(tauptime('mod', 'ak135', ...
            'dep', obs_struct.metadata.EVDP(ic(ii_cases(ii))), ...
            'ph', 'p,P,PKP,PKIKP', ...
            'deg', obs_struct.metadata.GCARC(ic(ii_cases(ii))), ...
            'stdp', -z(1)/1000), 1).time;
        tSPECFEM = t - t0 + (hdr_o.USER8 - hdr_o.T0) + tTauP;
        fs_SPECFEM = (length(t) - 1) / (t(end) - t(1));

        % interpolate the seismogram to align samples with the observed
        % seismogram
        x = lowpass(detrend(x) .* shanning(length(x), 0.2), fs_SPECFEM, 10, 2, 2, 'butter', 'linear');
        x = shannon(tSPECFEM, x, t_relative(wh));

        % apply bandpass filter
        xf = bandpass(detrend(x) .* shanning(length(x), 0.2), fs_o, fc(1), fc(2), 4, 2, ...
            'butter', 'linear');

        % read z-displacement at the ocean bottom from SPECFEM3D
        [tb,xb] = read_seismogram(fullfile(ddir, 'OUTPUT_FILES/AA.OBS01.CXZ.semd'));

        % downsample down to 20 hz
        xb = lowpass(detrend(xb) .* shanning(length(xb), 0.2), fs_SPECFEM, 10, 2, 2, 'butter', 'linear');
        xb = downsample(xb, floor(fs_SPECFEM/20));
        tb = downsample(tb, floor(fs_SPECFEM/20));

        % apply bandpass filter
        xbf = bandpass(detrend(xb) .* shanning(length(xb), 0.2), fs_SPECFEM / floor(fs_SPECFEM/20), ...
            fc(1), fc(2), 4, 2, 'butter', 'linear');

        % read Instaseis z-displacement at the ocean bottom
        synfile = cindeks(ls2cell(sprintf('%s%d/*_%s_*.sac', synmasterdir, eventid, stationid), 1), 1);
        [seis_s, hdr_s, ~, ~, tims_s] = readsac(synfile);
        [~, ~, ~, fs_s] = gethdrinfo(hdr_s);
        
        % apply bandpass filter
        seis_sf = bandpass(detrend(seis_s), fs_s, fc(1), fc(2), 4, ...
            2, 'butter', 'linear');
        
        % amplitude scaling for SPECFEM seismogram
        amp_xbf = max(abs(xbf(and(tb >= t0-2, tb <= t0+2))));
        amp_seis_sf = max(abs(seis_sf(and(tims_s >= hdr_s.T0 - 10, ...
            tims_s <= hdr_s.T0 + 30))));

        % rescale
        xf = xf * amp_seis_sf / amp_xbf;

        % compute the cross-correlation of the envelope
        seis_o3 = seis_o2(wh);
        [t_shift1, CCmax1, lags1, cc1, s1] = ccscale(seis_o3, xf, datetime('today'), datetime('today'), fs_o, seconds(10), 'hard', true);
        [t_shift2, CCmax2, lags2, cc2, s2] = ccscale(seis_o3, xf, datetime('today'), datetime('today') + seconds(t_shift1), fs_o, seconds(5), 'soft', false);
        t_shift3D = t_shift1 + t_shift2;

        plot(t_relative(wh) + t_shift3D, xf / max(abs(xf)) - 2.5, 'LineWidth', 0.75, 'Color', [0 0.2 0.9])
        title(sprintf('cc: %.2f | \\Delta\\tau: %.2f s | s: %.2g', CCmax2, t_shift3D, s2))

        xlabel('time since picked P-wave arrival (s)')
        legend('observed', ...
            sprintf('2D (cc: %.2f, \\Delta\\tau: %.2f s)', ...
            obs_struct.CCmaxs(ic(ii_cases(ii)), 2), ...
            obs_struct.t_shifts(ic(ii_cases(ii)), 2)), ...
            sprintf('3D (cc: %.2f, \\Delta\\tau: %.2f s)', ...
            CCmax2, t_shift3D), 'Location', 'southoutside')
        yticks([-2.5 0 2.5]);
        nolabels(gca, 2);
        set(gca, 'Position', [0.038+xposshift(ii)  0.14 0.15 0.26])
    end
    set(gcf, 'Renderer', 'painters')
    figdisp(sprintf('%s_%s_%02d-%02d', mfilename, keyfolders{kk}, ii_cases(1), ii_cases(end)), [], [], 2, [], 'epstopdf');
end
end
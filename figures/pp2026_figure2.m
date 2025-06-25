function pp2026_figure2(obs_struct, obsmasterdir, presmasterdir)
% PP2026_FIGURE2
%
% Makes figure 2A and 2B for Pipatprathanporn+2026 paper.
%
% Last modified by sirawich-at-princeton.edu, 06/19/2025

keywords = {'ascend', 'descend'};
for kk = 1:2
    [~, ic] = sort(obs_struct.CCmaxs(:,2), keywords{kk});
    
    figure;
    set(gcf, 'Units', 'inches', 'Position', [0 1 12 7])
    
    xposshift = [0 0.2 0.4 0.6 0.8];
    
    for ii = 1:5
        [ll, tt, zz, llons, llats] = ...
            bathymetryprofile2d([20000 20000], [81 81], ...
            [obs_struct.metadata.STLO(ic(ii)) ...
            obs_struct.metadata.STLA(ic(ii))], ...
            obs_struct.metadata.BAZ(ic(ii))+180);
    
        subplot('Position', [0.040+xposshift(ii) 0.65 0.16 0.25])
        imagesc([-10000 10000], [-10000 10000], rot90(zz'));
        colormap(kelicol);
        cb = colorbar('Location', 'southoutside');
        cb.Label.String = 'elevation (m)';
        cb.TickDirection = 'out';
        grid on;
        box on;
        hold on;
        plot([-10000 10000], [0 0], 'LineWidth', 2, 'Color', 'k')
        axis image;
        axis xy;
        xticks(-10000:5000:10000)
        xticklabels(-10000:5000:10000)
        yticks(-10000:5000:10000)
        yticklabels(-10000:5000:10000)
        xlabel('radial position (m)')
        ylabel('traverse position (m)')
        set(gca, 'TickDir', 'out', 'FontSize', 9, ...
            'Position', [0.043+xposshift(ii) 0.68 0.16 0.25])
        title(sprintf('%d-%s | CC: %.2f', obs_struct.metadata.USER7(ic(ii)), ...
            obs_struct.metadata.KSTNM{ic(ii)}, ...
            obs_struct.CCmaxs(ic(ii), 2)))
    
        subplot('Position', [0.048+xposshift(ii) 0.32 0.15 0.18])
        hold on
        for jj = 1:size(zz,1)
            hold on
            if jj ~= 41
                plot(ll(1,:), zz(jj,:), 'LineWidth', 0.25, 'Color', ...
                    [0.7 0.7 0.7]);
            end
        end
        plot(ll(1,:), zz(41,:), 'LineWidth', 1, 'Color', 'k');
        grid on
        box on
        xticks(-10000:5000:10000)
        xticklabels(-10000:5000:10000)
        xlabel('radial position (m)')
        ylabel('elevation (m)')
        
        subplot('Position', [0.048+xposshift(ii) 0.07 0.15 0.16])
    
        % identify observed pressure file
        eventid = obs_struct.metadata.USER7(ic(ii));
        stationid = indeks(obs_struct.metadata.KSTNM{ic(ii)}, '2:5');
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
        t_shift = obs_struct.t_shifts(ic(ii), 2);
        wh = and(t_relative + t_shift >= -10, t_relative + t_shift <= 30);
        seis_o2 = bandpass(seis_o, fs_o, obs_struct.fcorners(ic(ii),1), ...
            obs_struct.fcorners(ic(ii), 2), 4, 2, 'butter', 'linear');
        seis_o2 = seis_o2 / max(abs(seis_o2(wh)));
    
        % identify synthetic pressure file
        stationid = indeks(obs_struct.metadata.KSTNM{ic(ii)}, '2:5');
        presfile = cindeks(ls2cell(sprintf('%s/%d/*_%s_0_*.sac', ...
            presmasterdir, eventid, stationid), 1), 1);
        
        % read synthetic pressure file
        [seis_p, hdr_p] = readsac(presfile);
        [dt_ref_p, ~, ~, fs_p, ~, dts_p] = gethdrinfo(hdr_p);
        seis_p2 = bandpass(seis_p, fs_p, obs_struct.fcorners(ic(ii),1), ...
            obs_struct.fcorners(ic(ii), 2), 4, 2, 'butter', 'linear');
        seis_p2 = seis_p2 / max(abs(seis_p2(wh)));
    
        plot(t_relative, seis_o2, 'LineWidth', 1, 'Color', 'k')
        hold on
        plot(t_relative + t_shift, seis_p2, 'LineWidth', 1, 'Color', 'r')
        xlim([-10 30])
        grid on
        hold on
        set(gca, 'FontSize', 9, 'TickDir', 'out')
    end
end
end
%%
out_dir = 'results';
sub_chs = cell2table({ ...
    '01',   '01_20210825-2',    'L1', 'R10'; ... 
    '02',   '02_20210529-1',    'Q1', 'N9'; ...
    '03',   '03_20240525-3',    'B1', 'D4'; ... 
    '04',   'DR7031X9',         'T''1', 'S''7'; ... 
    '05',   '05_20210509',      'A2', 'D14'; ...
    '06',   '06_20200610',      'E1', 'F3'; ...
    '07',   '07_20220215',      'C1', 'D8'; ...
    '08',   '08_20241018-7',    'X2', 'B2'; ...
    '09',   '09_20220705',      'J16', 'C5'; ...
    }, 'VariableNames', {'Subject', 'Session', 'ANT', 'SOZ'});

%%
max_lag = 0.5;
max_nlag = round(max_lag * fs);
t_lag = (-max_nlag:max_nlag)' / fs;
t_lag_high = (-0.05:0.001:0.05)';
h = design(fdesign.bandpass('N,F3dB1,F3dB2', 4, 3, 100, fs), 'butter');
for nn = 1:numel(ds)
    t0 = tic; 
    fn = sprintf('%s_%s', ds(nn).subject, ds(nn).fname);    
    subid = strcmp(sub_chs.Subject, ds(nn).subject);

    load(fullfile(out_dir, 'prep_sig', [fn, '.mat']), 'sig_ant', 'sig_zone');
    sig_ant = sig_ant(:, strcmp(ds(nn).elec_ant, sub_chs.ANT{subid}));
    sig_ant = filtfilt(h.sosMatrix, h.ScaleValues, sig_ant);
    sig_soz = sig_zone(:, strcmp(ds(nn).elec_zone, sub_chs.SOZ{subid}));
    sig_soz = filtfilt(h.sosMatrix, h.ScaleValues, sig_soz);
    ts = (0:ds(nn).nt-1)' / fs;
    
    ch1 = strcmp(ds(nn).elec_ant(ds(nn).is_ant), sub_chs.ANT{subid});
    ch2 = strcmp(ds(nn).elec_zone(ds(nn).is_zone), sub_chs.SOZ{subid});
    % tev = evToVec(max(round(ds(nn).ant_ev_ied{ch1} * fs), 1), [ds(nn).nt, 1]) | ...
    %     evToVec(max(round(ds(nn).soz_ev_ied{ch2} * fs), 1), [ds(nn).nt, 1]);
    % tev = [find(diff([false; tev]) == 1), find(diff([tev; false]) == -1)] / fs;
    tev = ds(nn).ant_ev_ied{ch1};

    nev = size(tev, 1);
    V = zeros(nev, 2);
    C = zeros(nev, numel(t_lag_high));
    for kk = 1:nev
        tin = ts >= tev(kk, 1) & ts <= tev(kk, 2);
        cc = xcov(sig_ant(tin), sig_soz(tin), round(max_lag * fs), 'coeff');
        cc = interp1(t_lag, cc, t_lag_high, 'spline');
        C(kk, :) = cc;
        [vmax, imax] = max(cc);
        V(kk, :) = [vmax, t_lag_high(imax)];
    end
    ds(nn).tlag_ied = V;
    ds(nn).ied_xc = C;
    ds(nn).ev_ied = tev;

    figure; set(gcf, 'Position', [50, 50, 800, 300], 'Visible', 'off');
    subplot(1, 2, 1);
    histogram(V(:, 2), t_lag_high);
    subplot(1, 2, 2);
    plot(V(:, 2), V(:, 1), '.');
    set(gca, 'XLim', minmax(t_lag_high(:)'));
    sgtitle(fn, 'Interpreter', 'none');
    print(gcf, fullfile(out_dir, 'check_ied_tlag', fn), '-djpeg');
    close;
end
%%
nsub = size(sub_chs, 1);
stats = struct;
for ii = 1:nsub
    stats(ii, 1).subject = sub_chs.Subject{ii};    
    ids = arrayfun(@(x)(strcmp(x.subject, sub_chs.Subject{ii}) & ...
        strcmp(x.fname, sub_chs.Session{ii})), ds);
    V = cat(1, ds(ids).tlag_ied);
    tev = cat(1, ds(ids).ev_ied);
    dur = diff(tev, 1, 2);
    rt = sum(cat(1, ds(ids).duration));
    
    stats(ii).record_hr = rt / 60 / 60;
    stats(ii).num_events = size(V, 1);
    stats(ii).median_duration = median(dur);
    stats(ii).total_time_hr = sum(dur) / stats(ii).record_hr;
    stats(ii).times_per_hr = size(tev, 1) / stats(ii).record_hr;
    stats(ii).num_ant_only = sum(V(:, 1) <= 0.4);
    stats(ii).num_ant_syn = sum(V(:, 1) > 0.4);
    stats(ii).num_ant_early = sum(V(:, 2) < 0 & V(:, 1) > 0.4);
    stats(ii).num_ant_late = sum(V(:, 2) >= 0 & V(:, 1) > 0.4);
    stats(ii).perc_ant_only = stats(ii).num_ant_only / stats(ii).num_events;
    stats(ii).perc_ant_syn = stats(ii).num_ant_syn / stats(ii).num_events;
    stats(ii).perc_ant_early = stats(ii).num_ant_early / stats(ii).num_events;
    stats(ii).perc_ant_late = stats(ii).num_ant_late / stats(ii).num_events;
    stats(ii).perc_ant_syn_early = stats(ii).num_ant_early / stats(ii).num_ant_syn;
    stats(ii).perc_ant_syn_late = stats(ii).num_ant_late / stats(ii).num_ant_syn;
end
%%
fig_name = 'seeg_ied_ant_syn_stats_onefile';
cmap = [ ...
    1, 0.4, 0.4; ... early
    0.2, 0.5, 1; ... late
    .7, .7, .7; ... only
    ];
sub_order = [1, 6, 2, 4, 3, 5, 9, 7, 8];

figure; set(gcf, 'Position', [50, 50, 1000, 780]);

for ii = 1:nsub
    subplot(3, 5, ii+5);
    p = pie([stats(ii).num_ant_early, stats(ii).num_ant_late, stats(ii).num_ant_only]);
    colormap(cmap);
    for cc = 1:3
        p(cc*2-1).LineWidth = 2;
        p(cc*2-1).LineStyle = 'none';
        % p(cc*2-1).EdgeColor = 'w';
        vs = p(cc*2-1).Vertices(2:end-1, :);
        [x, y] = pol2cart(mean(unwrap(cart2pol(vs(:, 1), vs(:, 2)))), 0.6);
        p(cc*2).Position = [x, y, 0];
        p(cc*2).FontSize = 10;
        p(cc*2).HorizontalAlignment = 'center';
        p(cc*2).VerticalAlignment = 'middle';
    end
    title(stats(ii).subject, 'interpreter', 'none');
end
legend({'ANT-early', 'ANT-late', 'ANT-independent'}, 'Box', 'off', ...
    'Position', [0.82, 0.18, 0.1, 0.1], 'FontSize', 10);

subplot(3, 5, [1, 2]);
vmat = cat(2, cat(1, stats(:).num_ant_early), cat(1, stats(:).num_ant_late), ...
    cat(1, stats(:).num_ant_only)) / stats(ii).record_hr;
h = bar(1:nsub, vmat(sub_order, :)', 'grouped');
for cc = 1:3
    h(cc).FaceColor = cmap(cc, :);
    h(cc).EdgeColor = 'w';
end
set(gca, 'TickDir', 'out', 'XTickLabel', strrep(sub_chs.Subject(sub_order), '_', ''), ...
    'XTickLabelRotation', 30);
ylabel('Number (/hr)');

subplot(3, 5, [4, 5]);
vmat = cat(2, cat(1, stats(:).perc_ant_early), ...
    cat(1, stats(:).perc_ant_late), cat(1, stats(:).perc_ant_only)) * 100;
h = bar(1:nsub, vmat(sub_order, :)', 'stacked', 'BarWidth', 0.6);
for cc = 1:3
    h(cc).FaceColor = cmap(cc, :);
    h(cc).EdgeColor = 'w';
end
set(gca, 'TickDir', 'out', 'XTickLabel', strrep(sub_chs.Subject(sub_order), '_', ''), ...
    'XTickLabelRotation', 30);
ylabel('Percentage (%)');

print(gcf, fullfile(out_dir, fig_name), '-djpeg');
print(gcf, fullfile(out_dir, fig_name), '-dpdf', '-r300', '-vector');









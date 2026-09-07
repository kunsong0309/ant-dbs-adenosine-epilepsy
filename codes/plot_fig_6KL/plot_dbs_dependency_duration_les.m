%%
raw_dir = '..\WT_ANT_GFAP_Cas3_GAdo\';

dbs_list = { ...
    'WT_ANT_GFAP_Cas3_GAdo_M1_DBS3_260202', [10, 60, 300, 600, 900]; ...
    'WT_ANT_GFAP_Cas3_GAdo_M2_DBS3_260202', [10, 60, 300, 600, 900]; ...
    'WT_ANT_GFAP_Cas3_GAdo_M4_DBS3_260203', [10, 60, 300, 600, 900]; ...
    'WT_ANT_GFAP_Cas3_GAdo_M7_DBS3_260203', [10, 60, 300, 600, 900]; ...
    'WT_ANT_GFAP_Cas3_GAdo_M9_DBS3_260202', [60, 300, 600, 900, 10]; ...
    'WT_ANT_GFAP_Cas3_GAdo_M10_DBS3_260204', [10, 60, 300, 600, 900]; ...
    'WT_ANT_GFAP_Cas3_GAdo_M1_DBS4_260206', [10, 60, 300, 600, 900]; ...
    'WT_ANT_GFAP_Cas3_GAdo_M2_DBS4_260206', [10, 60, 300, 600, 900]; ...
    'WT_ANT_GFAP_Cas3_GAdo_M4_DBS4_260206', [10, 60, 300, 600, 900]; ...
    'WT_ANT_GFAP_Cas3_GAdo_M7_DBS4_260209', [10, 60, 300, 600, 900]; ...
    };

%%
ds = struct;
for nn = 1:size(dbs_list, 1)
    fls = dir(fullfile(raw_dir, dbs_list{nn, 1}, '*.csv'));
    fname = dbs_list{nn, 1};
    fpath = fullfile(fls(1).folder, fls(1).name);
    tbl = readtable(fpath, 'VariableNamingRule', 'preserve');

    ts0 = seconds(tbl.Time - tbl.Time(1));
    sig0 = [tbl.('ROI-1_sig_470'), tbl.('ROI-1_sig_410')];
    iso = any(isnan(sig0), 2);
    ts0(iso) = [];
    sig0(iso, :) = [];
    mrk = ts0(find(diff(tbl.('Input_DIO-2')) == 1) + 1);

    fs = round(1/median(diff(ts0)));
    ts = (0:(1/fs):ts0(end))';
    sig0 = interp1(ts0, sig0, ts, 'linear');

    sig0 = movmean(sig0, round(fs));
    sig0 = sig0(1:fs:end, :);
    ts = ts(1:fs:end);

    dsig = detrend(sig0, 2);
    sig1 = dsig ./ (sig0 - dsig) * 100;

    sig = sig1(:, 1) / (sqrt(diff(prctile(sig1(:, 1), [1, 99]))) * 0.1);
    % x = sig1 .* [0, 1] + [1, 0];
    % b = x\sig1(:, 1);
    % sig = sig1(:, 1) - x * b;

    ds(nn, 1).fname = fname;
    ds(nn).fpath = fpath;
    ds(nn).ts = ts;
    ds(nn).sig = sig;
    ds(nn).mrk = mrk;
    ds(nn).amp = dbs_list{nn, 2};

    % figure; set(gcf, 'Position', [50, 300, 1200, 600]); 
    % subplot(2, 1, 1); hold on;
    % plot(ts, sig0);
    % plot(ts, sig0 - dsig);
    % subplot(2, 1, 2);
    % plot(ts, sig);
    % xline(mrk, ':k');
    % title(fname, 'Interpreter', 'none');
    % print(gcf, fname, '-djpeg'); close;
end

%%
dur_list = [60, 300, 600, 900];
ndur = length(dur_list);
tx = (-30:0.5:1100);
for nn = 1:numel(ds)
    sig = ds(nn).sig;
    ts = ds(nn).ts;    

    [traces, peaks, aucs] = deal(cell(ndur, 1));
    for kk = 1:ndur
        isin = ds(nn).amp==dur_list(kk);
        mrk = ds(nn).mrk(isin);

        evt = mrk + [-30, 0];
        [bsl, ~] = timestampsInterpEvents(ds(nn).sig, ds(nn).ts, evt, 'avgrg');

        evt = mrk + tx;
        [omat, ~] = timestampsInterpEvents(ds(nn).sig, ds(nn).ts, evt, 'trace');
        traces{kk} = mean(omat - bsl, ndims(omat));

        evt = mrk + tx(tx >= 0 & tx <= (dur_list(kk) + 0));
        [omat, ~] = timestampsInterpEvents(ds(nn).sig, ds(nn).ts, evt, 'trace');
        aucs{kk} = mean(mean(omat - bsl) * dur_list(kk), ndims(omat));

        evt = mrk + [0, dur_list(kk) + 0];
        [omat, ~] = timestampsInterpEvents(ds(nn).sig, ds(nn).ts, evt, 'max');
        peaks{kk} = mean(omat - bsl, ndims(omat));
    end
    ds(nn).traces = cat(3, traces{:});
    ds(nn).peaks = cat(3, peaks{:});
    ds(nn).aucs = cat(3, aucs{:});
end
traces = cat(2, ds(:).traces);
peaks = cat(2, ds(:).peaks);
aucs = cat(2, ds(:).aucs);

%%
cmap = tab10;
dur_label = arrayfun(@(x)(sprintf('%d sec', x)), dur_list, 'UniformOutput', false);
figure; set(gcf, 'Position', [50, 50, 780, 780]);

subplot('Position', [0.1, 0.5, 0.25, 0.25]); hold on;
hs = nan(kk, 1);
for kk = 1:ndur
    m = mean(traces(:, :, kk), 2);
    sem = std(traces(:, :, kk), 1, 2) / sqrt(size(traces, 2));
    m(tx > (dur_list(kk) + 300)) = nan;
    sem(tx > (dur_list(kk) + 300)) = nan;
    shadedErrorBar(tx, m, sem, 'lineProps', ...
        {'Color', cmap(kk, :)});
    hs(kk) = plot(tx, m, 'Color', cmap(kk, :), 'LineWidth', 1.5);
end
set(gca, 'TickDir', 'out', 'XLim', [tx(1), tx(end)], 'YLim', [-10, 40]);
xline(0, ':k');
yline(0, ':k');
xlabel('Time (sec)');
ylabel('dF/F (%)');
h = legend(hs, dur_label, 'Box', 'off');
h.Position = [0.4, 0.6, 0.1, 0.1];

print(gcf, 'dbs_duration_traces', '-djpeg');
print(gcf, 'dbs_duration_traces', '-dpdf', '-r300', '-vector');

%%
writetable(array2table(permute(peaks, [2, 3, 1]), "RowNames", ...
    arrayfun(@(x)(x.fname), ds, 'UniformOutput', false), "VariableNames", ...
    dur_label), 'dbs_duration_peaks.xlsx', 'WriteMode', 'replacefile');

writetable(array2table(permute(aucs, [2, 3, 1]), "RowNames", ...
    arrayfun(@(x)(x.fname), ds, 'UniformOutput', false), "VariableNames", ...
    dur_label), 'dbs_duration_aucs.xlsx', 'WriteMode', 'replacefile');

%%
stats = arrayExportStats({permute(aucs, [2, 3, 1])}, ...
    dur_label, {''}, 'AUC', ...
    'dbs_duration_aucs.xlsx', 'ranova');

figure; set(gcf, 'Position', [50, 50, 780, 780]);

subplot('Position', [0.4, 0.5, 0.20, 0.25]); hold on;
for kk = 1:ndur
    vmat = reshape(aucs(:, :, kk), [], 1);
    bar(kk, mean(vmat), 'FaceColor', 'w', 'EdgeColor', 'k', ...
        'LineWidth', 1.3, 'BarWidth', 0.67);
    errorbar(kk, mean(vmat), std(vmat)/sqrt(length(vmat)), ...
        'Color', 'k', 'LineWidth', 1.3, 'CapSize', 8);
    swarmchart(jitterX(vmat, 400) * 0.3 + kk, vmat, 15, cmap(kk, :), ...
        'filled', 'XJitter', 'none');
end
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 1], 'XTick', 1:ndur, ...
    'XTickLabel', dur_label, 'XLim', [0.3, ndur+0.7], 'FontSize', 10, ...
    'LineWidth', 1.3, 'YTickLabel', get(gca, 'YTick')/100);
ylabel(sprintf('\\Sigma_{Floure.} (\\DeltaF/F) (x10^{2} %%)'), 'FontSize', 10);
print(gcf, 'dbs_duration_aucs', '-djpeg');
print(gcf, 'dbs_duration_aucs', '-dpdf', '-r300', '-vector');


%%
raw_dir = '..\WT_ANT_GFAP_Cas3_GAdo';

dbs_list = { ...
    'WT_ANT_GFAP_Cas3_GAdo_M2_DBS1_260127', [50, 100, 150, 200, 250, 300]; ...
    'WT_ANT_GFAP_Cas3_GAdo_M4_DBS1_260127', [50, 50, 100, 150, 200, 250, 300]; ... 
    'WT_ANT_GFAP_Cas3_GAdo_M7_DBS1_260128', [50, 50, 100, 150, 200, 250, 300]; ...
    'WT_ANT_GFAP_Cas3_GAdo_M9_DBS1_260129', [50, 100, 150, 200, 250, 300]; ... #
    'WT_ANT_GFAP_Cas3_GAdo_M10_DBS1_260128', [50, 100, 150, 200, 250, 300]; ...
    'WT_ANT_GFAP_Cas3_GAdo_M4_DBS2_260130', [50, 100, 150, 200, 250, 300]; ... #
    'WT_ANT_GFAP_Cas3_GAdo_M7_DBS2_260130', [50, 100, 150, 200, 250, 300]; ...
    'WT_ANT_GFAP_Cas3_GAdo_M9_DBS2_260131', [50, 100, 150, 200, 250, 300]; ... #
    'WT_ANT_GFAP_Cas3_GAdo_M10_DBS2_260131', [50, 100, 150, 200, 250, 300]; ...
    };

%%
ds = struct;
bases = nan(size(dbs_list, 1), 1);
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
    % sig = sig / (sqrt(diff(prctile(sig, [1, 99]))) * 0.1);

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
    % print(gcf, fname, '-djpeg'); %close;
end

%%
amp_list = [50, 100, 150, 200, 250, 300];
namp = length(amp_list);
tx = (-30:0.5:500);
for nn = 1:numel(ds)
    sig = ds(nn).sig;
    ts = ds(nn).ts;    

    [traces, peaks, aucs] = deal(cell(namp, 1));
    for kk = 1:namp
        isin = ds(nn).amp==amp_list(kk);
        mrk = ds(nn).mrk(isin);

        evt = mrk + [-30, 0];
        [bsl, ~] = timestampsInterpEvents(ds(nn).sig, ds(nn).ts, evt, 'avgrg');

        evt = mrk + tx;
        [omat, ~] = timestampsInterpEvents(ds(nn).sig, ds(nn).ts, evt, 'trace');
        traces{kk} = mean(omat - bsl, ndims(omat));

        evt = mrk + tx(tx >= 0 & tx <= 10);
        [omat, ~] = timestampsInterpEvents(ds(nn).sig, ds(nn).ts, evt, 'trace');
        aucs{kk} = mean(mean(omat - bsl) * 10, ndims(omat));

        evt = mrk + [0, 60];
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
amp_label = arrayfun(@(x)(sprintf('%d \\muA', x)), amp_list, 'UniformOutput', false);
figure; set(gcf, 'Position', [50, 50, 780, 780]);

subplot('Position', [0.1, 0.5, 0.25, 0.25]); hold on;
hs = nan(kk, 1);
for kk = 1:namp
    m = mean(traces(:, :, kk), 2);
    sem = std(traces(:, :, kk), 1, 2) / sqrt(size(traces, 2));
    shadedErrorBar(tx, m, sem, 'lineProps', ...
        {'Color', cmap(kk, :)});
    hs(kk) = plot(tx, m, 'Color', cmap(kk, :), 'LineWidth', 1.5);
end
set(gca, 'TickDir', 'out', 'XLim', [tx(1), tx(end)], 'YLim', [-10, 40]);
xline(0, ':k');
yline(0, ':k');
xlabel('Time (sec)');
ylabel('\DeltaF/F (%)');
h = legend(hs, amp_label, 'Box', 'off');
h.Position = [0.4, 0.6, 0.1, 0.1];

print(gcf, 'dbs_amplitude_traces', '-djpeg');
print(gcf, 'dbs_amplitude_traces', '-dpdf', '-r300', '-vector');

%%
writetable(array2table(permute(peaks, [2, 3, 1]), "RowNames", ...
    arrayfun(@(x)(x.fname), ds, 'UniformOutput', false), "VariableNames", ...
    strrep(amp_label, '\mu', 'u')), 'dbs_amplitude_peaks.xlsx', 'WriteMode', 'replacefile');

writetable(array2table(permute(aucs, [2, 3, 1]), "RowNames", ...
    arrayfun(@(x)(x.fname), ds, 'UniformOutput', false), "VariableNames", ...
    strrep(amp_label, '\mu', 'u')), 'dbs_amplitude_aucs.xlsx', 'WriteMode', 'replacefile');

%%
stats = arrayExportStats({permute(peaks, [2, 3, 1])}, ...
    strrep(amp_label, '\mu', 'u'), {''}, 'Peak', ...
    'dbs_amplitude_peaks.xlsx', 'ranova');

figure; set(gcf, 'Position', [50, 50, 780, 780]);

subplot('Position', [0.4, 0.5, 0.25, 0.25]); hold on;
for kk = 1:namp
    vmat = reshape(peaks(:, :, kk), [], 1);
    bar(kk, mean(vmat), 'FaceColor', 'w', 'EdgeColor', 'k', ...
        'LineWidth', 1.3, 'BarWidth', 0.67);
    errorbar(kk, mean(vmat), std(vmat)/sqrt(length(vmat)), ...
        'Color', 'k', 'LineWidth', 1.3, 'CapSize', 8);
    swarmchart(jitterX(vmat, 2) * 0.3 + kk, vmat, 15, cmap(kk, :), ...
        'filled', 'XJitter', 'none');
end
set(gca, 'TickDir', 'out', 'XTick', 1:namp, 'XTickLabel', amp_label, ...
    'XLim', [0.3, namp+0.7], 'FontSize', 10, 'LineWidth', 1.3, 'YLim', [0, 60]);
ylabel(sprintf('Peak \\DeltaF/F (%%)'), 'FontSize', 10);

print(gcf, 'dbs_amplitude_peaks', '-djpeg');
print(gcf, 'dbs_amplitude_peaks', '-dpdf', '-r300', '-vector');




















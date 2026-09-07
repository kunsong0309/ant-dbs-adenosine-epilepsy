%%
prep_dir = 'prep_sig';
out_dir = 'results';

%%
load('ds_ado_dbs_ptz_seizure_20240808.mat');

for nn = 1:length(ds)
    fpath = fullfile(prep_dir, sprintf('%s_%s.mat', ds(nn).uid, ds(nn).fid));
    load(fpath, 'ff470', 'ff410', 'ts_ff', 'ts');

    if contains(ds(nn).subject, 'EGFP')
        p = polyfit(ff410, ff470, 1);
        sff = ff470 - polyval(p, ff410);
    else
        sff = ff470;
    end

    h = design(fdesign.highpass('N,F3db', 4, 0.01, 1/median(diff(ts_ff))), 'butter');
    sig = ff470;
    sig(1:100) = mean(sig(100:500));
    sig = filtfilt(h.sosMatrix, h.ScaleValues, sig);
    sigma = mad(sig(ts_ff < ds(nn).t_med), 1) * 1.4826;

    if ~isempty(ds(nn).t_dbs)
        trg = ds(nn).t_dbs(1, 1) + (-120:0.5:1800);

        sig = movmean(ds(nn).ff_dbs_on_norm, 10);
        tx = trg - ds(nn).t_dbs(1, 1);
        [pkv, pkt] = findpeaks(sig(tx > 0), tx(tx > 0));%, 'MinPeakHeight', max(sig)/2);
        if isempty(pkv); pkv = nan; pkt = nan; end
        ds(nn).ado_dbs_pk = pkv(1);
        ds(nn).ado_dbs_pkt = pkt(1);
        ds(nn).ado_dbs_slope = pkv(1) / pkt(1);

    end
end

in_dbs = arrayfun(@(x)(contains(x.subject, 'WT_ANT_M')), ds);
in_egfp = arrayfun(@(x)(contains(x.subject, 'WT_ANT_EGFP_DBS_M')), ds);
in_ptz = arrayfun(@(x)(contains(x.subject, 'WT_ANT_C_M')), ds);

grp_idx = {find(in_dbs), find(in_egfp), find(in_ptz)};
grp_name = {'PTZ+DBS', 'EGFP', 'PTZ+sham'};
grp_tag = cellfun(@(x, y)(sprintf('%s (N=%d)', x, numel(y))), ...
    grp_name, grp_idx, 'UniformOutput', false);
grp_clr = [0.96, 0.55, 0.39; 0.99, 0.85, 0.18; 0.96, 0.25, 0.39];

ds2 = ds;

%%
fig_name = 'ff_trace_egfp_dbs_on';
xlim = [-20, 100];
tx2 = (-60:0.5:1200);
tx1 = (-120:0.5:1800);

dat = cell(1, 2);
hs = [];
figure; set(gcf, 'Position', [100, 50, 480, 360]); hold on;

vmat = arrayfun(@(x)(x.ff_dbs_on'), ds2(in_egfp), 'UniformOutput', false);
vmat = cat(1, vmat{:});
dat{1} = vmat;
ic1 = grp_clr(2, :);
plot(tx1, vmat, 'Color', ic1 + (1 - ic1) * 0.7, 'LineWidth', 0.1);

vmat = arrayfun(@(x)(x.ado_dbs_on_norm(:, ...
    find(x.dbs_cur == cur_list(1), 1, 'last'))'), ds1(in_cur), 'UniformOutput', false);
vmat = cat(1, vmat{:});
dat{2} = vmat;
ic2 = cur_clr(1, :);
plot(tx2, vmat, 'Color', ic2 + (1 - ic2) * 0.7, 'LineWidth', 0.1);

hs(1) = plot(tx1, mean(dat{1}, 1, 'omitnan'), 'Color', ic1, 'LineWidth', 1.5);
hs(2) = plot(tx2, mean(dat{2}, 1, 'omitnan'), 'Color', ic2, 'LineWidth', 1.5);

yline(0, ':k');
xline([0, 10], ':k');
set(gca, 'XLim', xlim, 'TickDir', 'out', 'FontSize', 12);
legend(hs, {'EGFP', 'DBS'}, 'Box', 'off', 'Location', 'best', 'FontSize', 10);
xlabel('Time from DBS on (s)', 'FontSize', 12);
ylabel('\DeltaF/F (%)', 'FontSize', 12);
title('ANT DBS Ado effect', 'FontSize', 10);
print(gcf, fullfile(out_dir, fig_name), '-djpeg');
print(gcf, fullfile(out_dir, fig_name), '-dpdf', '-r300', '-vector');

%%
vmat = { ...
    arrayfun(@(x)(x.ado_dbs_pk(find(x.dbs_cur == cur_list(1), 1, 'last'))'), ds1(in_cur)), ...
    cat(1, ds2(in_egfp).ado_dbs_pk)};
vname = 'Peak \DeltaF/F (%)';
ttl_name = 'Ado Peak';
fig_name = 'ff_peak_dbs_on';

figure; set(gcf, 'Position', [100, 100, 360, 360]);
arrayPlotBarScatter(vmat, {''}, {'DBS', 'EGFP'}, [ic2; ic1], [], vname, ttl_name);
print(gcf, fig_name, '-djpeg');
print(gcf, fig_name, '-dpdf', '-r300', '-vector');

arrayExportStats(vmat, {''}, {'DBS', 'EGFP'}, 'AUC', fullfile(out_dir, [fig_name, '.xlsx']));

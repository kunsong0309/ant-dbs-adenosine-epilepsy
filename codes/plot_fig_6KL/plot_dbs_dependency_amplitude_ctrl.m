%%
prep_dir = '..\prep_sig';
out_dir = 'results';

load('ds_info.mat', 'ds');

inc_ls = { ...
    'WT_ANT_ADO_RASTRO_M1_DBS1_1118'; ...
    'WT_ANT_ADO_RASTRO_M4_DBS1_1118'; ...
    'WT_ANT_ADO_RASTRO_M5_DBS1_1213'; ...
    'WT_ANT_ADO_RCAMP_M1_DBS1_0903'; ...
    'WT_ANT_ADO_RCAMP_M5_DBS1_1101'; ...
    'WT_ANT_HPMCA_ADO_M1_DBS2_1106'; ...
    'WT_ANT_HPMCA_ADO_M8_DBS1_0119'; ...
    'WT_ANT_HPMCA_ADO_M9_DBS1_0114'; ...
    };

inc_ls = unique(inc_ls, 'stable');

ds = ds(arrayfun(@(x)(ismember(x.fid, inc_ls)), ds));
%%
cur_list = [50, 100, 150, 200, 250, 300];
cur_name = arrayfun(@(x)(strcat(num2str(x), ' uA')), cur_list, 'UniformOutput', false);
cmap = colormap(rainbow);
cur_clr = cmap([1, 45, 80, 150, 200, 256], :);

rec_dbs_current = { ...
    'WT_ANT_ADO_RASTRO_M1_DBS1_1118', [50, 100, 150, 200, 250, 300]; ...
    'WT_ANT_ADO_RASTRO_M4_DBS1_1118', [50, 100, 150, 200, 250, 300]; ...
    'WT_ANT_ADO_RASTRO_M5_DBS1_1213', [50, 100, 150, 200, 250, 300]; ...
    'WT_ANT_ADO_RCAMP_M1_DBS1_0903', [50, 100, 150, 200, 250, 300]; ...
    'WT_ANT_ADO_RCAMP_M5_DBS1_1101', [50, 100, 150, 200, 250, 300]; ...
    'WT_ANT_HPMCA_ADO_M1_DBS2_1106', [50, 100, 150, 200, 250, 300]; ...
    'WT_ANT_HPMCA_ADO_M8_DBS1_0119', [50, 100, 150, 200, 250, 300]; ...
    'WT_ANT_HPMCA_ADO_M9_DBS1_0114', [50, 100, 150, 200, 250, 300]; ...
    };

for nn = 1:length(ds)
    idx = find(strcmp(ds(nn).fid, rec_dbs_current(:, 1)));
    if ~any(idx)
        ds(nn).dbs_cur = nan;
        continue; 
    end
    cur = rec_dbs_current{idx, 2}(:);
    if numel(rec_dbs_current{idx, 2}) == 1
        cur = repmat(cur, size(ds(nn).t_dbs, 1), 1);
    end
    ds(nn).dbs_cur = cur;

    if ~isempty(ds(nn).t_dbs)
        dur = ds(nn).t_dbs(:, 2) - ds(nn).t_dbs(:, 1);
        dur = round(dur / 10) * 10;
    else
        dur = nan;
    end
    ds(nn).dbs_dur = dur;
end

%%
for nn = 1:length(ds)
    fpath = fullfile(prep_dir, sprintf('%s_%s.mat', ds(nn).uid, ds(nn).fid));
    load(fpath, 'ff470', 'ts_ff');
    sado = sigGetValid(ff470(:, end));
    if isempty(sado); continue; end 
    if isempty(ds(nn).t_dbs); continue; end

    trg = ds(nn).t_dbs(:, 1) + [-30, 0];
    ds(nn).ado_dbs_on_bsl = timestampsInterpEvents(sado, ts_ff, trg, ...
        'avgrg', 'linear', 2);

    trg = ds(nn).t_dbs(:, 1) + (-60:0.5:1200);
    ds(nn).ado_dbs_on_norm = timestampsInterpEvents(sado, ts_ff, trg, ...
        'trace', 'linear', 2) - ds(nn).ado_dbs_on_bsl;

    % trg = ds(nn).t_dbs;
    % ds(nn).auado_dbs_on_off_norm = ( ...
    %     timestampsInterpEvents(sado, ts_ff, trg, 'avgrg', 'linear', 2) ...
    %     - ds(nn).ado_dbs_on_bsl) .* diff(trg, 1, 2)';

    trg = ds(nn).t_dbs(:, 1) + [0, 300];
    ds(nn).auado_dbs_on_off_norm = ( ...
        timestampsInterpEvents(sado, ts_ff, trg, 'avgrg', 'linear', 2) ...
        - ds(nn).ado_dbs_on_bsl) .* diff(trg, 1, 2)';

    trg = ds(nn).t_dbs(:, 1) + (0:0.5:60);
    ds(nn).ado_dbs_pk = max(movmean(timestampsInterpEvents(sado, ts_ff, trg, ...
        'trace', 'linear', 2) - ds(nn).ado_dbs_on_bsl, 10, 1, 'SamplePoints', trg(1, :)));
end

%%
ttl_name = 'Ado Change with DBS Current';
fig_name = 'dbs_dependency_current_ado_trace';
xlim = [-30, 500];
tx = (-60:0.5:1200);
ncond = length(cur_list);

figure; set(gcf, 'Position', [100, 100, 600, 400]); hold on; 
for cc = 1:ncond
    vmat = arrayfun(@(x)(x.ado_dbs_on_norm(:, ...
        find(x.dbs_cur == cur_list(cc), 1, 'last'))'), ds(in_cur), 'UniformOutput', false);
    vmat = cat(1, vmat{:});
    vmean = mean(vmat, 1, 'omitnan');
    vsem = std(vmat, 0, 1, 'omitnan') / sqrt(size(vmat, 1));

    shadedErrorBar(tx, vmean, vsem, ...
        'lineProps', {'Color', cur_clr(cc, :), 'LineWidth', 1.5});
    % plot(tx, vmat(:, :, cc), 'Color', cur_clr(cc, :));
end
yline(0, ':k');
xline([0, 10], ':k');
set(gca, 'XLim', xlim, 'TickDir', 'out', 'FontSize', 12);
legend(cur_name, 'Box', 'off', 'Location', 'eastoutside', 'FontSize', 10);
xlabel('Time from DBS on (s)', 'FontSize', 12);
ylabel('\DeltaF/F (%)', 'FontSize', 12);
title(sprintf('%s', ttl_name));
print(gcf, fullfile(out_dir, fig_name), '-djpeg');
print(gcf, fullfile(out_dir, fig_name), '-dpdf', '-r300', '-vector');

%%
ttl_name = 'Ado Total Change with DBS Current';
fig_name = 'dbs_dependency_current_ado_auc';
ncond = length(cur_list);
vmat = cell(1, ncond);

figure; set(gcf, 'Position', [100, 100, 400, 400]); hold on; 
for cc = 1:ncond
    imat = arrayfun(@(x)(x.auado_dbs_on_off_norm( ...
        find(x.dbs_cur == cur_list(cc), 1, 'last'))'), ds(in_cur), 'UniformOutput', false);
    vmat{cc} = cat(1, imat{:});
    vmean = mean(vmat{cc}, 1, 'omitnan');
    vsem = std(vmat{cc}, 0, 1, 'omitnan') / sqrt(size(vmat{cc}, 1));

    bar(cc, vmean, 'FaceColor', cur_clr(cc, :), 'EdgeColor', 'w');
    % errorbar(cc, vmean, vsem, 'ok', 'MarkerSize', 3, ...
    %     'MarkerFaceColor', 'k', 'LineWidth', 1.3);
end
plot(1:ncond, cat(2, vmat{:}), 'o-', 'Color', [.3, .3, .3], 'MarkerSize', 3, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
set(gca, 'TickDir', 'out', 'TickLength', [0.005, 0.1], 'FontSize', 12, ...
    'XTick', 1:ncond, 'XTickLabel', cur_name, 'XTickLabelRotation', 45, ...
    'XLim', [0.3, ncond+0.7]);
ylabel('Area under \DeltaF/F (%)', 'FontSize', 12);
title(sprintf('%s', ttl_name));
print(gcf, fullfile(out_dir, fig_name), '-djpeg');
print(gcf, fullfile(out_dir, fig_name), '-dpdf', '-r300', '-vector');

cnames = arrayfun(@(x)(sprintf('y%d', x)), 1:length(cur_name), 'UniformOutput', false);
between = array2table(cat(2, vmat{:}), 'VariableNames', cnames);
within = table(cnames(:), 'VariableNames', {'Current'});
rm = fitrm(between, [strjoin(cnames, ','), ' ~ 1'], 'WithinDesign', within);
out_name = fullfile(out_dir, [fig_name, '.xlsx']);

between = array2table(cat(2, vmat{:}), 'VariableNames', cur_name);
writetable(between, out_name, 'Sheet', 'Data', 'WriteMode', 'replacefile');
writetable(ranova(rm), out_name, 'WriteRowNames', true, ...
    'Sheet', 'Stats', 'Range', 'A1');
tbl = multcompare(rm, 'Current');
tbl.Current_1 = replace(tbl.Current_1, cnames, cur_name);
tbl.Current_2 = replace(tbl.Current_2, cnames, cur_name);
writetable(tbl, out_name, 'WriteRowNames', true, 'Sheet', 'Stats', 'Range', 'A10');

%%
ttl_name = 'Ado Peak Change with DBS Current';
fig_name = 'dbs_dependency_current_ado_peak';
ncond = length(cur_list);
vmat = cell(1, ncond);

figure; set(gcf, 'Position', [100, 100, 400, 400]); hold on; 
for cc = 1:ncond
    imat = arrayfun(@(x)(x.ado_dbs_pk( ...
        find(x.dbs_cur == cur_list(cc), 1, 'last'))'), ds(in_cur), 'UniformOutput', false);
    vmat{cc} = cat(1, imat{:});
    vmean = mean(vmat{cc}, 1, 'omitnan');
    vsem = std(vmat{cc}, 0, 1, 'omitnan') / sqrt(size(vmat{cc}, 1));
end
arrayPlotBarScatterPaired(vmat, {''}, cur_name, cur_clr, [], 'Peak \DeltaF/F (%)', ttl_name);
print(gcf, fullfile(out_dir, fig_name), '-djpeg');
print(gcf, fullfile(out_dir, fig_name), '-dpdf', '-r300', '-vector');

vmat = {cat(2, vmat{:})};
arrayExportStats(vmat, cur_name, {''}, 'Current', fullfile(out_dir, [fig_name, '.xlsx']));





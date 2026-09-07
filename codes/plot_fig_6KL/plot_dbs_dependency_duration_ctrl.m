%%
prep_dir = '..\prep_sig';
out_dir = 'results';

load('ds_info.mat', 'ds');

dur_list = [10, 60, 300, 600, 900];
dur_name = arrayfun(@(x)(strcat(num2str(x), ' s')), dur_list, 'UniformOutput', false);
cmap = rainbow();
dur_clr = cmap([1, 45, 80, 150, 200, 256], :);
%%
for nn = 1:length(ds)
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

    trg = ds(nn).t_dbs(:, 1) + (-60:0.5:1600);
    ds(nn).ado_dbs_on_norm = timestampsInterpEvents(sado, ts_ff, trg, ...
        'trace', 'linear', 2) - ds(nn).ado_dbs_on_bsl;

    trg = ds(nn).t_dbs;
    ds(nn).auado_dbs_on_off_norm = ( ...
        timestampsInterpEvents(sado, ts_ff, trg, 'avgrg', 'linear', 2) ...
        - ds(nn).ado_dbs_on_bsl) .* diff(trg, 1, 2)';
end

%%
ttl_name = 'Ado Change with DBS Duration';
fig_name = 'ado_trace_dbs_duration';
xlim = [-30, 1600];
tx = (-60:0.5:1600);
ncond = length(dur_list);

figure; set(gcf, 'Position', [100, 100, 600, 400]); hold on; 
for cc = 1:ncond
    if dur_list(cc) > 100
        vmat = arrayfun(@(x)(x.ado_dbs_on_norm(:, ...
            find(x.dbs_dur(1) == dur_list(cc)))'), ds(in_dur), 'UniformOutput', false);
    else
        vmat = arrayfun(@(x)(x.ado_dbs_on_norm(:, ...
            find(x.dbs_dur == dur_list(cc), 1))'), ds(in_dur), 'UniformOutput', false);
    end
    vmat = cat(1, vmat{:});
    vmean = mean(vmat, 1, 'omitnan');
    vsem = std(vmat, 0, 1, 'omitnan') / sqrt(size(vmat, 1));

    tin = double(tx < dur_list(cc) + 500);
    tin(tin == 0) = nan;
    shadedErrorBar(tx, vmean.*tin, vsem.*tin, ...
        'lineProps', {'Color', dur_clr(cc, :), 'LineWidth', 1.5});
    % plot(tx, vmat(:, :, cc).*tin, 'Color', cond_clr(cc, :));
end
yline(0, ':k');
set(gca, 'XLim', xlim, 'TickDir', 'out', 'FontSize', 12);
legend(dur_name, 'Box', 'off', 'Location', 'eastoutside', 'FontSize', 10);
xlabel('Time from DBS on (s)', 'FontSize', 12);
ylabel('\DeltaF/F (%)', 'FontSize', 12);
title(sprintf('%s', ttl_name));
% print(gcf, fullfile(out_dir, fig_name), '-djpeg');
% print(gcf, fullfile(out_dir, fig_name), '-dpdf', '-r300', '-vector');

%%
ttl_name = 'Ado Total Change with DBS Duration';
fig_name = 'ado_auc_dbs_duration';
ncond = length(dur_list);
vmat = cell(1, ncond);

figure; set(gcf, 'Position', [100, 100, 400, 400]); hold on; 
for cc = 1:ncond
    if dur_list(cc) > 100
        imat = arrayfun(@(x)(x.auado_dbs_on_off_norm(:, ...
            find(x.dbs_dur(1) == dur_list(cc)))'), ds(in_dur), 'UniformOutput', false);
    else
        imat = arrayfun(@(x)(x.auado_dbs_on_off_norm(:, ...
            find(x.dbs_dur == dur_list(cc), 1))'), ds(in_dur), 'UniformOutput', false);
    end
    subm(:, cc-1) = cellfun(@(x)(~isempty(x)), imat);
    vmat{cc} = cat(1, imat{:});
    vmean = mean(vmat{cc}, 1, 'omitnan');
    vsem = std(vmat{cc}, 0, 1, 'omitnan') / sqrt(size(vmat{cc}, 1));

    bar(cc, vmean, 'FaceColor', dur_clr(cc, :), 'EdgeColor', 'w');
    errorbar(cc, vmean, vsem, 'ok', 'MarkerSize', 3, ...
        'MarkerFaceColor', 'k', 'LineWidth', 1.3);
end
% plot(1:ncond, vmat, 'o-', 'Color', [.3, .3, .3], 'MarkerSize', 3, ...
%         'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
set(gca, 'TickDir', 'out', 'TickLength', [0.005, 0.1], 'FontSize', 12, ...
    'XTick', 1:ncond, 'XTickLabel', dur_name, 'XTickLabelRotation', 45, ...
    'XLim', [0.3, ncond+0.7]);
ylabel('Area under \DeltaF/F (%)', 'FontSize', 12);
title(sprintf('%s', ttl_name));
print(gcf, fullfile(out_dir, fig_name), '-djpeg');
print(gcf, fullfile(out_dir, fig_name), '-dpdf', '-r300', '-vector');

grp = cellfun(@(x, y)(repmat(string(x), size(y, 1), 1)), dur_name, vmat, 'UniformOutput', false);
between = table(cat(1, grp{:}), cat(1, vmat{:}), 'VariableNames', {'Group', 'Value'});
rm = anova(between, 'Value ~ Group');
out_name = fullfile(out_dir, [fig_name, '.xlsx']);

writetable(between, out_name, 'Sheet', 'Data', 'WriteMode', 'replacefile');
writetable(stats(rm), out_name, 'WriteRowNames', true, ...
    'Sheet', 'Stats', 'Range', 'A1');
tbl = multcompare(rm, 'Group');
writetable(tbl, out_name, 'WriteRowNames', true, 'Sheet', 'Stats', 'Range', 'A10');





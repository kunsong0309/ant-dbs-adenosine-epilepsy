%%
out_dir = 'results';
load('20240914_ds_stats_bs_se.mat', 'stats');
fds = {'group', 'batch', 'mouse', ...
    'spk_zt_day', 'ea_zt_day', 'spk_frequency', 'spk_ph_bs', ...
    'ea_ph_bs', 'ea_frequency', 'ea_duration', 'ea_per_hour', ...
    'se_frequency', 'se_duration', 'se_per_day', 'seizure_ts'};
stats = rmfield(stats, setdiff(fields(stats), fds));
stats(1).seizure_ts = [];

load('20250119_ds_seizure_3_100.mat', 'ds_ant', 'ds_ca1');
ds = ds_ant;
for nn = 1:length(ds)
    ds(nn).group = 'KO';
end
stats = cat(1, stats, rmfield(ds, setdiff(fields(ds), fds)));

% load('ds_seizure_ca1_20250106.mat', 'ds');
ds = ds_ca1;
for nn = 1:length(ds)
    ds(nn).group = 'CA1_KO';
end
stats = cat(1, stats, rmfield(ds, setdiff(fields(ds), fds)));

ds = load('20250119_ds_ka1_seizure_stats_highfreq.mat', 'stats');
ds = ds.stats;
stats = cat(1, stats, rmfield(ds, setdiff(fields(ds), fds)));

ds = load('20250119_ds_ka_info_stats_highfreq.mat', 'ds');
ds = ds.ds;
for nn = 1:length(ds)
    ds(nn).group = 'KA';
end
stats = cat(1, stats, rmfield(ds, setdiff(fields(ds), fds)));

load('20250124_ds_oe_info_stats_highfreq.mat', 'ds_oe_1', 'ds_oe_2');
ds = ds_oe_1;
for nn = 1:length(ds)
    ds(nn).group = 'OE';
end
stats = cat(1, stats, rmfield(ds, setdiff(fields(ds), fds)));

ds = ds_oe_2;
for nn = 1:length(ds)
    ds(nn).group = 'OE';
end
stats = cat(1, stats, rmfield(ds, setdiff(fields(ds), fds)));

%%
fig_name = 'fig5h_ka_oe';
var_info = cell2table({ ... 
    'se_frequency', 'Seizure (per day)', 'Seizure Frequency'; ...
    'se_duration', 'Seizure duration (s)', 'Seizure Duration'; ...
    'se_per_day', 'Seizure time (s/day)', 'Seizure Time'; ...
    'ea_frequency', 'EA (per hour)', 'EA Frequency'; ...
    'ea_duration', 'EA (s)', 'EA Duration'; ...
    'ea_per_hour', 'EA (s/hr)', 'EA Time'; ...
    'spk_frequency', 'IED (per hour)', 'IED Frequency'; ...
    }, 'VariableNames', {'field', 'ylabel', 'title'});

grp_ls = [5, 6];
nvar = size(var_info, 1);
for vv = 1:nvar
    vmat = arrayfun(@(x)(x.(var_info.field{vv})), stats, 'UniformOutput', false);
    vname = var_info.ylabel{vv};
    ttl_name = var_info.title{vv};
    var_name = var_info.field{vv};

    figure; set(gcf, 'Position', [100, 100, 480, 360]); hold on;
    vmat = cellfun(@(x)(cat(1, vmat{x})), grp_idx(grp_ls), 'UniformOutput', false);
    arrayPlotBarScatter(vmat, {''}, grp_name(grp_ls), ...
        grp_clr(grp_ls, :), [], vname, ttl_name);
    set(gca, 'TickDir', 'out');
    ylabel(vname, 'FontSize', 12);
    title(sprintf('%s', ttl_name), 'FontSize', 12);
    print(gcf, fullfile(out_dir, [fig_name, '_', var_name]), '-djpeg');
    print(gcf, fullfile(out_dir, [fig_name, '_', var_name]), '-dpdf', '-r300', '-vector');

    arrayExportStats(vmat, {'value'}, grp_name(grp_ls), vname, ...
        fullfile(out_dir, [fig_name, '_', var_name, '.xlsx']));
end

%%
in_pump = arrayfun(@(x)(strcmp(x.group, 'pump')), stats);
in_bsl = arrayfun(@(x)(strcmp(x.group, 'baseline')), stats);
in_ko = arrayfun(@(x)(strcmp(x.group, 'KO')), stats);
in_ca1 = arrayfun(@(x)(strcmp(x.group, 'CA1_KO')), stats);
in_ka = arrayfun(@(x)(strcmp(x.group, 'KA')), stats);
in_oe = arrayfun(@(x)(strcmp(x.group, 'OE')), stats);

grp_idx = {find(in_bsl), find(in_pump), find(in_ko), find(in_ca1), find(in_ka), find(in_oe)};
grp_name = {'Control', 'Ado', 'ANT_KO', 'CA1_KO', 'KA', 'OE'};
grp_clr = [0.70, 0.70, 0.70; 0.30, 0.40, 0.98; 0.99, 0.40, 0.50; 0.80, 0.90, 0.40; 0.30, 0.40, 0.98; 0.70, 0.70, 0.70];




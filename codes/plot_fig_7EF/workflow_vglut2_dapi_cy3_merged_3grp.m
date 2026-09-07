%%
tot_mm3 = 318.2 * 318.2 * 6 / 1e9; % mm3
tot_mm2 = 318.2 * 318.2 / 1e6; % mm2
ds = dir('mages\Confocal\40x_PTZ_Only_6um\*.tif');
ds = rmfield(ds, {'date', 'bytes', 'isdir', 'datenum'});
for nn = 1:length(ds)
    ds(nn).fpath = fullfile(ds(nn).folder, ds(nn).name);
    ds(nn).fid = regexp(ds(nn).name, '\S+(?=[.])', 'match', 'once');

    ds(nn).M = str2double(regexp(ds(nn).fid, '(?<=_M)\d+(?=_)', 'match', 'once'));
    ds(nn).P = str2double(regexp(ds(nn).fid, '(?<=_P)\d+(?=_)', 'match', 'once'));
    if isnan(ds(nn).P)
        ds(nn).P = str2double(regexp(ds(nn).fid, '(?<=_)00\d', 'match', 'once'));
    end

    ds(nn).group = 'ptz';

    nch = length(imfinfo(ds(nn).fpath));
    im = cell(nch, 1);
    for kk = 1:nch
        im{kk} = imread(ds(nn).fpath, kk);
    end
    ds(nn).im = cat(3, im{:});
end

%%
for nn = 1:length(ds)
    im = double(ds(nn).im) / 1000;

    ds(nn).dapi_ind = imMaskDAPI(im(:, :, 1));
    % figure; set(gcf, 'Position', [100, 100, 600, 600], 'Visible', 'off'); 
    % imshow(labeloverlay(im(:, :, 1), ds(nn).dapi_ind > 0, 'Transparency', 0.5));
    % title(ds(nn).fid, 'Interpreter', 'none');
    % print(gcf, sprintf('%.3d_%s_dapi', nn, ds(nn).fid), '-djpeg');
    % close;

    [ds(nn).vglut2_ind, ds(nn).vglut2_area, ds(nn).vglut2_pow] = imMaskVglut2(im(:, :, 2));
    % figure; set(gcf, 'Position', [100, 100, 600, 600], 'Visible', 'off'); 
    % imshow(labeloverlay(im(:, :, 2), ds(nn).vglut2_ind > 0, 'Transparency', 0.5));
    % title(ds(nn).fid, 'Interpreter', 'none');
    % print(gcf, sprintf('%.3d_%s_vglut2', nn, ds(nn).fid), '-djpeg');
    % close;

    ds(nn).dapi_nroi = max(ds(nn).dapi_ind, [], 'all');

    thr = median(im(:, :, 3), 'all') + mad(im(:, :, 3), 1, 'all') * 1.4826 * 5;
    [ds(nn).cy3_ndot, ds(nn).cy3_dots_ind] = imSegmentDots(im(:, :, 3), thr);

    ds(nn).cy3_ratio = mean(ds(nn).cy3_dots_ind > 0, 'all') * 100;
    ds(nn).cy3_dots_dens = ds(nn).cy3_ndot / tot_mm2;
    ds(nn).cy3_ndot_per_cell = ds(nn).cy3_ndot / ds(nn).dapi_nroi;

    ds(nn).vglut2_dapi_nroi = 0;
    msk0 = ds(nn).vglut2_ind > 0;
    for kk = 1:max(ds(nn).dapi_ind, [], 'all')
        mskk = ds(nn).dapi_ind == kk;
        if mean(msk0(mskk), 'all') == 1
            ds(nn).vglut2_dapi_nroi = ds(nn).vglut2_dapi_nroi + 1;
        end
    end

    msk = msk0 & ds(nn).cy3_dots_ind > 0;
    ds(nn).vglut2_cy3_ndot = length(unique(ds(nn).cy3_dots_ind(msk)));
    ds(nn).vglut2_cy3_ratio = mean(msk(msk0), 'all') * 100;
    ds(nn).vglut2_cy3_dots_dens = ds(nn).vglut2_cy3_ndot / (tot_mm2 * mean(msk0, 'all'));
    ds(nn).vglut2_cy3_ndot_per_cell = ds(nn).vglut2_cy3_ndot / ds(nn).vglut2_dapi_nroi;

    % figure; set(gcf, 'Position', [100, 100, 600, 600], 'Visible', 'off');
    % imshow(cat(3, ds(nn).cy3_dots_ind > 0, ds(nn).vglut2_ind > 0, ds(nn).dapi_ind > 0)*1);
    % title(ds(nn).fid, 'Interpreter', 'none');
    % print(gcf, sprintf('%.3d_%s_merged', nn, ds(nn).fid), '-djpeg');
    % close;
end

%%
ds0 = load('20240626_ds_dbs_info_res.mat', 'ds');
ds = cat(1, ds, ds0.ds);

%%
out_dir = 'results';

in_ptz = arrayfun(@(x)(strcmp(x.group, 'ptz')), ds);
in_dbs = arrayfun(@(x)(strcmp(x.group, 'dbs')), ds);
in_ctrl = arrayfun(@(x)(strcmp(x.group, 'ctrl')), ds);

grp_idx = {find(in_ctrl), find(in_dbs), find(in_ptz)};
grp_name = {'Ctrl', 'DBS', 'PTZ'};
grp_tag = cellfun(@(x, y)(sprintf('%s (N=%d)', x, numel(y))), ...
    grp_name, grp_idx, 'UniformOutput', false);
grp_clr = [0.56, 0.40, 0.67; 0.29, 0.68, 0.29; 0.55, 0.34, 0.29];

%%
ttl_name = 'CY3 + Vglut2';
var_name = 'Dots per cell (#)';
fd_name = 'vglut2_cy3_ndot_per_cell';
fig_name = 'dbs_ctrl_vglut2_cy3_ndot_per_cell';

grp = [1, 2];
ngrp = length(grp);
vmat = cell(1, ngrp);
figure; set(gcf, 'Position', [100, 100, 400, 400]); hold on;
for gg = 1:ngrp 
    % pid = find(arrayfun(@(x)(strcmp(x.group, 'ctrl') & x.M == ds(nn).M ...
    %     & x.P == ds(nn).P), ds));
    vmat{gg} = arrayfun(@(x)(x.(fd_name)), ds(grp_idx{grp(gg)}));
end
arrayPlotBarScatterPaired(vmat, {''}, grp_name(grp), grp_clr(grp, :), [], var_name, ttl_name);
print(gcf, fullfile(out_dir, fig_name), '-djpeg');
print(gcf, fullfile(out_dir, fig_name), '-dpdf', '-r300', '-vector');

arrayExportStats({cat(2, vmat{:})}, grp_name(grp), {''}, 'Condition', ...
    fullfile(out_dir, [fig_name, '.xlsx']));

%%
ttl_name = 'CY3 + Vglut2';
var_name = 'Dots per cell (#)';
fd_name = 'vglut2_cy3_ndot_per_cell';
fig_name = 'ptz_dbs_vglut2_cy3_ndot_per_cell';

grp = [3, 2];
ngrp = length(grp);
vmat = cell(1, ngrp);
figure; set(gcf, 'Position', [100, 100, 400, 400]); hold on;
for gg = 1:ngrp 
    % pid = find(arrayfun(@(x)(strcmp(x.group, 'ctrl') & x.M == ds(nn).M ...
    %     & x.P == ds(nn).P), ds));
    vmat{gg} = arrayfun(@(x)(x.(fd_name)), ds(grp_idx{grp(gg)}));
end
arrayPlotBarScatter(vmat, {''}, grp_name(grp), grp_clr(grp, :), [], var_name, ttl_name);
print(gcf, fullfile(out_dir, fig_name), '-djpeg');
print(gcf, fullfile(out_dir, fig_name), '-dpdf', '-r300', '-vector');

arrayExportStats(vmat, {''}, grp_name(grp), 'Value', ...
    fullfile(out_dir, [fig_name, '.xlsx']));


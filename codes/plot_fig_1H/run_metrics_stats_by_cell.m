%%
chunk_a = 1024;
vx_log2 = (0:0.05:8)';
vx = 2 .^ vx_log2 - 1;
precx = (0:5:100)';

nv = length(vx);
ns = length(ds);

im_fun = @(f, h1, h2, w1, w2)(cat(3, ...
    imread(f, 1, 'PixelRegion', {[h1, h2], [w1, w2]}), ...
    imread(f, 3, 'PixelRegion', {[h1, h2], [w1, w2]}), ...
    imread(f, 4, 'PixelRegion', {[h1, h2], [w1, w2]})));

%% set threshold according slice distribution
for nn = 1:ns
    isin = unique(ds(nn).tile_idx);

    % D1: 1 (g1())
    % dx = mean(ds(nn).chan_dist.gfap_th(isin, :), 1);
    % X = arrayfun(@(x)(ones(round(dx(x)), 1) * vx_log2(x)), 
    % 2:nv, 'UniformOutput', false);

    % D1: 50 (g3(m+2sd)); D2: 25 (g3(m-2sd))
    dx = mean(ds(nn).chan_dist.gfap_filt(isin, :), 1);
    X = arrayfun(@(x)(ones(round(dx(x)), 1) * vx_log2(x)), ...
        22:nv, 'UniformOutput', false);
    X = cat(1, X{:});
    ds(nn).gfap_thr = 2 ^ (median(X) + mad(X, 1) * 1.4826 * 1) - 1;
    
    % D1: 0, D2: m+1sd
    dx = mean(ds(nn).chan_dist.nt5c2_filt(isin, :), 1);
    X = arrayfun(@(x)(ones(round(dx(x)), 1) * vx_log2(x)), ...
        2:nv, 'UniformOutput', false);
    X = cat(1, X{:});
    ds(nn).nt5c2_log2_m = median(X);
    ds(nn).nt5c2_log2_sd = mad(X, 1) * 1.4826;
    ds(nn).nt5c2_thr = 2 ^ (median(X) + mad(X, 1) * 1.4826 * 1) - 1;

    % D1: 0, D2: m+1sd
    dx = mean(ds(nn).chan_dist.adora1_filt(isin, :), 1);
    X = arrayfun(@(x)(ones(round(dx(x)), 1) * vx_log2(x)), ...
        2:nv, 'UniformOutput', false);
    X = cat(1, X{:});
    ds(nn).adora1_log2_m = median(X);
    ds(nn).adora1_log2_sd = mad(X, 1) * 1.4826;
    ds(nn).adora1_thr = 2 ^ (median(X) + mad(X, 1) * 1.4826 * 1) - 1;

    % D1: 1, D2: m+5sd, D3: gm2 (dmin), D4: 10
    dx = mean(ds(nn).chan_dist.adk_filt(isin, :), 1);
    X = arrayfun(@(x)(ones(round(dx(x)), 1) * vx_log2(x)), ...
        2:nv, 'UniformOutput', false);
    X = cat(1, X{:});
    rng(6535)
    gm = fitgmdist(X, 2);
    [mu, oo] = sort(gm.mu); %mu'
    sigma = gm.Sigma(oo);
    gx = (mu(1):0.05:mu(2))';
    gy = pdf(gm, gx);
    [~, idx] = min(gy);
    % ds(nn).adk_thr = 2 ^ (gx(idx)) - 1;
    ds(nn).adk_thr = 10;
    % ds(nn).adk_thr = 2 ^ (median(X) + mad(X, 1) * 1.4826 * 3) - 1;
    X = X(X > ds(nn).adk_thr);
    ds(nn).adk_log2_m = median(X);
    ds(nn).adk_log2_sd = mad(X, 1) * 1.4826;

    % D1: 2, D2: m+1sd, D3: 2gm, m+1sd
    dx = mean(ds(nn).chan_dist.ppat_filt(isin, :), 1);
    X = arrayfun(@(x)(ones(round(dx(x)), 1) * vx_log2(x)), ...
        2:nv, 'UniformOutput', false);
    X = cat(1, X{:});
    ds(nn).ppat_thr = 2 ^ (median(X) + mad(X, 1) * 1.4826 * 1) - 1;    
    ds(nn).ppat_log2_m = median(X);
    ds(nn).ppat_log2_sd = mad(X, 1) * 1.4826;

    % D1: 0, D2: m+2sd, D3: m+1sd
    dx = mean(ds(nn).chan_dist.cd73_filt(isin, :), 1);
    X = arrayfun(@(x)(ones(round(dx(x)), 1) * vx_log2(x)), ...
        2:nv, 'UniformOutput', false);
    X = cat(1, X{:});
    ds(nn).cd73_log2_m = median(X);
    ds(nn).cd73_log2_sd = mad(X, 1) * 1.4826;
    ds(nn).cd73_thr = 2 ^ (median(X) + mad(X, 1) * 1.4826 * 1) - 1;

    % D1: m+1sd
    dx = mean(ds(nn).chan_dist.gfap_th(isin, :), 1);
    X = arrayfun(@(x)(ones(round(dx(x)), 1) * vx_log2(x)), ...
        2:nv, 'UniformOutput', false);
    X = cat(1, X{:});
    ds(nn).gfap_th_thr = 2 ^ (median(X) + mad(X, 1) * 1.4826 * 1) - 1;

    % rng(6535)
    % gm = fitgmdist(X, 2);
    % [mu, oo] = sort(gm.mu); %mu'
    % sigma = gm.Sigma(oo);
    % gx = (mu(1):0.05:mu(2))';
    % gy = pdf(gm, gx);
    % [~, idx] = min(gy);
    % % mu(2) = median(X);
    % % sigma(2) = mad(X, 1) * 1.4826;
    % % fprintf('%d %.2f (%.2f), thr: %.2f (%.2f).\n', nn, ...
    % %     mu(1) + sigma(1) * 1, 2^(mu(1) + sigma(1) * 1)-1, ...
    % %     mu(2) + sigma(2) * 1, 2^(mu(2) + sigma(2) * 1)-1);
    % fprintf('%d %.2f (%.2f), thr: %.2f (%.2f).\n', nn, ...
    %     mu(1), sigma(1), mu(2), sigma(2));
    % figure; histogram(X, vx_log2);
    % xline(mu(2) + sigma(2) * 1, ':r');
    % xline(median(X) + mad(X, 1) * 1.4826 * 1, ':m');
end

%% set 95th-perc threshold
adora1_thr = [4.58, 4.39, 4.16, 4.00, 4.40, 3.04, 4.57, 3.85, 4.06];
cd73_thr = [5.35, 5.20, 4.55, 5.35, 5.40, 5.10, 4.15, 4.85, 5.40]; 
ppat_thr = [4.95, 4.88, 4.6, 5.4, 4.9, 5.05, 4.13, 4.4, 5.7];
nt5c2_thr = [3.8, 4.40, 3.3, 6.2, 4.80, 4.4, 5.35, 4.25, 3.65];
adk_thr = [0.53, 5.8, 0.43, 5.12, 6, 1.14, 6.34, 5.55, 0.65];

for nn = 1:ns
    ds(nn).adora1_thr = 2^adora1_thr(nn)-1;
    ds(nn).cd73_thr = 2^cd73_thr(nn)-1;
    ds(nn).ppat_thr = 2^ppat_thr(nn)-1;
    ds(nn).nt5c2_thr = 2^nt5c2_thr(nn)-1;
    ds(nn).adk_thr = 2^adk_thr(nn)-1;
end

%% metrics
precx_hi = (0:0.5:100)';
precx_md = precx == 50;
for nn = 1:ns
    t0 = tic;
    nc = length(ds(nn).cc_dapi.PixelIdxList);

    vmat = zeros(nc, 19);
    for ii = 1:nc
        vmat(ii, 1) = ds(nn).circ_area{ii}(1);
        X = ds(nn).misc_metrics{ii};
        X_ = interp1(precx, X, precx_hi, 'spline');
        % neun
        vmat(ii, 2) = mean(X_(:, 2, 1) > 0.5);
        vmat(ii, 3) = X(precx_md, 2, 1);
        % gfap filt hole
        vmat(ii, 4) = mean(X_(:, 3, 1) > ds(nn).gfap_thr);
        vmat(ii, 5) = X(precx_md, 3, 1);
        % gfap filt ring
        vmat(ii, 6) = mean(X_(:, 3, 2) > ds(nn).gfap_thr);
        vmat(ii, 7) = X(precx_md, 3, 2);
        % gfap tophat ring
        vmat(ii, 8) = mean(X_(:, 4, 2) > ds(nn).gfap_th_thr);
        vmat(ii, 9) = trapz(precx, log(X(:, 4, 2)+1)/log(2)) / 100;
        % adora1 hole
        vmat(ii, 10) = mean(X_(:, 5, 1) > ds(nn).adora1_thr);
        vmat(ii, 11) = trapz(precx, log(X(:, 5, 1)+1)/log(2)) / 100;
        % vmat(ii, 11) = log(X(precx_md, 5, 1)+1)/log(2)/ds(nn).adora1_log2_m;
        % cd73 ring
        vmat(ii, 12) = mean(X_(:, 6, 2) > ds(nn).cd73_thr);
        vmat(ii, 13) = trapz(precx, log(X(:, 6, 2)+1)/log(2)) / 100;
        % ppat
        vmat(ii, 14) = mean(X_(:, 7, 2) > ds(nn).ppat_thr);
        vmat(ii, 15) = trapz(precx, log(X(:, 7, 2)+1)/log(2)) / 100;
        % nt5c2
        vmat(ii, 16) = mean(X_(:, 8, 2) > ds(nn).nt5c2_thr);
        vmat(ii, 17) = trapz(precx, log(X(:, 8, 2)+1)/log(2)) / 100;
        % adk
        vmat(ii, 18) = mean(X_(:, 9, 2) > ds(nn).adk_thr);
        vmat(ii, 19) = trapz(precx, log(X(:, 9, 2)+1)/log(2)) / 100;
    end
    vmat = array2table(vmat, 'VariableNames', ...
        {'dapi_area', 'neun_thr_ratio', 'neun_median', ...
        'gfap_filt_hole_ratio', 'gfap_filt_hole_median', ...
        'gfap_filt_ring_ratio', 'gfap_filt_ring_median', ...
        'gfap_th_ring_ratio', 'gfap_th_ring_median', ...
        'adora1_hole_ratio', 'adora1_hole_median', ...
        'cd73_ring_ratio', 'cd73_ring_median', ...
        'ppat_ring_ratio', 'ppat_ring_median', ...
        'nt5c2_ring_ratio', 'nt5c2_ring_median', ...
        'adk_ring_ratio', 'adk_ring_median', ...
        });

    ds(nn).vmat = vmat;
    fprintf('%d finish in %gs.\n', nn, toc(t0));
end

%% identify neuron & astrocyte v3
for nn = 1:ns
    vmat = ds(nn).vmat;
    is_dapi = vmat.dapi_area >= 256 & vmat.dapi_area <= 1024;

    v1 = log(vmat.neun_median+1)/log(2);
    v3 = vmat.gfap_th_ring_median;
    v2 = vmat.gfap_th_ring_ratio;

    rng(4488);
    gm = fitgmdist(v1(is_dapi & v1 > 0), 2);
    [mu, oo] = sort(gm.mu); %mu'
    sigma = gm.Sigma(oo);
    gx = (mu(1):0.05:mu(2))';
    gy = pdf(gm, gx);
    [~, idx] = min(gy);
    % fprintf('%d %.2f (%.2f), %.2f (%.2f)\n', nn, mu(1)+sigma(1)*5, 2^(mu(1)+sigma(1)*5)-1, gx(idx), 2^(gx(idx))-1);
    thr1 = gx(idx);

    gm = fitgmdist(v2(is_dapi & v1 <= thr1 & v2 > 0), 2);
    [mu, oo] = sort(gm.mu); %mu'
    sigma = gm.Sigma(oo);
    gx = (mu(1):0.05:mu(2))';
    gy = pdf(gm, gx);
    [~, idx] = min(gy);
    % fprintf('%d %.2f (%.2f), %.2f (%.2f)\n', nn, mu(1)+sigma(1)*5, 2^(mu(1)+sigma(1)*5)-1, gx(idx), 2^(gx(idx))-1);
    thr2 = gx(idx);

    X = v1(is_dapi & v1 > thr1 & v2 <= thr2);
    thr3 = median(X) + mad(X, 1) * 1.4826*1;

    is_neun = is_dapi & v1 > thr3 & v2 <= thr2;
    is_astro = is_dapi & v2 > thr2 & v1 <= thr1;

    % [mu, oo] = sort(gm.mu); mu'
    % sigma = gm.Sigma(oo);
    % gx = (mu(1):0.01:mu(3))';
    % gy = pdf(gm, gx);
    % [~, idx] = min(gy);
    figure; set(gcf, 'Position', [150, 200, 1200, 800]); hold on;
    subplot(2, 3, 1); hold on;
    histogram(v1(is_dapi), vx_log2);
    xline(thr1);
    subplot(2, 3, 4); hold on; 
    histogram(v2(is_dapi), (0:0.02:1));
    histogram(v2(is_dapi & v1 <= thr1 & v2 > 0), (0:0.02:1));
    xline(thr2);
    subplot(2, 3, [2, 3, 5, 6]);    
    scatter(v1(is_dapi), v2(is_dapi), 5, round(v3(is_dapi) * 20));
    yline(thr2);
    xline(thr1);
    xline(thr3);

    ds(nn).is_dapi = is_dapi;
    ds(nn).is_neun = is_neun;
    ds(nn).is_astro = is_astro;
    ds(nn).neun_thr = 2.^[thr1, thr3] -1;
    ds(nn).gfap_th_ratio_thr = thr2;
end

%%
outname = 'neun_gfap_fraction_all_thr';
cmap = tab10;
sname = arrayfun(@(x)(x.fid), ds, 'UniformOutput', false);
np = 4;
fd_name = { ...
    'adora1_hole_ratio', 'cd73_ring_ratio', 'ppat_ring_ratio', ...
    'nt5c2_ring_ratio', 'adk_ring_ratio'};
fd_label = {'ADORA1', 'CD73', 'PPAT', 'NT5C2', 'ADK'};
grp_idx = {[2, 7, 8], [4, 6, 9], [1, 3, 5]};
grp_label = {'Ctrl', 'TLE', 'FCD'};
grp_color = {[0.6, 0.6, 0.6], [0.2, 0.7, 0.8], [0.9, 0.4, 0.8]};
ngrp = length(grp_idx);
rng(76);

figure; set(gcf, 'Position', [50, 50, 780, 780]);
for ff = 1:length(fd_name)
    vmat = cell(ngrp, 1);
    [vm, vsd] = deal(nan(ngrp, 2));
    for gg = 1:ngrp
        mv = cell(np, 1);
        for kk = 1:np
            tmp = arrayfun(@(x)([ ...
                mean(x.vmat.(fd_name{ff})(x.is_neun & ismember(x.tile_idx, (kk:np:9000)')) > 0.1), ...
                mean(x.vmat.(fd_name{ff})(x.is_astro & ismember(x.tile_idx, (kk:np:9000)')) > 0.1)] * 100), ...
                ds(grp_idx{gg}), 'UniformOutput', false);
            mv{kk} = cat(1, tmp{:});
        end
        vmat{gg} = cat(1, mv{:});
        vm(gg, :) = mean(vmat{gg}, 1);
        vsd(gg, :) = std(vmat{gg}, 0, 1)/sqrt(size(vmat{gg}, 1));
    end
    s = arrayExportStats(vmat', {'NeuN+', 'GFAP+'}, grp_label, 'Fraction', ...
        'test.xlsx', 'manova');

    subplot('Position', [ff * 0.2 - 0.15, 0.5, 0.14, 0.14]); hold on; 
    h = bar([1, 2], vm', 'BarWidth', 0.8);
    for gg = 1:ngrp
        set(h(gg), 'EdgeColor', grp_color{gg}, 'FaceColor', 'w', ...
            'LineWidth', 1.2);
        errorbar(h(gg).XEndPoints, vm(gg, :), vsd(gg, :), '.', ...
            'Color', grp_color{gg}, 'CapSize', 5, 'LineWidth', 1.2, ...
            'MarkerSize', 1);
        plot(h(gg).XEndPoints + (rand(size(vmat{gg}, 1), 1) - 0.5) * 0.1, ...
            vmat{gg}, 'o', 'MarkerSize', 3, 'LineWidth', 1.2, ...
            'MarkerEdgeColor', grp_color{gg});
    end
    ylim = max(cat(1, vmat{:}), [], 1);

    p = s.groups.pValue(5);
    if p < 0.001; p = '***'; elseif p < 0.01; p = '**'; elseif p < 0.05; p = '*'; else; p = ''; end
    if ~isempty(p)
        plot([h(1).XEndPoints(1), h(1).XEndPoints(1), ...
            h(2).XEndPoints(1), h(2).XEndPoints(1)]-0.05, ...
            [0, 0.03, 0.03, 0]*max(ylim)+max(ylim)*0.05+ylim(1), '-k', 'LineWidth', 1);
        text((h(1).XEndPoints(1)+h(2).XEndPoints(1))/2-0.05, ...
            max(ylim)*0.07+ylim(1), p, ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'baseline', 'FontSize', 8);
    end

    p = s.groups.pValue(4);
    if p < 0.001; p = '***'; elseif p < 0.01; p = '**'; elseif p < 0.05; p = '*'; else; p = ''; end
    if ~isempty(p)
        plot([h(2).XEndPoints(1), h(2).XEndPoints(1), ...
            h(3).XEndPoints(1), h(3).XEndPoints(1)]+0.05, ...
            [0, 0.03, 0.03, 0]*max(ylim)+max(ylim)*0.05+ylim(1), '-k', 'LineWidth', 1);
        text((h(2).XEndPoints(1)+h(3).XEndPoints(1))/2+0.05, ...
            max(ylim)*0.07+ylim(1), p, ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'baseline', 'FontSize', 8);
    end

    p = s.groups.pValue(3);
    if p < 0.001; p = '***'; elseif p < 0.01; p = '**'; elseif p < 0.05; p = '*'; else; p = ''; end
    if ~isempty(p)
        plot([h(1).XEndPoints(1), h(1).XEndPoints(1), ...
            h(3).XEndPoints(1), h(3).XEndPoints(1)], ...
            [0, 0.03, 0.03, 0]*max(ylim)+max(ylim)*0.18+ylim(1), '-k', 'LineWidth', 1);
        text((h(1).XEndPoints(1)+h(3).XEndPoints(1))/2, ...
            max(ylim)*0.20+ylim(1), p, ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'baseline', 'FontSize', 8);
    end

    p = s.groups.pValue(11);
    if p < 0.001; p = '***'; elseif p < 0.01; p = '**'; elseif p < 0.05; p = '*'; else; p = ''; end
    if ~isempty(p)
        plot([h(1).XEndPoints(2), h(1).XEndPoints(2), ...
            h(2).XEndPoints(2), h(2).XEndPoints(2)]-0.05, ...
            [0, 0.03, 0.03, 0]*max(ylim)+max(ylim)*0.05+ylim(2), '-k', 'LineWidth', 1);
        text((h(1).XEndPoints(2)+h(2).XEndPoints(2))/2-0.05, ...
            max(ylim)*0.07+ylim(2), p, ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'baseline', 'FontSize', 8);
    end

    p = s.groups.pValue(10);
    if p < 0.001; p = '***'; elseif p < 0.01; p = '**'; elseif p < 0.05; p = '*'; else; p = ''; end
    if ~isempty(p)
        plot([h(2).XEndPoints(2), h(2).XEndPoints(2), ...
            h(3).XEndPoints(2), h(3).XEndPoints(2)]+0.05, ...
            [0, 0.03, 0.03, 0]*max(ylim)+max(ylim)*0.05+ylim(2), '-k', 'LineWidth', 1);
        text((h(2).XEndPoints(2)+h(3).XEndPoints(2))/2+0.05, ...
            max(ylim)*0.07+ylim(2), p, ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'baseline', 'FontSize', 8);
    end

    p = s.groups.pValue(9);
    if p < 0.001; p = '***'; elseif p < 0.01; p = '**'; elseif p < 0.05; p = '*'; else; p = ''; end
    if ~isempty(p)
        plot([h(1).XEndPoints(2), h(1).XEndPoints(2), ...
            h(3).XEndPoints(2), h(3).XEndPoints(2)], ...
            [0, 0.03, 0.03, 0]*max(ylim)+max(ylim)*0.18+ylim(2), '-k', 'LineWidth', 1);
        text((h(1).XEndPoints(2)+h(3).XEndPoints(2))/2, ...
            max(ylim)*0.20+ylim(2), p, ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'baseline', 'FontSize', 8);
    end

    set(gca, 'YLim', [0, max(ylim) * 1.25], 'XLim', [0.5, 2.5], ...
        'XTick', [1, 2], 'XTickLabels', {'NeuN+', 'GFAP+'});
    ylabel('Fraction (%)');
    title(fd_label{ff});

    if ff == 1; wm = 'replacefile'; else; wm = 'overwritesheet'; end
    writetable(s.values, [outname, '.xlsx'], 'Sheet', fd_label{ff}, ...
        'Range', 'A1', 'WriteMode', wm);
    writetable(s.main, [outname, '.xlsx'], 'Sheet', fd_label{ff}, ...
        'Range', 'A58');
    writetable(s.groups, [outname, '.xlsx'], 'Sheet', fd_label{ff}, ...
        'Range', 'A68');
end   
h = legend(grp_label, 'Box', 'off', 'FontSize', 9, 'NumColumns', 3);
h.Position = [0.6, 0.7, 0.2, 0.02];

print(gcf, outname, '-djpeg', '-r300');
print(gcf, outname, '-dpdf', '-r300', '-vector');

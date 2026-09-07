%%
chunk_a = 1024;
vx_log2 = (0:0.05:8)';
vx = 2 .^ vx_log2 - 1;
precx = (0:5:100)';

nv = length(vx);
ns = length(ds);

im_fun = @(f, h1, h2, w1, w2)(cat(3, ...
    imread(f, 1, 'PixelRegion', {[h1, h2], [w1, w2]}), ...
    imread(f, 3, 'PixelRegion', {[h1, h2], [w1, w2]})));

%% set threshold according slice distribution
for nn = 1:ns
    isin = unique(ds(nn).tile_idx);

    % % D1: 2gm(c2, m+1sd), D2: m+1sd
    % dx = mean(ds(nn).chan_dist.neun_filt(isin, :), 1);
    % X = arrayfun(@(x)(ones(round(dx(x)), 1) * vx_log2(x)), ...
    %     2:nv, 'UniformOutput', false);
    % X = cat(1, X{:});
    % ds(nn).neun_thr = 2 ^ (median(X) + mad(X, 1) * 1.4826 * 1) - 1;

    % D1: 2gm(c2, m+1sd), D2: m+1sd
    dx = mean(ds(nn).chan_dist.adora1_filt(isin, :), 1);
    X = arrayfun(@(x)(ones(round(dx(x)), 1) * vx_log2(x)), ...
        2:nv, 'UniformOutput', false);
    X = cat(1, X{:});
    ds(nn).adora1_thr = 2 ^ (median(X) + mad(X, 1) * 1.4826 * 1) - 1;

    % D1: 2gm(dmin), D2: 2gm(c2, m+1sd)
    dx = mean(ds(nn).chan_dist.gcsh_filt(isin, :), 1);
    X = arrayfun(@(x)(ones(round(dx(x)), 1) * vx_log2(x)), ...
        2:nv, 'UniformOutput', false);
    X = cat(1, X{:});
    rng(6535)
    gm = fitgmdist(X, 2);
    [mu, oo] = sort(gm.mu); 
    gx = (mu(1):0.05:mu(2))';
    gy = pdf(gm, gx);
    [~, idx] = min(gy);
    ds(nn).gcsh_thr = 2 ^ (gx(idx)) - 1;
    ds(nn).gcsh_log2_m = gx(idx);

    % D1: m+1sd, D2: 2gm(c2, m+1sd)
    dx = mean(ds(nn).chan_dist.rtn4_filt(isin, :), 1);
    X = arrayfun(@(x)(ones(round(dx(x)), 1) * vx_log2(x)), ...
        2:nv, 'UniformOutput', false);
    X = cat(1, X{:});
    ds(nn).rtn4_thr = 2 ^ (median(X) + mad(X, 1) * 1.4826 * 1) - 1;
    ds(nn).rtn4_log2_m = median(X);

    % D1: m+1sd, D2: 2gm(c2, m+1sd)
    dx = mean(ds(nn).chan_dist.bmper_filt(isin, :), 1);
    X = arrayfun(@(x)(ones(round(dx(x)), 1) * vx_log2(x)), ...
        2:nv, 'UniformOutput', false);
    X = cat(1, X{:});
    ds(nn).bmper_thr = 2 ^ (median(X) + mad(X, 1) * 1.4826 * 1) - 1;
    ds(nn).bmper_log2_m = median(X);

    % rng(6535)
    % gm = fitgmdist(X, 2);
    % [mu, oo] = sort(gm.mu); %mu'
    % sigma = gm.Sigma(oo);
    % gx = (mu(1):0.05:mu(2))';
    % gy = pdf(gm, gx);
    % [~, idx] = min(gy);
    % % mu(2) = median(X);
    % % sigma(2) = mad(X, 1) * 1.4826;
    % fprintf('%d %.2f (%.2f), thr: %.2f (%.2f).\n', nn, ...
    %     median(X), mad(X, 1) * 1.4826, mu(2), sigma(2));
    % figure; histogram(X, vx_log2);
    % xline(mu, ':k');
    % xline(mu(2) + sigma(2) * 1, ':r');
    % xline(median(X) + mad(X, 1) * 1.4826 * 1, ':m');
end

%% metrics
precx_hi = (0:0.1:100)';
precx_md = precx == 50;
for nn = 1:ns
    t0 = tic;
    nc = length(ds(nn).cc_dapi.PixelIdxList);

    vmat = zeros(nc, 9);
    for ii = 1:nc
        vmat(ii, 1) = ds(nn).circ_area{ii}(1);
        Ar = ds(nn).circ_area{ii} / sum(ds(nn).circ_area{ii});
        X = ds(nn).misc_metrics{ii};
        X_ = interp1(precx, X, precx_hi, 'spline');
        % neun
        vmat(ii, 2) = X(precx_md, 2, 1);
        % adora1 hole
        vmat(ii, 3) = X(precx_md, 3, 1);
        % gcsh ring
        vmat(ii, 4) = mean(X_(:, 4, 2) > ds(nn).gcsh_thr);
        vmat(ii, 5) = trapz(precx, log(X(:, 4, 2)+1)/log(2)) / 100;
        % rtn4 ring
        vmat(ii, 6) = mean(X_(:, 5, 2) > ds(nn).rtn4_thr);
        vmat(ii, 7) = trapz(precx, log(X(:, 5, 2)+1)/log(2)) / 100;
        % bmper
        vmat(ii, 8) = Ar * squeeze(mean(X_(:, 6, :) > ds(nn).bmper_thr));
        vmat(ii, 9) = Ar * squeeze(trapz(precx, log(X(:, 6, :)+1)/log(2))) / 100;
    end
    vmat = array2table(vmat, 'VariableNames', ...
        {'dapi_area', 'neun_median', 'adora1_hole_median', ...
        'gcsh_ring_ratio', 'gcsh_ring_median', ...
        'rtn4_ring_ratio', 'rtn4_ring_median', ...
        'bmper_ring_ratio', 'bmper_ring_median', ...
        });

    ds(nn).vmat = vmat;
    fprintf('%d finish in %gs.\n', nn, toc(t0));
end

%% identify neuron & adora1 high/low v1
for nn = 1:ns
    vmat = ds(nn).vmat;
    is_dapi = vmat.dapi_area >= 256 & vmat.dapi_area <= 1024;
 
    rng(4488);
    X = log(vmat.neun_median(is_dapi & vmat.neun_median > 1 & vmat.adora1_hole_median > 7)+1)/log(2);
    neun_gm = fitgmdist(X, 1);
    [~, imax] = max(neun_gm.mu);
    neun_thr = 2^(neun_gm.mu(imax) + neun_gm.Sigma(imax))-1;
    is_neun = is_dapi & vmat.neun_median > neun_thr & vmat.adora1_hole_median > 7;

    X = log(vmat.adora1_hole_median(is_neun)+1)/log(2);
    adora1_gm = fitgmdist(X, 2);
    gx = (min(adora1_gm.mu):0.01:max(adora1_gm.mu))';
    gy = pdf(adora1_gm, gx);
    [~, imin] = min(gy);
    adora1_thr = 2^gx(imin)-1;
    is_adora1_hi = is_neun & vmat.adora1_hole_median > adora1_thr;
    is_adora1_lo = is_neun & vmat.adora1_hole_median <= adora1_thr;

    figure; set(gcf, 'Position', [50, 50, 1000, 780]);
    subplot(3, 4, [5, 6, 9, 10]);
    plot(vmat.dapi_area(is_dapi)+1, vmat.neun_median(is_dapi)+1, '.');
    % xline(adora1_thr, ':k');
    yline(neun_thr, ':k');
    set(gca, 'TickDir', 'out', 'XScale', 'log', 'YScale', 'log', 'YLim', [0, 256]);
    xlabel('Adora1 median');  
    ylabel('NeuN median');

    subplot(3, 4, [7, 8, 11, 12]);
    plot(vmat.adora1_hole_median(is_dapi)+1, vmat.neun_median(is_dapi)+1, '.');
    xline(adora1_thr, ':k');
    yline(neun_thr, ':k');
    set(gca, 'TickDir', 'out', 'XScale', 'log', 'YScale', 'log', 'YLim', [0, 256]);
    xlabel('Adora1 median');  
    ylabel('NeuN median');

    ds(nn).is_dapi = is_dapi;
    ds(nn).is_neun = is_neun;
    ds(nn).is_adora1_hi = is_adora1_hi;
    ds(nn).is_adora1_lo = is_adora1_lo;
end

%%
outname = 'test.xlsx';
cmap = tab10;
sname = arrayfun(@(x)(x.fid), ds, 'UniformOutput', false);
figure; set(gcf, 'Position', [50, 50, 780, 780]);

subplot('Position', [0.1, 0.5, 0.14, 0.14]); hold on;
vmat = arrayfun(@(x)([ ...
    mean(x.vmat.gcsh_ring_ratio(x.is_adora1_lo)), ...
    mean(x.vmat.gcsh_ring_ratio(x.is_adora1_hi))]), ...
    ds, 'UniformOutput', false);
vmat = cat(1, vmat{:});
[~, p] = ttest(diff(vmat, 1, 2));
for ii = 1:length(ds)
plot(vmat(ii, :), 'Color', cmap(ii, :));
end
text(1.5, max(get(gca, 'YLim')), sprintf('p = %.3f', p), 'FontSize', 8, ...
    'HorizontalAlignment', 'center');
set(gca, 'XLim', [0.7, 2.3], 'XTick', [1, 2], 'XTickLabel', {'Low', 'High'});
xlabel('Adora1');
ylabel('Fraction (%)');
title('GCSH');
writetable(table(sname, vmat(:, 1), vmat(:, 2), 'VariableNames', {'Slice', 'Low', 'High'}), ...
    outname, 'Sheet', 'GCSH', 'WriteMode', 'replacefile');

subplot('Position', [0.3, 0.5, 0.14, 0.14]); hold on;
vmat = arrayfun(@(x)([ ...
    mean(x.vmat.rtn4_ring_median(x.is_adora1_lo)), ...
    mean(x.vmat.rtn4_ring_median(x.is_adora1_hi))]/100), ...
    ds, 'UniformOutput', false);
vmat = cat(1, vmat{:});
[~, p] = ttest(diff(vmat, 1, 2));
for ii = 1:length(ds)
plot(vmat(ii, :), 'Color', cmap(ii, :));
end
text(1.5, max(get(gca, 'YLim')), sprintf('p = %.3f', p), 'FontSize', 8, ...
    'HorizontalAlignment', 'center');
set(gca, 'XLim', [0.7, 2.3], 'XTick', [1, 2], 'XTickLabel', {'Low', 'High'});
xlabel('Adora1');
ylabel('Fraction (%)');
title('RTN4');
writetable(table(sname, vmat(:, 1), vmat(:, 2), 'VariableNames', {'Slice', 'Low', 'High'}), ...
    outname, 'Sheet', 'RTN4');

subplot('Position', [0.5, 0.5, 0.14, 0.14]); hold on;
vmat = arrayfun(@(x)([ ...
    mean(x.vmat.bmper_ring_median(x.is_adora1_lo)), ...
    mean(x.vmat.bmper_ring_median(x.is_adora1_hi))]), ...
    ds, 'UniformOutput', false);
vmat = cat(1, vmat{:});
[~, p] = ttest(diff(vmat, 1, 2));
for ii = 1:length(ds)
plot(vmat(ii, :), 'Color', cmap(ii, :));
end
text(1.5, max(get(gca, 'YLim')), sprintf('p = %.3f', p), 'FontSize', 8, ...
    'HorizontalAlignment', 'center');
set(gca, 'XLim', [0.7, 2.3], 'XTick', [1, 2], 'XTickLabel', {'Low', 'High'});
xlabel('Adora1');
ylabel('Fraction (%)');
title('BMPER');
writetable(table(sname, vmat(:, 1), vmat(:, 2), 'VariableNames', {'Slice', 'Low', 'High'}), ...
    outname, 'Sheet', 'BMPER');

print(gcf, strrep(outname, '.xlsx', ''), '-djpeg');


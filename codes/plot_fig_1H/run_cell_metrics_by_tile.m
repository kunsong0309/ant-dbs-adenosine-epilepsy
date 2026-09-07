%%
outpath = 'results';
rawpaths = dir('..\20260312 上海 纪小帅 9张荧光7标扫描（9-9-9) 黄中庆 #SWR03106636\*.svs');
%%
chunk_a = 1024;
vx_log2 = (0:0.05:8)';
vx = 2 .^ vx_log2 - 1;

nv = length(vx);
ns = length(rawpaths);

im_fun = @(f, h1, h2, w1, w2)(cat(3, ...
    imread(f, 1, 'PixelRegion', {[h1, h2], [w1, w2]}), ...
    imread(f, 3, 'PixelRegion', {[h1, h2], [w1, w2]}), ...
    imread(f, 4, 'PixelRegion', {[h1, h2], [w1, w2]})));

ds = struct;
%%
for nn = 1:ns
    t0 = tic;
    fprintf('[%d] is processing >>>>>>>>>\n', nn);
    rp = fullfile(rawpaths(nn).folder, rawpaths(nn).name);
    fid = regexp(rawpaths(nn).name, '^[0-9-]+(?= )', 'match', 'once');
    fp = fullfile(outpath, [fid, '.mat']);
    % load(fp, 'cc_dapi', 'chan_dist');

    info = imfinfo(rp);
    w = info(1).Width;
    h = info(1).Height;

    [x1, y1] = meshgrid( ...
        (1:floor(w/chunk_a))*chunk_a-chunk_a/2, ...
        (1:floor(h/chunk_a))*chunk_a-chunk_a/2);
    [x2, y2] = meshgrid( ...
        (1:floor((w-chunk_a)/chunk_a))*chunk_a, ...
        (1:floor((h-chunk_a)/chunk_a))*chunk_a);
    tilex = cat(1, x1(:), x2(:));
    tiley = cat(1, y1(:), y2(:));
    [tilex, oo] = sort(tilex);
    tiley = tiley(oo);
    nt = length(tilex);

    t = tic;
    mask_dapi = tileProcess(rp, ...
        @(f, h1, h2, w1, w2)(imread(f, 1, 'PixelRegion', {[h1, h2], [w1, w2]})), ...
        @(m)(chanMaskCell(m(:, :, 1), 1.5)), ...
        [h, w], 0, 0, chunk_a, 16);
    fprintf('\t>>>>>> round 1 in %gs.\n', toc(t));

    mask_dapi = mask_dapi | tileProcess(rp, ...
        @(f, h1, h2, w1, w2)(imread(f, 1, 'PixelRegion', {[h1, h2], [w1, w2]})), ...
        @(m)(chanMaskCell(m(:, :, 1), 1.5)), ...
        [h, w], chunk_a/2, chunk_a/2, chunk_a, 16);
    fprintf('\t>>>>>> round 2 in %gs.\n', toc(t));

    cc_dapi = bwconncomp(mask_dapi, 8); % B = bwboundaries(dapi_mask, 'noholes');

    nc = length(cc_dapi.PixelIdxList);
    tile_idx = zeros(nc, 1);
    cell_cnt = zeros(nc, 2);
    for kk = 1:nc
        [ch, cw] = ind2sub([h, w], cc_dapi.PixelIdxList{kk});
        hc = median(ch);
        wc = median(cw);
        [~, tile_idx(kk)] = min((hc - tiley).^2 + (wc - tilex).^2);
        cell_cnt(kk, :) = [hc, wc];
    end
    fprintf('finish identify dapi in %gs.\n', toc(t));

    t = tic;
    cdist = struct;
    parfor (tt = 1:nt, 16)
        h1 = tiley(tt) - chunk_a / 2 + 1;
        h2 = tiley(tt) + chunk_a / 2;
        w1 = tilex(tt) - chunk_a / 2 + 1;
        w2 = tilex(tt) + chunk_a / 2;
        m = im_fun(rp, h1, h2, w1, w2);
        m = single(m);

        tmp = m(:, :, 8);
        cdist(tt).neun = histcounts(tmp, [vx; 256]);
        cdist(tt).neun_filt = histcounts(imgaussfilt(tmp, 5), [vx; 256]);

        tmp = m(:, :, 7);
        cdist(tt).gfap = histcounts(tmp, [vx; 256]);
        tmp = imgaussfilt(tmp, 3);
        cdist(tt).gfap_filt = histcounts(tmp, [vx; 256]);        
        cdist(tt).gfap_th = histcounts(imtophat(tmp, strel('disk', 8)), [vx; 256]);

        tmp = m(:, :, 6);
        cdist(tt).adora1 = histcounts(tmp, [vx; 256]);
        cdist(tt).adora1_filt = histcounts(imgaussfilt(tmp, 5), [vx; 256]);

        tmp = m(:, :, 2);
        cdist(tt).cd73 = histcounts(tmp, [vx; 256]);
        cdist(tt).cd73_filt = histcounts(imgaussfilt(tmp, 3), [vx; 256]);

        tmp = m(:, :, 3);
        cdist(tt).ppat = histcounts(tmp, [vx; 256]);
        cdist(tt).ppat_filt = histcounts(imgaussfilt(tmp, 3), [vx; 256]);

        tmp = m(:, :, 4);
        cdist(tt).nt5c2 = histcounts(tmp, [vx; 256]);
        cdist(tt).nt5c2_filt = histcounts(imgaussfilt(tmp, 3), [vx; 256]);

        tmp = m(:, :, 5);
        cdist(tt).adk = histcounts(tmp, [vx; 256]);
        cdist(tt).adk_filt = histcounts(imgaussfilt(tmp, 5), [vx; 256]);

        if mod(tt, 500) == 0
            fprintf('\t>>>>> %d tiles in %gs.\n', tt, toc(t));
        end
    end

    fds = fields(cdist);
    chan_dist = struct;
    for jj = 1:length(fds)
        chan_dist.(fds{jj}) = cat(1, cdist(:).(fds{jj}));
    end
    fprintf('finish channel distribution in %gs.\n', toc(t));

    ds(nn, 1).fid = fid;
    ds(nn).rp = rp;
    ds(nn).fp = fp;
    ds(nn).im_size = [h, w];
    ds(nn).tilex = tilex;
    ds(nn).tiley = tiley;
    ds(nn).tile_idx = tile_idx;
    ds(nn).cell_cnt = cell_cnt;
    ds(nn).chan_dist = chan_dist;
    ds(nn).cc_dapi = cc_dapi;

    save(fp, 'cc_dapi', 'chan_dist', '-v7.3');
    fprintf('[%d] finish in %gs.\n', nn, toc(t0));
end

%%
for nn = 1:ns
    figure; set(gcf, 'Position', [50, 50, 1200, 600]);
    isin = unique(ds(nn).tile_idx);

    subplot(2, 4, 1); hold on;
    plot(vx_log2, mean(ds(nn).chan_dist.neun(isin, :)));
    plot(vx_log2, mean(ds(nn).chan_dist.neun_filt(isin, :)));
    set(gca, 'XLim', [0, 8], 'YScale', 'log', 'YLim', [1e-3, 1e6]); grid on;
    title('NeuN');

    subplot(2, 4, 2); hold on;
    plot(vx_log2, mean(ds(nn).chan_dist.gfap(isin, :)));
    plot(vx_log2, mean(ds(nn).chan_dist.gfap_filt(isin, :)));
    plot(vx_log2, mean(ds(nn).chan_dist.gfap_th(isin, :)));
    set(gca, 'XLim', [0, 8], 'YScale', 'log', 'YLim', [1e-3, 1e6]); grid on;
    title('GFAP');

    subplot(2, 4, 3); hold on;
    plot(vx_log2, mean(ds(nn).chan_dist.adora1(isin, :)));
    plot(vx_log2, mean(ds(nn).chan_dist.adora1_filt(isin, :)));
    set(gca, 'XLim', [0, 8], 'YScale', 'log', 'YLim', [1e-3, 1e6]); grid on;
    title('ADORA1');

    subplot(2, 4, 4); hold on;
    plot(vx_log2, mean(ds(nn).chan_dist.cd73(isin, :)));
    plot(vx_log2, mean(ds(nn).chan_dist.cd73_filt(isin, :)));
    set(gca, 'XLim', [0, 8], 'YScale', 'log', 'YLim', [1e-3, 1e6]); grid on;
    title('CD73');

    subplot(2, 4, 5); hold on;
    plot(vx_log2, mean(ds(nn).chan_dist.ppat(isin, :)));
    plot(vx_log2, mean(ds(nn).chan_dist.ppat_filt(isin, :)));
    set(gca, 'XLim', [0, 8], 'YScale', 'log', 'YLim', [1e-3, 1e6]); grid on;
    title('PPAT');

    subplot(2, 4, 6); hold on;
    plot(vx_log2, mean(ds(nn).chan_dist.nt5c2(isin, :)));
    plot(vx_log2, mean(ds(nn).chan_dist.nt5c2_filt(isin, :)));
    set(gca, 'XLim', [0, 8], 'YScale', 'log', 'YLim', [1e-3, 1e6]); grid on;
    title('NT5C2');

    subplot(2, 4, 7); hold on;
    plot(vx_log2, mean(ds(nn).chan_dist.adk(isin, :)));
    plot(vx_log2, mean(ds(nn).chan_dist.adk_filt(isin, :)));
    set(gca, 'XLim', [0, 8], 'YScale', 'log', 'YLim', [1e-3, 1e6]); grid on;
    title('ADK');

    sgtitle(ds(nn).fid);
    print(gcf, fullfile(outpath, 'chan_dist', ds(nn).fid), '-djpeg', '-r300');
end

%%
precx = (0:5:100)';
dom_rad = 10;

for nn = 1:ns
    t = tic;
    fprintf('%d processing %s >>>>>>>\n', nn, ds(nn).fid);
    rp = ds(nn).rp;
    h = ds(nn).im_size(1);
    w = ds(nn).im_size(2);
    tilex = ds(nn).tilex;
    tiley = ds(nn).tiley;
    nt = length(ds(nn).tilex);
    nc = length(ds(nn).tile_idx);

    mask_dapi = false(ds(nn).im_size);
    mask_dapi(cat(1, ds(nn).cc_dapi.PixelIdxList{:})) = true;

    ds(nn).misc_metrics = cell(nc, 1);
    ds(nn).circ_area = cell(nc, 1);

    t0 = 0;
    for tt = 1:50:nt
        cidx = find(ismember(ds(nn).tile_idx, tt + (0:49)));
        % cidx = find(ds(nn).tile_idx == tt);
        ntc = length(cidx);
        if ntc == 0; continue; end
        misc_metrics = cell(ntc, 1);
        circ_area = cell(ntc, 1);
        cell_pixels = ds(nn).cc_dapi.PixelIdxList(cidx);

        h1 = tiley(tt) - chunk_a / 2 + 1;
        h2 = tiley(tt) + chunk_a / 2;
        w1 = tilex(tt) - chunk_a / 2 + 1;
        w2 = tilex(tt) + chunk_a / 2;
        b = mask_dapi(h1:h2, w1:w2);
        m = single(im_fun(rp, h1, h2, w1, w2));

        parfor (ii = 1:ntc, 32)
        % for ii = 1:ntc
            [ch, cw] = ind2sub([h, w], cell_pixels{ii});
            h1 = min(max(min(ch) - dom_rad * 2, 1), h);
            h2 = min(max(max(ch) + dom_rad * 2, 1), h);
            w1 = min(max(min(cw) - dom_rad * 2, 1), w);
            w2 = min(max(max(cw) + dom_rad * 2, 1), w);
            cbi = sub2ind([h2-h1+1, w2-w1+1], ch-h1+1, cw-w1+1);
            cb = false(h2-h1+1, w2-w1+1);
            cb(cbi) = true;

            cbg = imdilate(cb, strel('disk', dom_rad)) & ~cb;
            cbgi = find(cbg);

            cm = single(im_fun(rp, h1, h2, w1, w2));
            cm = cat(3, ...
                imgaussfilt(cm(:, :, 1), 5), ...
                imgaussfilt(cm(:, :, 2), 3), ...
                imgaussfilt(cm(:, :, 3), 3), ...
                imgaussfilt(cm(:, :, 4), 3), ...
                imgaussfilt(cm(:, :, 5), 5), ...
                imgaussfilt(cm(:, :, 6), 5), ...
                imgaussfilt(cm(:, :, 7), 3), ...
                imgaussfilt(cm(:, :, 8), 5));
            cm = cat(3, cm, imtophat(cm(:, :, 7), strel('disk', 8)));
            cm = cm(:, :, [1, 8, 7, 9, 6, 2, 3, 4, 5]);

            % [misc_metrics{ii}, circ_area{ii}] = ...
            %     cellFitBoundary(cm, cb, dom_rad);
            
            cm = reshape(cm, [], size(cm, 3));
            misc_metrics{ii} = cat(3, ...
                prctile(cm(cbi, :), precx, 1), ...
                prctile(cm(cbgi, :), precx, 1));
            circ_area{ii} = [length(cbi), length(cbgi)];
        end

        ds(nn).misc_metrics(cidx) = misc_metrics;
        ds(nn).circ_area(cidx) = circ_area;

        t1 = sum(cellfun(@(x)(~isempty(x)), ds(nn).misc_metrics));
        if t1 - t0 > 10000
            fprintf('finish %d in %gs.\n', t1, toc(t));
            t0 = t1;
        end
    end
    fprintf('%d finish in %gs.\n', nn, toc(t));
end

%%
out_dir = '..';
load('..\ds.mat', 'ds', 'fs');
ds(arrayfun(@(x)(~any(x.ich_ant(x.is_ant)) | ~any(x.ich_zone)), ds)) = [];
%%
h = design(fdesign.bandpass('N,F3dB1,F3dB2', 4, 0.1, 100, fs), 'butter');
for nn = 1:numel(ds)
    t0 = tic; 
    fn = sprintf('%s_%s', ds(nn).subject, ds(nn).fname);
    load(fullfile(out_dir, 'prep_sig', [fn, '.mat']), 'sig_ant', 'sig_zone');
    
    sig = sig_ant(:, ds(nn).is_ant);
    sig(isnan(sig)) = 0;
    sig = filtfilt(h.sosMatrix, h.ScaleValues, sig);
    chs = ds(nn).elec_ant(ds(nn).is_ant);
    [sxx, txx, fxx] = mtspecgramc(sig, [4, 4], ...
        struct('tapers', [3, 5], 'Fs', fs, 'fpass', [0, 50]));
    txx = txx / 60 / 60;
    nch = size(sig, 2);
    for kk = 1:nch
        figure; set(gcf, 'Position', [100, 100, 800, 300], 'Visible', 'off');
        imagesc(txx, fxx, 10 * log10(sxx(:, :, kk)'), [0, 40]);
        colormap(jet);
        ylabel(chs{kk});
        set(gca, 'YDir', 'normal', 'TickDir', 'out', 'XLim', [0, 4]);
        bstr = sprintf('%s (ANT) %s', fn, chs{kk});
        title(bstr, 'Interpreter', 'none');
        print(fullfile(out_dir, 'check_misc', bstr), '-djpeg');
        close;
    end

    sig = sig_zone;
    sig(isnan(sig)) = 0;
    sig = filtfilt(h.sosMatrix, h.ScaleValues, sig);
    chs = ds(nn).elec_zone;
    [sxx, txx, fxx] = mtspecgramc(sig, [4, 4], ...
        struct('tapers', [3, 5], 'Fs', fs, 'fpass', [0, 50]));
    txx = txx / 60 / 60;
    nch = size(sig, 2);
    for kk = 1:nch
        figure; set(gcf, 'Position', [100, 100, 800, 300], 'Visible', 'off');
        imagesc(txx, fxx, 10 * log10(sxx(:, :, kk)'), [0, 40]);
        colormap(jet);
        ylabel(chs{kk});
        set(gca, 'YDir', 'normal', 'TickDir', 'out', 'XLim', [0, 4]);
        bstr = sprintf('%s (Zone) %s', fn, chs{kk});
        title(bstr, 'Interpreter', 'none');
        print(fullfile(out_dir, 'check_misc', bstr), '-djpeg');
        close;
    end

    fprintf('[%.2d] %s >>> %gs.\n', nn, fn, toc(t0));
end


%% for IEDs
freq_low = 3;
freq_high = 100;
h = design(fdesign.bandpass('N,F3dB1,F3dB2', 4, freq_low, freq_high, fs), 'butter');

for nn = 1:length(ds)
    t0 = tic; 
    fn = sprintf('%s_%s', ds(nn).subject, ds(nn).fname);
    ts = (0:ds(nn).nt-1)' / fs;
    load(fullfile(out_dir, 'prep_sig', [fn, '.mat']), 'sig_ant', 'sig_zone');

    sig = sig_ant(:, ds(nn).is_ant);
    sig(isnan(sig)) = 0;
    sig = filtfilt(h.sosMatrix, h.ScaleValues, sig);
    sigma = ds(nn).ant_sigma;
    nch = size(sig, 2);
    ds(nn).ant_ev_ied = deal(cell(1, nch));
    for kk = 1:nch
        [~, ev_hamp] = eegDetectHighAmplitude(sig(:, kk), ts, [], sigma(kk), 1, [3, 2]);
        [~, ev_sspk] = eegDetectSharpSpikes(sig(:, kk), ts, [], sigma(kk), 1, 3);
        [~, ev_swd] = eegDetectSpikeWave(sig(:, kk), ts, [], sigma(kk), 1, [5, 3]);
        [~, ev_se] = eegCalibrateSeizure(sig(:, kk) / sigma(kk), ts, 1, ...
            ev_hamp, [-inf, inf], ev_sspk, [-inf, inf], ev_swd, [-inf, inf]);
        ds(nn).ant_ev_ied{kk} = [ts(ev_se(:, 1)), ts(ev_se(:, 2))];    
    end

    sig = sig_zone(:, ds(nn).is_zone);
    sig(isnan(sig)) = 0;
    sig = filtfilt(h.sosMatrix, h.ScaleValues, sig);
    sigma = ds(nn).soz_sigma(ds(nn).is_zone');
    nch = size(sig, 2);
    ds(nn).soz_ev_ied = deal(cell(1, nch));
    for kk = 1:nch
        [~, ev_hamp] = eegDetectHighAmplitude(sig(:, kk), ts, [], sigma(kk), 1, [3, 2]);
        [~, ev_sspk] = eegDetectSharpSpikes(sig(:, kk), ts, [], sigma(kk), 1, 3);
        [~, ev_swd] = eegDetectSpikeWave(sig(:, kk), ts, [], sigma(kk), 1, [5, 3]);
        [~, ev_se] = eegCalibrateSeizure(sig(:, kk) / sigma(kk), ts, 0.5, ...
            ev_hamp, [-inf, inf], ev_sspk, [-inf, inf], ev_swd, [-inf, inf]);
        ds(nn).soz_ev_ied{kk} = [ts(ev_se(:, 1)), ts(ev_se(:, 2))];    
    end

    fprintf('[%.2d] %s >>> %gs.\n', nn, fn, toc(t0));
end




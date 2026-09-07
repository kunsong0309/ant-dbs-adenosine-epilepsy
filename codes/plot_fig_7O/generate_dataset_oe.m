%%
workdir = 'prep_sig';
fdir = '24_sleep_seizure\WT_ANT_A1ROE_CA1_KA_20250122';
fls = dir(fullfile(fdir, 'M*'));
fls = fls(arrayfun(@(x)(x.isdir), fls));
%%
fs = 250;
h = design(fdesign.bandpass('N,F3dB1,F3dB2', 4, 0.1, 100, fs), 'butter');

ds = struct;
for nn = 1:numel(fls)
    to = tic;
    fprintf('[%.3d] %s >>> ', nn, fls(nn).name);
    mid = str2double(regexp(fls(nn).name, '(?<=M)\d+', 'match', 'once'));
    mls = dir(fullfile(fls(nn).folder, fls(nn).name, '*.nex5'));

    ds(nn, 1).gtype = 'WT';
    ds(nn).model = 'A1ROE';
    % ds(nn).batch = 1;
    % ds(nn).date = '20250117';
    ds(nn).batch = 2;
    ds(nn).date = '20250122';
    ds(nn).mouse = mid;
    ds(nn).ztshift = 0; 
    ds(nn).numblk = 18;
    ds(nn).fid = sprintf('%s_%s_batch#%d_%s_m%.1d', ...
        ds(nn).gtype, ds(nn).model, ds(nn).batch, ds(nn).date, ds(nn).mouse);

    nses = ds(nn).numblk;
    sig = cell(nses, 1);
    for ii = 1:nses
        da = readNex5File(fullfile(mls(ii).folder, mls(ii).name));
        chnames = cellfun(@(x)(lower(x.name)), da.contvars, 'UniformOutput', false);
        sfreq = da.freq;
        nt = round(fs * da.tend);
        da = cellfun(@(x)(single(x.data)), da.contvars, 'UniformOutput', false);
        dsx = sfreq / fs;

        [~, ich] = ismember({'eeg1_raw', 'eeg1', 'eeg_raw', 'eeg'}, chnames);
        ich = ich(find(ich, 1));
        if ~isempty(ich)
            eeg1 = filtfilt(h.sosMatrix, h.ScaleValues, da{ich}(1:dsx:end));
        else
            eeg1 = zeros(nt, 1);
        end

        [~, ich] = ismember({'eeg2_raw', 'eeg2'}, chnames);
        ich = ich(find(ich, 1));
        if ~isempty(ich)
            eeg2 = filtfilt(h.sosMatrix, h.ScaleValues, da{ich}(1:dsx:end));
        else
            eeg2 = zeros(nt, 1);
        end

        [~, ich] = ismember({'emg1_raw', 'emg1', 'emg_raw', 'emg'}, chnames);
        ich = ich(find(ich, 1));
        if ~isempty(ich)
            emg1 = filtfilt(h.sosMatrix, h.ScaleValues, da{ich}(1:dsx:end));
        else
            emg1 = zeros(nt, 1);
        end

        [~, ich] = ismember({'emg2_raw', 'emg2'}, chnames);
        ich = ich(find(ich, 1));
        if ~isempty(ich)
            emg2 = filtfilt(h.sosMatrix, h.ScaleValues, da{ich}(1:dsx:end));
        else
            emg2 = zeros(nt, 1);
        end

        sig{ii} = cat(2, eeg1, eeg2, emg1, emg2);
    end

    sig = cat(1, sig{:});
    scale = log10(max(mad(sig, 1, 1) * 1.4826));
    ds(nn).scale = scale;
    if scale < 1; sig = sig * 1e3; end

    [eeg1, eeg2, emg1, emg2] = deal(sig(:, 1), sig(:, 2), sig(:, 3), sig(:, 4));
    save(fullfile(workdir, 'prep_sig', [ds(nn).fid, '.mat']), ...
        'eeg1', 'eeg2', 'emg1', 'emg2', 'fs');
    fprintf('%gs.\n', toc(to));      
end

%%
for nn = 1:numel(ds)
    fp = fullfile(workdir, 'prep_sig', [ds(nn).fid, '.mat']);
    if ~exist(fp, 'file'); continue; end
    to = tic;
    fprintf('[%.3d] %s >>> ', nn, ds(nn).fid);
    load(fp, 'eeg1', 'eeg2', 'emg1', 'emg2', 'fs');

    [sxx, txx, fxx] = mtspecgramc([eeg1, eeg2], [4, 4], ...
        struct('tapers', [3, 5], 'Fs', fs, 'fpass', [0, 30]));
    h = design(fdesign.bandpass('N,F3dB1,F3dB2', 4, 10, 50, fs), 'butter');
    emg = filtfilt(h.sosMatrix, h.ScaleValues, [emg1, emg2]);
    ztxx = txx / 60 / 60 + ds(nn).ztshift;
    zts = (0:numel(emg1)-1)' / fs / 60 / 60 + ds(nn).ztshift;

    for ii = 1:ds(nn).numblk
        t1 = (ii - 1) * 4 + ds(nn).ztshift;
        t2 = t1 + 4;
        istxx = ztxx >= t1 & ztxx < t2;
        ists = zts >= t1 & zts < t2;  
        figure; set(gcf, 'Position', [50, 50, 1000, 800], 'Visible', 'off');
        
        subplot(4, 1, 1); hold on;
        imagesc(ztxx(istxx), fxx, 10 * log10(sxx(istxx, :, 1)'), [0, 40]);
        % plot(ds(nn).zts_se{1}' / 60 / 60 + ds(nn).ztshift, ...
        %     ds(nn).zts_se{1}' * 0 + 20, '-k', 'LineWidth', 1.5);
        colormap(jet); 
        ylabel('eeg1');
        set(gca, 'YDir', 'normal', 'TickDir', 'out', 'XLim', [t1, t2]);

        subplot(4, 1, 2); hold on;
        imagesc(ztxx(istxx), fxx, 10 * log10(sxx(istxx, :, 2)'), [0, 40]);
        % plot(ds(nn).zts_se{2}' / 60 / 60 + ds(nn).ztshift, ...
        %     ds(nn).zts_se{2}' * 0 + 20, '-k', 'LineWidth', 1.5);
        colormap(jet); 
        ylabel('eeg2');
        set(gca, 'YDir', 'normal', 'TickDir', 'out', 'XLim', [t1, t2]);

        subplot(4, 1, 3); hold on;
        plot(zts(ists), emg(ists, 1), '-k', 'LineWidth', 0.1);
        ylabel('emg1');
        set(gca, 'TickDir', 'out', 'XLim', [t1, t2], 'YLim', [-300, 300]);

        subplot(4, 1, 4); hold on;
        plot(zts(ists), emg(ists, 2), '-k', 'LineWidth', 0.1);
        ylabel('emg2');
        set(gca, 'TickDir', 'out', 'XLim', [t1, t2], 'YLim', [-300, 300]);

        bstr = sprintf('%s_bk%.2d', ds(nn).fid, ii);
        sgtitle(bstr, 'Interpreter', 'none');
        print(fullfile(workdir, 'check_misc', bstr), '-djpeg');
        close;
    end
    fprintf('%gs.\n', toc(to));
end

%%
for nn = 1:numel(ds)
    fp = fullfile(workdir, 'prep_sig', [ds(nn).fid, '.mat']);
    to = tic;
    fprintf('[%.3d] %s >>> ', nn, ds(nn).fid);
    load(fp, 'eeg1', 'eeg2', 'emg1', 'emg2', 'fs');

    ts = (0:1:numel(eeg1)-1)' / fs;
    h = design(fdesign.bandpass('N,F3dB1,F3dB2', 4, 3, 100, fs), 'butter');
    sig = filtfilt(h.sosMatrix, h.ScaleValues, [eeg1, eeg2]);

    emg = log(movmean([emg1, emg2] .^ 2, fs * 60, 1) .^ 0.5);
    [nt, nch] = size(emg);
    thr = zeros(2, nch);
    for kk = 1:nch        
        gm = fitgmdist(emg(round(rand(100000, 1) * nt), kk), 2, 'Replicates', 5);
        dx = (min(gm.mu):0.01:max(gm.mu))';
        [~, dpdf] = min(pdf(gm, dx));
        thr(1, kk) = dx(dpdf);
        thr(2, kk) = max(gm.mu) * 2 - dx(dpdf);
    end
    % disp(thr);
    is_bsl = all(emg > thr(1, :) & emg < thr(2, :), 2);
    sigma = mad(sig(is_bsl, :), 1, 1) * 1.4826;
    ds(nn).sigma = sigma;

    [ds(nn).ev_se, ds(nn).zts_se, ds(nn).se_att, ds(nn).zts_spk] = deal(cell(1, 2));
    for kk = 1:2
        [is_hamp, ev_hamp] = eegDetectHighAmplitude(sig(:, kk), ts, [], sigma(kk), 1, [3, 2]);
        [is_sspk, ev_sspk] = eegDetectSharpSpikes(sig(:, kk), ts, [], sigma(kk), 1, 3);
        [is_swd, ev_swd, idx_spk] = eegDetectSpikeWave(sig(:, kk), ts, [], sigma(kk), 1, [5, 3]);
        [is_se, ev_se, se_att] = eegCalibrateSeizure(sig(:, kk) / sigma(kk), ts, 0.5, ...
            ev_hamp, [-inf, inf], ev_sspk, [-inf, inf], ev_swd, [-inf, inf]);

        ds(nn).ev_se{kk} = ev_se;
        ds(nn).zts_se{kk} = [ts(ev_se(:, 1)), ts(ev_se(:, 2))];
        ds(nn).se_att{kk} = se_att;
        ds(nn).zts_spk{kk} = [idx_spk, ts(idx_spk)];
    end
    fprintf('%gs.\n', toc(to));
end

%%
seizure_duration = 10;
seizure_gap = 10;
spike_interval = 1;
fs_bs = 1/5;
tx = (0:1:24)';
for nn = 1:length(ds)
    kk = 2;
    zts = (0.5:1:(60 * 60 * 24 * 3 * fs_bs))' / fs_bs;
    bs = ones(size(zts));
    bs_duration = sum(bs == [1, 0, -1]) / fs_bs / 60 / 60;

    ea_att = ds(nn).se_att{kk};
    if isempty(ea_att); continue; end
    is_se = ea_att.Duration >= seizure_duration;

    ea_zts = ds(nn).zts_se{kk};

    se_zts = ea_zts(is_se, :);
    while size(se_zts, 1) > 1
        se_gap = se_zts(2:end, 1) - se_zts(1:end-1, 2);
        if ~any(se_gap <= seizure_gap)
            break;
        else
            ids = find(se_gap <= seizure_gap);
            t1 = se_zts(:, 1);
            t2 = se_zts(:, 2);
            t1(ids + 1) = [];
            t2(ids) = [];
            se_zts = [t1, t2];
        end
    end
    
    ea_zts(is_se, :) = [];  %

    spk_zts = ds(nn).zts_spk{kk}(:, 2);
    spk_zts = spk_zts(diff([-inf; spk_zts]) > spike_interval); % interval: 1s

    in_se = any(spk_zts - se_zts(:, 1)' >= 0 & spk_zts - se_zts(:, 2)' <= 0, 2);    
    spk_zts(in_se) = [];    %

    spk_ztx = histcounts(mod(spk_zts / 60 / 60, 24), tx);
    ds(nn).spk_zt_day = spk_ztx;

    ea_ztx = histcounts(mod(ea_zts / 60 / 60, 24), tx);
    ds(nn).ea_zt_day = ea_ztx;

    ds(nn).spk_frequency = numel(spk_zts) / sum(bs_duration);

    spk_bs = interp1(zts, bs, spk_zts, 'nearest');
    ds(nn).spk_ph_bs = sum(spk_bs == [1, 0, -1]) ./ bs_duration; 

    ea_bs = interp1(zts, bs, ea_zts(:, 1), 'nearest');
    ds(nn).ea_ph_bs = sum(ea_bs == [1, 0, -1]) ./ bs_duration; 

    ds(nn).ea_frequency = size(ea_zts, 1) / sum(bs_duration); % per hour
    ds(nn).ea_duration = median(ea_att.Duration(~is_se));
    ds(nn).ea_per_hour = sum(ea_att.Duration(~is_se)) / sum(bs_duration);

    ds(nn).se_frequency = size(se_zts, 1) / (sum(bs_duration) / 24); % per day
    ds(nn).se_duration = median(ea_att.Duration(is_se));
    if isnan(ds(nn).se_duration); ds(nn).se_duration = 0; end
    ds(nn).se_per_day = sum(ea_att.Duration(is_se)) / (sum(bs_duration) / 24);
    ds(nn).seizure_ts = se_zts;
end

%%





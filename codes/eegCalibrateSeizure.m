function [is_se, ev_se, se_att] = eegCalibrateSeizure(sig, ts, ev_perc, ...
    ev_hamp, len_hamp, ev_sspk, len_sspk, ev_swd, len_swd)
sz_sig = size(ts);
is_se = false(sz_sig);

if ~isempty(ev_hamp)
    ev_hamp = split_ev_gap(ev_hamp, ts);
    ev_hamp = filter_ev_att(ev_hamp, ts(ev_hamp(:, 2)) - ts(ev_hamp(:, 1)), len_hamp);
end
if ~isempty(ev_sspk)
    ev_sspk = split_ev_gap(ev_sspk, ts);
    ev_sspk = filter_ev_att(ev_sspk, ts(ev_sspk(:, 2)) - ts(ev_sspk(:, 1)), len_sspk);
end
if ~isempty(ev_swd)
    ev_swd = split_ev_gap(ev_swd, ts);
    ev_swd = filter_ev_att(ev_swd, ts(ev_swd(:, 2)) - ts(ev_swd(:, 1)), len_swd);
end

is_sspk = evToVec(ev_sspk, sz_sig);
is_swd = evToVec(ev_swd, sz_sig);

nhamp = size(ev_hamp, 1);
ratios = zeros(nhamp, 1);
for ii = 1:nhamp
    ints = ev_hamp(ii, 1):ev_hamp(ii, 2);
    ratios(ii) = mean(is_sspk(ints) | is_swd(ints));
    if ratios(ii) >= ev_perc
        is_se(ints) = true;
    end
end
ratios = ratios(ratios >= ev_perc);

ibeg = find(diff([false; is_se]) == 1);
iend = find(diff([is_se; false]) == -1);
ev_se = [ibeg, iend];
ev_se = split_ev_gap(ev_se, ts);

nev = size(ev_se, 1);
[se_dur, se_amp, se_rms] = deal(nan(nev, 1));
for kk = 1:nev
    se_dur(kk) = diff(ts(ev_se(kk, :)));
    isig = sig(ev_se(kk, 1):ev_se(kk, 2));
    se_amp(kk) = max(abs(isig));
    se_rms(kk) = rms(isig);
end
se_att = table(ratios, se_dur, se_amp, se_rms, ...
    'VariableNames', {'SWD', 'Duration', 'Amplitude', 'RMS'});
end

function evt = filter_ev_att(evt, att, rg)
isin = att > min(rg) & att < max(rg);
evt = evt(isin, :);
end

function evt = split_ev_gap(evt, ts)
dts = diff(ts);
itl = median(dts);
tgap = itl * 10;
isgap = find(dts >= tgap);
ngap = numel(isgap);
for kk = 1:ngap
    idx = find(evt(:, 1) <= isgap(kk) & evt(:, 2) >= isgap(kk) + 1);
    if ~any(idx); continue; end
    evt = cat(1, evt(1:idx-1, :), [evt(idx, 1), isgap(kk)], ...
        [isgap(kk) + 1, evt(idx, 2)], evt(idx+1:end, :));
end
end

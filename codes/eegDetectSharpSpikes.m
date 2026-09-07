function [is_sspk, ev_sspk] = eegDetectSharpSpikes(sig, ts, freqr, isbsl, sw, xfold)
fs = 1/median(diff(ts));
nt = length(sig);
nsw = round(fs * sw);
nhsw = round(fs * sw / 2);

if ~isempty(freqr)
    h = design(fdesign.bandpass('N,F3dB1,F3dB2', 6, ...
        min(freqr), max(freqr), fs), 'butter');
    sig = filtfilt(h.sosMatrix, h.ScaleValues, sig);
end
if islogical(isbsl)
    sigma = mad(sig(isbsl), 1) * 1.4826;
else
    sigma = isbsl;
end

ishigh = sig > sigma * xfold;
ibeg = find(diff([false; ishigh]) == 1);
iend = find(diff([ishigh; false]) == -1);

idip = [1; find(sig(1:end-1) .* sig(2:end) <= 0) + 1; nt];
[ibeg, iend] = evExtendMerge(ibeg, iend, idip);
dur = ts(iend) - ts(ibeg) + 0.5 / fs;
% ibeg = ibeg(dur < 0.2);
% iend = iend(dur < 0.2);
wid = spikeHalfWidth(sig, ts, ibeg, iend);
ibeg = ibeg(wid < 0.1 & wid > 0.01 & dur < 0.2);
iend = iend(wid < 0.1 & wid > 0.01 & dur < 0.2);

spkv = nan(nt, 1);
for ii = 1:numel(ibeg)
    spkv(ibeg(ii):iend(ii)) = ii;
end
spk_max = movmax(spkv, nsw);
spk_min = movmin(spkv, nsw);
spkv = (spk_max - spk_min + 1) / sw;

ishigh = spkv >= 3;
ibeg = max(find(diff([false; ishigh]) == 1) - nhsw, 1);
iend = min(find(diff([ishigh; false]) == -1) + nhsw, nt);
[ibeg, iend] = evExtendMerge(ibeg, iend);
is_sspk = evToVec([ibeg, iend], [nt, 1]);
ev_sspk = [ibeg, iend];
end
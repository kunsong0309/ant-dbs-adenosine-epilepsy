function [is_swd, ev_swd, ipk, hws] = eegDetectSpikeWave(sig, ts, freqr, isbsl, sw, xfold, hw)
if nargin < 7; hw = [0.01, 0.1]; end
fs = 1/median(diff(ts));
nt = length(sig);
nsw = round(fs * sw);

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
amp = movmean(abs(hilbert(sig)), nsw);

ishigh = -sig > sigma * max(xfold);
ibeg = find(diff([false; ishigh]) == 1);
iend = find(diff([ishigh; false]) == -1);

idip = [1; find(sig(1:end-1) .* sig(2:end) <= 0) + 1; nt];
[ibeg, iend] = evExtendMerge(ibeg, iend, idip);
% dur = ts(iend) - ts(ibeg) + 0.5 / fs;
% ibeg = ibeg(dur < 0.1);
% iend = iend(dur < 0.1);
[wid, ipk] = spikeHalfWidth(-sig, ts, ibeg, iend);
is_val = wid < max(hw) & wid > min(hw);
ibeg = ibeg(is_val);
iend = iend(is_val);
ipk = ipk(is_val);
hws = wid(is_val);
% ibeg = ibeg(wid < 0.075);
% iend = iend(wid < 0.075);

dsig = diff(amp);
idip = [1; find(dsig(1:end-1) < 0 & dsig(2:end) >= 0 & ...
    amp(2:end-1) < sigma * min(xfold)) + 1; nt];
[ibeg, iend] = evExtendMerge(ibeg, iend, idip);
is_swd = evToVec([ibeg, iend], [nt, 1]);
ev_swd = [ibeg, iend];
end
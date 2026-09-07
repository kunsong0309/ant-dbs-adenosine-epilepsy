function [is_hamp, ev_hamp] = eegDetectHighAmplitude(sig, ts, freqr, isbsl, sw, xfold)
fs = 1/median(diff(ts));
nt = length(ts);
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

ishigh = amp > sigma * max(xfold);
ibeg = find(diff([false; ishigh]) == 1);
iend = find(diff([ishigh; false]) == -1);

dsig = amp - sigma * min(xfold);
idip = [1; find(dsig(1:end-1) .* dsig(2:end) <= 0) + 1; nt];
for ii = 1:numel(ibeg)
    imin = idip - ibeg(ii) <= 0;
    ibeg(ii) = idip(find(imin, 1, 'last'));
    imax = idip - iend(ii) >= 0;
    iend(ii) = idip(find(imax, 1, 'first'));
end
iiv = find(ibeg(2:end) - iend(1:end-1) <= 1);
ibeg(iiv + 1) = [];
iend(iiv) = [];
is_hamp = evToVec([ibeg, iend], [nt, 1]);
ev_hamp = [ibeg, iend];
end
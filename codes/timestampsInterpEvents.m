function [omat, omat2] = timestampsInterpEvents(vmat, ts, evt, outtype, method, evdim)
[nev, nts] = size(evt);
evt = evt';
ndim = ndims(vmat);

if nargin < 4 || isempty(outtype)
    if nts == 2
        outtype = 'avgrg';
    else
        outtype = 'trace';
    end
end

if nargin < 5 || isempty(method)
    method = 'linear';
end

if nargin < 6 || isempty(evdim)
    evdim = ndim + 1;
end

[omat, omat2] = deal(cell(nev, 1));
f = griddedInterpolant(ts, vmat, method, 'none');
switch outtype
    case 'avgrg'
        for kk = 1:nev
            xt = ts(ts > evt(1, kk) & ts < evt(2, kk));
            omat{kk} = mean(f(xt), 1, 'omitnan');
        end
    case 'trace'
        for kk = 1:nev
            xt = evt(:, kk);
            omat{kk} = f(xt);
        end
    case 'max'       
        for kk = 1:nev
            xt = ts(ts > evt(1, kk) & ts < evt(2, kk));
            [omat{kk}, omat2{kk}] = max(movmean(f(xt), 10, 1, 'omitnan'), [], 1);
        end
    case 'median'
        for kk = 1:nev
            xt = ts(ts > evt(1, kk) & ts < evt(2, kk));
            omat{kk} = median(f(xt), 1, 'omitnan');
        end
end

omat = cat(evdim, omat{:});
omat2 = cat(evdim, omat2{:});
end
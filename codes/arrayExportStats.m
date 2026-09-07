function s = arrayExportStats(vmat, gnames1, gnames2, vname, fpath, type)
if nargin < 6; type = ''; end
s = [];

ngrp1 = length(gnames1);
ngrp2 = length(gnames2);

vgrp = cellfun(@(x, y)(repmat({x}, size(y, 1), 1)), gnames2, vmat, 'UniformOutput', false);
vgrp = table(cat(1, vgrp{:}), 'VariableNames', {'Group'});

if ~isempty(type)
    switch type
        case 'ttest2'
            do_stats_ttest2(vmat, vgrp, vname, fpath);
        case 'ttest'
            do_stats_paired_ttest(vmat, gnames1, vgrp, fpath);
        case 'anova'
            s = do_stats_anova(vmat, vgrp, vname, fpath);
        case 'ranova'
            s = do_stats_ranova(vmat, gnames1, vgrp, vname, fpath);
        case 'ranova2'
            s = do_stats_ranova2(vmat, gnames1, gnames2, vgrp, vname, fpath);
        case 'manova'
            s = do_stats_manova(vmat, gnames1, vgrp, vname, fpath);
        case 'ranksum'
            do_stats_ranksum(vmat, vgrp, vname, fpath);
        case 'signrank'
            do_stats_signrank(vmat, gnames1, vgrp, fpath);
        case 'kwtest'
            do_stats_kwtest(vmat, gnames2, vgrp, vname, fpath);
        case 'friedman'
            do_stats_friedman(vmat, gnames1, vgrp, fpath);
        case 'mult-ranksum'
            do_stats_mult_ranksum(vmat, gnames1, vgrp, vname, fpath);
        otherwise
            disp('TODO.');
    end
else
    isnorm = false(ngrp2, ngrp1);
    for ii = 1:ngrp2
        for jj = 1:ngrp1
            isnorm(ii, jj) = lillietest(vmat{ii}(:, jj)) == 0;
        end
    end
    % disp(isnorm);
    isnorm = all(isnorm, 'all');

    if isnorm
        if ngrp1 == 1 && ngrp2 == 2
            do_stats_ttest2(vmat, vgrp, vname, fpath);
        elseif ngrp1 <= 2 && ngrp2 == 1
            do_stats_paired_ttest(vmat, gnames1, vgrp, fpath);
        elseif ngrp1 == 1 && ngrp2 > 2
            s = do_stats_anova(vmat, vgrp, vname, fpath);
        elseif ngrp1 > 2 && ngrp2 == 1
            s = do_stats_ranova(vmat, gnames1, vgrp, vname, fpath);
        elseif ngrp1 >= 2 && ngrp2 >= 2
            s = do_stats_manova(vmat, gnames1, vgrp, vname, fpath);
        else
            disp('TODO.');
        end
    else
        if ngrp1 == 1 && ngrp2 == 2
            do_stats_ranksum(vmat, vgrp, vname, fpath);
        elseif ngrp1 <= 2 && ngrp2 == 1
            do_stats_signrank(vmat, gnames1, vgrp, fpath);
        elseif ngrp1 == 1 && ngrp2 > 2
            do_stats_kwtest(vmat, gnames2, vgrp, vname, fpath);
        elseif ngrp1 > 2 && ngrp2 == 1 % ????
            do_stats_friedman(vmat, gnames1, vgrp, fpath);
        elseif ngrp1 >= 2 && ngrp2 == 2
            % do_stats_mult_ranksum(vmat, gnames1, vgrp, vname, fpath);
            s = do_stats_manova(vmat, gnames1, vgrp, vname, fpath);
        elseif ngrp1 >= 2 && ngrp2 > 2 % ????
            s = do_stats_manova(vmat, gnames1, vgrp, vname, fpath);
        else
            disp('TODO.');
        end
    end
end
end

function do_stats_ranksum(vmat, vgrp, vname, fpath)
T = cat(2, vgrp, array2table(cat(1, vmat{:}), 'VariableNames', {vname}));
writetable(T, fpath, 'Sheet', 'Data', 'WriteMode', 'replacefile');

[P, ~, S] = ranksum(vmat{1}, vmat{2});
S = cat(2, struct2table(S), table(P, 'VariableNames', {'pvalue'}));
writetable(S, fpath, 'WriteRowNames', true, 'Sheet', 'Stats', 'Range', 'A1');
end

function do_stats_mult_ranksum(vmat, gnames1, vgrp, vname, fpath)
T = cat(2, vgrp, array2table(cat(1, vmat{:}), 'VariableNames', gnames1));
writetable(T, fpath, 'Sheet', 'Data', 'WriteMode', 'replacefile');

ngrp1 = numel(gnames1);
[P, S] = deal(cell(ngrp1, 1));
for kk = ngrp1:-1:1
    [P{kk}, ~, S{kk}] = ranksum(vmat{1}(:, kk), vmat{2}(:, kk));
end
S = cat(1, S{:});
P = cat(1, P{:});
Padj = mafdr(P, 'BHFDR', true);
S = cat(2, cell2table(gnames1(:), 'VariableNames', {vname}), struct2table(S), ...
    array2table([P, Padj], 'VariableNames', {'pvalue', 'padjust'}));
writetable(S, fpath, 'WriteRowNames', true, 'Sheet', 'Stats', 'Range', 'A1');
end

function do_stats_signrank(vmat, gnames1, vgrp, fpath)
T = cat(2, vgrp, array2table(cat(1, vmat{:}), 'VariableNames', gnames1));
writetable(T, fpath, 'Sheet', 'Data', 'WriteMode', 'replacefile');

if size(vmat{1}, 2) == 1
    [P, ~, S] = signrank(vmat{1});
else
    [P, ~, S] = signrank(vmat{1}(:, 1), vmat{1}(:, 2));
end
S = cat(2, struct2table(S), table(P, 'VariableNames', {'pvalue'}));
writetable(S, fpath, 'WriteRowNames', true, 'Sheet', 'Stats', 'Range', 'A1');
end

function do_stats_ttest2(vmat, vgrp, vname, fpath)
T = cat(2, vgrp, array2table(cat(1, vmat{:}), 'VariableNames', {vname}));
writetable(T, fpath, 'Sheet', 'Data', 'WriteMode', 'replacefile');

[~, P, ~, S] = ttest2(vmat{1}, vmat{2});
S = cat(2, struct2table(S), table(P, 'VariableNames', {'pvalue'}));
writetable(S, fpath, 'WriteRowNames', true, 'Sheet', 'Stats', 'Range', 'A1');
end

function do_stats_paired_ttest(vmat, gnames1, vgrp, fpath)
T = cat(2, vgrp, array2table(cat(1, vmat{:}), 'VariableNames', gnames1));
writetable(T, fpath, 'Sheet', 'Data', 'WriteMode', 'replacefile');

if size(vmat{1}, 2) == 1
    [~, P, ~, S] = ttest(vmat{1});
else
    [~, P, ~, S] = ttest(vmat{1}(:, 1), vmat{1}(:, 2));
end
S = cat(2, struct2table(S), table(P, 'VariableNames', {'pvalue'}));
writetable(S, fpath, 'WriteRowNames', true, 'Sheet', 'Stats', 'Range', 'A1');
end

function s = do_stats_ranova(vmat, gnames1, vgrp, vname, fpath)
ngrp1 = numel(gnames1);
vtime = arrayfun(@(x)(sprintf('y%d', x)), 1:ngrp1, 'UniformOutput', false);
s = struct;

between = cat(2, vgrp, array2table(cat(1, vmat{:}), 'VariableNames', vtime));
within = table(vtime(:), 'VariableNames', {vname});
rm = fitrm(between, [strjoin(vtime, ','), ' ~ 1'], ...
    'WithinDesign', within, 'WithinModel', vname);

T = cat(2, vgrp, array2table(cat(1, vmat{:}), 'VariableNames', gnames1));
s.values = T;
writetable(T, fpath, 'Sheet', 'Data', 'WriteMode', 'replacefile');
s.main = ranova(rm);
writetable(ranova(rm), fpath, 'WriteRowNames', true, ...
    'Sheet', 'Stats', 'Range', 'A1');

tbl = multcompare(rm, vname);
tbl{:, [1, 2]} = replace(tbl{:, [1, 2]}, vtime, gnames1);
s.repeat = tbl;
writetable(tbl, fpath, 'WriteRowNames', true, 'Sheet', 'Stats', 'Range', 'A10');
end

function s = do_stats_ranova2(vmat, gnames1, gnames2, vgrp, vname, fpath)
ngrp1 = numel(gnames1);
ngrp2 = numel(gnames2);
ngrp_ = ngrp1 * ngrp2;
vtime = arrayfun(@(x)(sprintf('y%d', x)), 1:ngrp_, 'UniformOutput', false);
s = struct;

between = array2table(cat(2, vmat{:}), 'VariableNames', vtime);
within = combinations(gnames2, gnames1);
within.Properties.VariableNames = vname([2, 1]);
rm = fitrm(between, [strjoin(vtime, ','), ' ~ 1'], ...
    'WithinDesign', within, 'WithinModel', ...
    [vname{1} '+', vname{2}, '+', vname{1}, '*', vname{2}]);

T = cat(2, vgrp, array2table(cat(1, vmat{:}), 'VariableNames', gnames1));
s.values = T;
writetable(T, fpath, 'Sheet', 'Data', 'WriteMode', 'replacefile');
s.main = manova(rm);
writetable(manova(rm), fpath, 'WriteRowNames', true, ...
    'Sheet', 'Stats', 'Range', 'A1');

tbl = multcompare(rm, vname{2}, 'By', vname{1});
tbl = tidy_up_stats(tbl, {gnames1, gnames2, gnames2});
s.groups = tbl;
writetable(tbl, fpath, 'WriteRowNames', true, 'Sheet', 'Stats', 'Range', 'A20');

tbl = multcompare(rm, vname{1}, 'By', vname{2});
tbl = tidy_up_stats(tbl, {gnames2, gnames1, gnames1});
s.repeat = tbl;
writetable(tbl, fpath, 'WriteRowNames', true, 'Sheet', 'Stats', 'Range', ...
    sprintf('A%d', 23 + size(s.groups, 1)));
end

function s = do_stats_manova(vmat, gnames1, vgrp, vname, fpath)
ngrp1 = numel(gnames1);
vtime = arrayfun(@(x)(sprintf('y%.2d', x)), 1:ngrp1, 'UniformOutput', false);
s = struct;

between = cat(2, vgrp, array2table(cat(1, vmat{:}), 'VariableNames', vtime));
within = table(vtime(:), 'VariableNames', {vname});
rm = fitrm(between, [strjoin(vtime, ','), ' ~ Group'], ...
    'WithinDesign', within, 'WithinModel', vname);

T = cat(2, vgrp, array2table(cat(1, vmat{:}), 'VariableNames', gnames1));
s.values = T;
writetable(T, fpath, 'Sheet', 'Data', 'WriteMode', 'replacefile');
s.main = ranova(rm);
writetable(ranova(rm), fpath, 'WriteRowNames', true, ...
    'Sheet', 'Stats', 'Range', 'A1');
tbl = multcompare(rm, 'Group', 'By', vname);
for kk = 1:ngrp1
    ri = cellfun(@(x)(strcmp(x, vtime{kk})), tbl.(vname));
    tbl.(vname)(ri) = gnames1(kk);
end
s.groups = tbl;
writetable(tbl, fpath, 'WriteRowNames', true, 'Sheet', 'Stats', 'Range', 'A10');
nr = size(tbl, 1);

tbl = multcompare(rm, vname, 'By', 'Group');
for kk = 1:ngrp1
    ri = cellfun(@(x)(strcmp(x, vtime{kk})), tbl.([vname, '_1']));
    tbl.([vname, '_1'])(ri) = gnames1(kk);
    ri = cellfun(@(x)(strcmp(x, vtime{kk})), tbl.([vname, '_2']));
    tbl.([vname, '_2'])(ri) = gnames1(kk);
end
s.repeat = tbl;
writetable(tbl, fpath, 'WriteRowNames', true, 'Sheet', 'Stats', ...
    'Range', sprintf('A%d', nr+15));
end

function s = do_stats_anova(vmat, vgrp, vname, fpath)
s = struct;

T = cat(2, vgrp, array2table(cat(1, vmat{:}), 'VariableNames', {vname}));
s.values = T;
writetable(T, fpath, 'Sheet', 'Data', 'WriteMode', 'replacefile');

aov = anova(T, vname);
s.main = stats(aov);
writetable(stats(aov), fpath, 'WriteRowNames', true, 'Sheet', 'Stats', 'Range', 'A1');
s.groups = multcompare(aov);
writetable(multcompare(aov), fpath, 'WriteRowNames', true, 'Sheet', 'Stats', 'Range', 'A10');
end

function do_stats_kwtest(vmat, gnames2, vgrp, vname, fpath)
T = cat(2, vgrp, array2table(cat(1, vmat{:}), 'VariableNames', {vname}));
writetable(T, fpath, 'Sheet', 'Data', 'WriteMode', 'replacefile');

[~, tbl, stats] = kruskalwallis(cat(1, vmat{:}), vgrp.Group, 'off');
writetable(cell2table(tbl(2:end, :), 'VariableNames', tbl(1, :)), ...
    fpath, 'WriteRowNames', true, 'Sheet', 'Stats', 'Range', 'A1');
tbl = array2table(multcompare(stats, 'Display', 'off'), 'VariableNames', ...
    ["Group A", "Group B", "Lower Limit", "A-B", "Upper Limit", "P-value"]);
tbl.("Group A") = gnames2(tbl.("Group A"))';
tbl.("Group B") = gnames2(tbl.("Group B"))';
writetable(tbl, fpath, 'WriteRowNames', true, 'Sheet', 'Stats', 'Range', 'A10');
end

function do_stats_friedman(vmat, gnames1, vgrp, fpath)
ngrp1 = numel(gnames1);

T = cat(2, vgrp, array2table(cat(1, vmat{:}), 'VariableNames', gnames1));
writetable(T, fpath, 'Sheet', 'Data', 'WriteMode', 'replacefile');

[P, ~, S] = friedman(vmat{1}, 1, "off");
tbl = cat(2, struct2table(S), table(P, 'VariableNames', {'pvalue'}));
writetable(tbl, fpath, 'WriteRowNames', true, 'Sheet', 'Stats', 'Range', 'A1');

tbl = array2table(multcompare(S, 'Display', 'off'), 'VariableNames', ...
    {'Group_1', 'Group_2', 'Lower', 'Upper', 'Difference', 'pValue'});
tbl = convertvars(tbl, {'Group_1', 'Group_2'}, 'cell');
for kk = 1:ngrp1
    ri = cellfun(@(x)(strcmp(string(x), num2str(kk))), tbl.Group_1);
    tbl.Group_1(ri) = gnames1(kk);
    ri = cellfun(@(x)(strcmp(string(x), num2str(kk))), tbl.Group_2);
    tbl.Group_2(ri) = gnames1(kk);
end
writetable(tbl, fpath, 'WriteRowNames', true, 'Sheet', 'Stats', 'Range', 'A10');
end

function tbl = tidy_up_stats(tbl, vargin)
nv = length(vargin);
nr = size(tbl, 1);
oo = zeros(nr, nv);
for vv = 1:nv
    [~, oo(:, vv)] = ismember(tbl{:, vv}, vargin{vv});
end
[~, oo] = sortrows(oo);
tbl = tbl(oo, :);
end
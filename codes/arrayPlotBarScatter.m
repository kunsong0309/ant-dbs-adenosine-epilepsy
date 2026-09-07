function arrayPlotBarScatter(vmat, gnames1, gnames2, gclrs, ylim, vname, ttl)
ngrp1 = length(gnames1);
ngrp2 = length(gnames2);

vn = cellfun(@(x)(sum(~isnan(x), 1)), vmat, 'UniformOutput', false);
vn = cat(1, vn{:});
vmean = cellfun(@(x)(mean(x, 1, 'omitnan')), vmat, 'UniformOutput', false);
vmean = cat(1, vmean{:});
vsem = cellfun(@(x)(std(x, 0, 1, 'omitnan')), vmat, 'UniformOutput', false);
vsem = cat(1, vsem{:}) ./ sqrt(vn);

rng(42);
hold on;
hbars = bar((1:ngrp1)', vmean');
for gg = 1:ngrp2
    hbars(gg).FaceColor = gclrs(gg, :);
    hbars(gg).EdgeColor = gclrs(gg, :);
    hbars(gg).LineWidth = 1;
    plot(hbars(gg).XEndPoints + normrnd(0.05, 0.00, size(vmat{gg})), vmat{gg}, 'o', ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', ...
        'LineWidth', 1, 'MarkerSize', 1);
    errorbar(hbars(gg).XEndPoints - 0.05, vmean(gg, :), vsem(gg, :), '.', ...
        'Color', 'k', 'LineWidth', 1, 'CapSize', 2);
end
set(gca, 'TickDir', 'out', 'TickLength', [0.005, 0.1], 'FontSize', 10, ...
    'XTick', (1:ngrp1), 'XTickLabel', gnames1);
if ~isempty(ylim); set(gca, 'YLim', ylim); end
ylabel(vname, 'FontSize', 10);
legend(gnames2, 'Box', 'off', 'Location', 'best', ...
    'Interpreter', 'none', 'FontSize', 10);
title(ttl, 'FontSize', 10);

hold off;
end
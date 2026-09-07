%% for amplitude
da0 = readtable('dbs_dependency_current_ado_peak_ctrl.xlsx', ...
    'Sheet', 'Data', 'VariableNamingRule', 'preserve');
da0 = da0.("200 uA");
da1 = readtable('dbs_dependency_current_ado_peak_les.xlsx', ...
    'VariableNamingRule', 'preserve');
da1 = da1.("200 uA");
rng(42);

figure; set(gcf, 'Position', [50, 50, 780, 780]);
subplot('Position', [0.1, 0.5, 0.14, 0.25]); hold on;
bar(1, mean(da0), 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.3, ...
    'BarWidth', 0.67);
errorbar(1, mean(da0), std(da0)/sqrt(length(da0)), 'Color', 'k', ...
    'LineWidth', 1.3, 'CapSize', 10);
plot(randn(size(da0))*0.1 + 1, da0, '.', 'MarkerSize', 10, 'Color', 'r');
bar(2, mean(da1), 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.3, ...
    'BarWidth', 0.67);
errorbar(2, mean(da1), std(da1)/sqrt(length(da1)), 'Color', 'k', ...
    'LineWidth', 1.3, 'CapSize', 10);
plot(randn(size(da1))*0.1 + 2, da1, '.', 'MarkerSize', 10, 'Color', 'b');

[~, p] = ttest2(da0, da1);
plot([1, 2], [1, 1] * 1.05 * max([da0; da1]), '-k', 'LineWidth', 1.3);
text(1.5, 1.07 * max([da0; da1]), sprintf('\\itP\\rm = %.3f', p), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

set(gca, 'TickDir', 'out', 'XTick', [1, 2], 'XTickLabel', {'Control', 'Leision'}, ...
    'XLim', [0.3, 2.7], 'FontSize', 10, 'LineWidth', 1.3, 'YLim', [0, 65]);
ylabel(sprintf('Peak \\DeltaF/F (%%)'), 'FontSize', 10);
title('@ 10 s, 200 \muA', 'FontSize', 10);

%% for duration
da0 = readtable('dbs_dependency_duration_ado_auc_ctrl.xlsx', ...
    'Sheet', 'Data', 'VariableNamingRule', 'preserve');
da0 = da0.Value(cellfun(@(x)(strcmp(x, '300 s')), da0.Group));
da1 = readtable('dbs_dependency_duration_ado_auc_les.xlsx', ...
    'VariableNamingRule', 'preserve');
da1 = da1.("300 sec");
rng(42);

% figure; set(gcf, 'Position', [50, 50, 780, 780]);
subplot('Position', [0.35, 0.5, 0.14, 0.25]); hold on;
bar(1, mean(da0), 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.3, ...
    'BarWidth', 0.67);
errorbar(1, mean(da0), std(da0)/sqrt(length(da0)), 'Color', 'k', ...
    'LineWidth', 1.3, 'CapSize', 10);
plot(randn(size(da0))*0.1 + 1, da0, '.', 'MarkerSize', 10, 'Color', 'r');
bar(2, mean(da1), 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.3, ...
    'BarWidth', 0.67);
errorbar(2, mean(da1), std(da1)/sqrt(length(da1)), 'Color', 'k', ...
    'LineWidth', 1.3, 'CapSize', 10);
plot(randn(size(da1))*0.1 + 2, da1, '.', 'MarkerSize', 10, 'Color', 'b');

[~, p] = ttest2(da0, da1);
plot([1, 2], [1, 1] * 1.05 * max([da0; da1]), '-k', 'LineWidth', 1.3);
text(1.5, 1.08 * max([da0; da1]), sprintf('\\itP\\rm = %.3f', p), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

set(gca, 'TickDir', 'out', 'XTick', [1, 2], 'XTickLabel', {'Control', 'Leision'}, ...
    'XLim', [0.3, 2.7], 'FontSize', 10, 'LineWidth', 1.3, ...
    'YTickLabel', get(gca, 'YTick')/100);
ylabel(sprintf('\\Sigma_{Floure.} (\\DeltaF/F) (x10^{2} %%)'), 'FontSize', 10);

title('@ 150 \muA, 300 s', 'FontSize', 10);

%%
print(gcf, 'dbs_dependence_stats_ctrl', '-djpeg');
print(gcf, 'dbs_dependence_stats_ctrl', '-dpdf', '-r300', '-vector');













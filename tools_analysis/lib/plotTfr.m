function plotTfr(ftDataStruct, cfg, bfs, savename)
figure('Position', [100 100 1600 1600]);

cfg.figure = subplot(1,2,1);
ft_singleplotTFR(cfg, ftDataStruct);
title('TFR')
pbaspect([1 1 1]);  % Makes this subplot square

subplot(1,2,2);
load('bayes_colourmap.mat');
exponential_minmax = 10;
val_col_map = logspace(-exponential_minmax, exponential_minmax, size(colours, 1));
imagesc(log10(bfs)); % plotting the tfr
xticks(1:numel(ftDataStruct.time))
xticklabels(ftDataStruct.time)
xtickangle(45);
yticks(1:numel(ftDataStruct.freq))
yticklabels(floor(ftDataStruct.freq(end:-1:1)))
pbaspect([1 1 1]);  % Makes this subplot square
% f = gcf; f.Position = [10 10 1600 1600];
colormap(subplot(1,2,2),colours);
% adjusting the tick labels for the colourbar
cbh = colorbar;
cbh.Ticks = [-exponential_minmax, -6, -3, -1, -0.5, 0, 0.5, 1, 3, 6, exponential_minmax];
cbh.TickLabels = arrayfun(@(x) ['10^{' num2str(x) '}'], cbh.Ticks, 'UniformOutput', false);
cbh.TickLabels(strcmp(cbh.TickLabels, '10^{0}')) = {'Inconclusive'};
cbh.TickLabels(strcmp(cbh.TickLabels, '10^{0.5}') | strcmp(cbh.TickLabels, '10^{-0.5}')) = {'Moderate'};
cbh.TickLabels(strcmp(cbh.TickLabels, '10^{1}') | strcmp(cbh.TickLabels, '10^{-1}')) = {'Strong'};
title('Bayes Factors')

if exist('savename', 'var') && ~isempty(savename)
    print(savename, '-dpng');
end

return
end
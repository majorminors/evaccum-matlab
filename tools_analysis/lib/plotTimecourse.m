function plotTimecourse(layout,colours,averageDataStructs,dataDiff,subjectDataStructs,xlims,ylims,channels,sensorType,leg,legendloc,theTitle,savename,h)

if size(colours,1) ~= numel(averageDataStructs)
    error('not enough colours')
end
if numel(averageDataStructs) ~= numel(subjectDataStructs)
    error('mismatch between grand average and subjectwise data structures')
end

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
if ~isempty(xlims); cfg.xlim = xlims; end
if ~isempty(ylims); cfg.ylim = ylims; end
cfg.channel    = channels;
% cfg.channel    = 'EEG';
cfg.linecolor = colours;

for i = 1:numel(subjectDataStructs)
    theseErrs{i} = getErrorFromStructs(subjectDataStructs{i});
end; clear i

if ~isempty(dataDiff)
    subplot(2, 1, 1);
    set(gca, 'Position', [0.13 0.57 0.775 0.37]); % set position of the first subplot to make it taller
end

% ft_singleplotER(cfg, averageDataStructs{:});
for i = 1:numel(averageDataStructs)
    if size(averageDataStructs{i}.avg(find(ismember(averageDataStructs{i}.label,channels)),:),1) > 1 % if there is more than one channel
        thisData = mean(averageDataStructs{i}.avg(find(ismember(averageDataStructs{i}.label,channels)),:)); % average across them
    else
        thisData = averageDataStructs{i}.avg(find(ismember(averageDataStructs{i}.label,channels)),:);
    end
    plot(averageDataStructs{i}.time, thisData,...
        'Color',colours(i,:))
    hold on
end; clear i

hold on
plotErrorFromStructs(averageDataStructs,colours,channels,theseErrs)
line([0 0], get(gca,'YLim'), 'Color', 'k', 'LineStyle', '--')
if exist('h','var') % plot our anova results of the RSA
    ax = gca;
    yLow = ax.YLim(1);
    xLeft = ax.XLim(1);
    xRight = ax.XLim(2);
    ind = find(h.Interaction == 1);
    yVals = yLow * ones(size(ind));
    scatter(ind/1000+xLeft, yVals, [] , [1, 0.8196, 0.8627], 'filled');
    ind = find(h.Coherence == 1);
    yVals = (yLow + 5.0000e-07) * ones(size(ind));
    scatter(ind/1000+xLeft, yVals, [], [0.851, 0.8196, 0.93135], 'filled');
    ind = find(h.Rule == 1);
    yVals = (yLow + 1.0000e-06) * ones(size(ind));
    scatter(ind/1000+xLeft, yVals, [], [0.7020, 0.8196, 1], 'filled');
    xlim([xLeft xRight]); clear ax yLow xLeft xRight
end
if ~isempty(xlims); xlim(xlims); end
hold off;
ylabel(['Mean ' sensorType{1} ' Amplitude (' sensorType{2} ')'])
xlabel('Time (s)')
legend(leg, 'Location', legendloc);
title(theTitle, 'FontSize', 14)

if ~isempty(dataDiff)
    subplot(2, 1, 2);
    bf = dataDiff.bfs;
    load('bayes_colourmap.mat'); % in BFF repo
    exponential_minmax=6;
    val_col_map = logspace(-exponential_minmax,exponential_minmax,size(colours,1));
    scatter_colours = zeros(length(dataDiff.time), 3);  % preallocate for efficiency
    for t = 1:length(dataDiff.time)
        [~,idx] = min(abs(val_col_map-bf(t)));
        scatter_colours(t, :) = colours(idx,1:3);
    end
    scatter(dataDiff.time, bf, 30, scatter_colours, 'filled');
    line(get(gca,'XLim'),[1 1], 'Color', [0.7 0.7 0.7], 'LineStyle', '--')
    ax = gca;
    set(ax,'YScale','log','XLim',[dataDiff.time(1),dataDiff.time(end)], ...
            'YLim',[1e-6 1e6],'YTick',10.^(-6:2:6))
    if ~isempty(xlims); xlim(xlims); end
    xlabel('Time (s)')
    ylabel('BF (log scale)')
    colormap(colours)
    cbh = colorbar;
    caxis([-exponential_minmax,exponential_minmax])
    cbh.Units = 'normalized';
    cbh.Limits = [-exponential_minmax,exponential_minmax];
    cbh.Position(1) = 0.92;cbh.Position(3) = 0.01;cbh.Position(4) = ax.Position(4);cbh.Position(2) = ax.Position(2);
    cbh.Label.String = 'Bayes Factor';
    f = gcf; f.Position = [10 10 1600 1600];
    cbh.Ticks = [-6, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 6];
    cbh.TickLabels=arrayfun(@(x) ['10^{' num2str(x) '}'], cbh.Ticks, 'UniformOutput', false);
    cbh.TickLabels(strcmp(cbh.TickLabels,'10^{0}')) = {'Inconclusive'};
    cbh.TickLabels(strcmp(cbh.TickLabels,'10^{0.5}') | strcmp(cbh.TickLabels,'10^{-0.5}')) = {'Moderate'};
    cbh.TickLabels(strcmp(cbh.TickLabels,'10^{1}') | strcmp(cbh.TickLabels,'10^{-1}')) = {'Strong'};

end

if ~isempty(savename)
    print(savename, '-dpng');
end


return
end
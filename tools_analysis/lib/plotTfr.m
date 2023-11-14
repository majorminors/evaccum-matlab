function plotTfr(ftDataStruct, cfg, bfs, savename)
figure('Position', [100 100 1600 1600]);

cfg.figure = subplot(1,2,1);
ft_singleplotTFR(cfg, ftDataStruct);
title('TFR', 'FontSize', 12)
% pbaspect([1 1 1]);  % makes this subplot square

subplot(1,2,2);
load('bayes_colourmap.mat');
h=heatmap(ftDataStruct.time,floor(ftDataStruct.freq(end:-1:1)),bfs,'ColorScaling','log');
caxis([-5 5]) % this does not map properly to the colour values but [-5 5] gets +/-10^2
h.Colormap=colours;
h.GridVisible = 'off';
xlabel('time (ms)');
ylabel('frequency');
set(struct(h).Axes.Title,'String','Bayes Factors')
% % temporarily change axis units 
% originalUnits = h.Units;  % save original units (probaly normalized)
% h.Units = 'centimeters';  % any unit that will result in squares
% % make axes square
% h.Position(3:4) = min(h.Position(3:4))*[1,1]; 
% h.Units = originalUnits; 

% title('Bayes Factors', 'FontSize', 12)

if exist('savename', 'var') && ~isempty(savename)
    print(savename, '-dpng');
end

return
end
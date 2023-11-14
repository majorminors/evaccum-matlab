function doTopos(channels,neighbours,dataAll1,dataAll2,zlims,timestep,timebounds,subplots,totalplots,dataForColour,layout,figsavename,statsavename,datadir,figdir)
if subplots(1)*subplots(2) ~= totalplots; error('right now we need a square of subplots'); end
% first set up the cluster test settings
cfg = [];
cfg.channel = channels;
cfg.latency = 'all';
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT'; % dependent samples t test (also independent: indepsamplesT)
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05; % alpha level of the sample-specific test statistic that will be used for thresholding
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2; % min number of neighbours required to be considered for clustering
cfg.neighbours = neighbours;
cfg.tail = 0; % one or two sided test? (one would be -1 or 1)
cfg.clustertail = 0;
cfg.alpha = 0.025; % alpha for the permutation
cfg.numrandomization = 500; % how many draws?

% now let's set up the design
numSubjs  = numel(dataAll1);
design = zeros(2, numSubjs*2); % design is 2 x num subjects matrix
design(1,:) = [1:numSubjs 1:numSubjs]; % first row is which subject the data belongs to in the two datasets we are comparing
design(2,:) = [ones(1,numSubjs) ones(1,numSubjs)*2]; % second row is which dataset the data belongs to
% so we're contrasting two things, and so we repeat subjects twice in the first row, then do
% a corresponding set of ones and twos for the things
cfg.design = design;
cfg.uvar = 1; % uvar is the row that the subject definition is on
cfg.ivar = 2; % ivar is the row that the dataset (independent variables) definition is on

% get the stats
disp('getting stats')
statSaveFile = fullfile(datadir,[statsavename '.mat']);
if ~exist(statSaveFile,'file')
    disp('no stats saved, calculating and saving')
    stats = ft_timelockstatistics(cfg, dataAll1{:}, dataAll2{:});
    save(statSaveFile,'stats')
else
    disp('stats found, loading')
    load(statSaveFile,'stats')
end; clear statSaveFile

% define parameters for plotting
% timestep      = 0.05; %(in seconds)
sampling_rate = 250;
sample_count  = length(stats.time);
j = [timebounds(1):timestep:timebounds(2)];   % in secs, get some timings for each subplot for the xlims
m = [1:timestep*sampling_rate:sample_count];  % get the timings according to the M/EEG samples to index into the clusters
% get relevant values
pos_cluster_pvals = [stats.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.025);
pos       = ismember(stats.posclusterslabelmat, pos_clust);
neg_cluster_pvals = [stats.negclusters(:).prob];
neg_clust = find(neg_cluster_pvals < 0.025);
neg       = ismember(stats.negclusterslabelmat, neg_clust);% ensure the channels to have the same order in the average and in the statistical output
% which might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(dataForColour.label, stats.label);
% plot
for k = 1:totalplots
   cfg.figure     = subplot(subplots(1),subplots(2),k);
   cfg.xlim       = [j(k) j(k+1)]; % set the xlimits to timepoints
   cfg.zlim = zlims;
   % now we create an index to the channels that are in clusters to be plotted,
   % by timepoint
   % so ff a channel is in a to-be-plotted cluster, then
   % the element of pos/neg_int with an index equal to that channel
   % number will be set to 1 (otherwise 0).
   % so let's check which channels are in the clusters over the
   % time interval of interest
   pos_int        = zeros(numel(dataForColour.label),1);
   pos_int(i1)    = all(pos(i2, m(k):m(k+1)), 2);
   neg_int        = zeros(numel(dataForColour.label),1);
   neg_int(i1)    = all(neg(i2, m(k):m(k+1)), 2);
   cfg.highlight  = 'on';
   % now we pull the indices we created earlier for the channels to
   % highlight
   cfg.highlightchannel = find(pos_int | neg_int);
   cfg.comment    = 'xlim';
   cfg.commentpos = 'title';
   cfg.layout     = layout;
   cfg.interactive = 'no';
   cfg.figure     = 'gca';
   ft_topoplotER(cfg, dataForColour);
end

cbh = colorbar('Location', 'eastoutside', 'Position', [0.93 0.1 0.01 0.8]);
caxis([zlims(1),zlims(2)])
cbh.Units = 'normalized';
cbh.Limits = [zlims(1),zlims(2)];
% cbh.Label.String = 'Bayes Factor';
f = gcf; f.Position = [10 10 1600 1600];
disp('saving figure')
print([figdir filesep figsavename '.png'], '-dpng');
% close all
disp('done!')

end
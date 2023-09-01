clear all

% set up paths
addpath /hpc-software/matlab/cbu/
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/');

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';
datadir = fullfile(rootdir,'data','meg_pilot_4'); addpath(datadir);
toolsdir = fullfile(rootdir,'tools_analysis');
ftDir = fullfile(toolsdir,'..','..','..','Toolboxes','fieldtrip');
addpath(ftDir); ft_defaults;
addpath(genpath(fullfile(toolsdir,'lib')))

erpFigDir = fullfile(datadir, 'erpFigs');
if ~exist(erpFigDir); mkdir(erpFigDir); end

% we'll loop through subject data dirs, so just get the path from there
inputFileName = ['Preprocess' filesep 'timelocked_averages.mat'];

% since our folders are all e.g. S01, S12 etc, let's load information about
% the directories that have that pattern
tmp = dir([datadir filesep 'S*']);
subjectFolders = tmp([tmp.isdir]); clear tmp
% and loop through them
for subjectNum = 1:numel(subjectFolders)
    clear thisFile

    % full path to the inputfile
    thisFile = fullfile(subjectFolders(subjectNum).folder,subjectFolders(subjectNum).name,inputFileName);
    % skip this loop if we can't find one
    if ~exist(thisFile,'file'); continue; end
    
    disp('this is subject:')
    disp(subjectFolders(subjectNum).name);
    
     % grab info about all the meeg data files we care about
    theseFiles = dir([subjectFolders(subjectNum).folder filesep subjectFolders(subjectNum).name filesep inputFileName]);
    if isempty(theseFiles); continue; end
    
    %% so here, load what you care about and play with them!
    disp('loading')
    whichVars = {'ec_responseLockedAverage' 'hc_responseLockedAverage' 'er_responseLockedAverage' 'hr_responseLockedAverage'};
    
    data{subjectNum} = load(thisFile,whichVars{:});
    
end; clear theseFiles thisFile subjectFolders subjectNum

% analysis plan:
%
% RQ: is there a difference between ramping activity related to perceptual
%     processing and activity related to categorisation?
% is there ramping activity that scales with the difficulty of the
% perceptual difficulty for coherence as it does for categorisation
%%
% CPP or P300 (good justification here: https://www.nature.com/articles/s41562-020-0863-4 and here: https://www.sciencedirect.com/science/article/pii/S0006899319301386)
%        - select the Pz, or some averaged cluster centred on Pz (e.g. the two
%          closest neighbours)
% can also compare the motor prep using bilateral beta across conditions,
%     but the specifics of this I don't know yet->I would like to do at least
%     one frequency-based analysis
%
%
% 1. first subtract easy from hard seperately for coherence and rule to
%     determine which sensors are the most sensitive to this effect. So we
%     will have:
%         - coherence difficulty sensors
%         - rule difficulty sensors
%         - then we ask how are these two sensor groups different in space?
% 2. on these sensors and the classic CPP or P300, compare the amplitude and the slope across the four
%     conditions
%         - to compare amplitude, test mean amplitude within a window centred on the
%           response locked grand average peak
%         - to compare slope, test values from roughly 400ms prior to response
%             - alternatively, if I have time, I can fit a straight line to
%               the response locked waveform, but this seems hard
%         - so just compare at every timepoint across the conditions with a
%           bayes analysis
%             - to do that, just write about them OR do the difference of
%             the differences (coherence effect-categorisation effect) BUT
%             ONLY IF THEY LOOK DIFFERENT
%         - can select MEG equivs and also google it and do the exact same
%         analysis but with the GRADs
%         - so now we have are cat and coh comparable in terms of
%         timecourse
%% now is the topography equivalent across the two

% 3. compare easy vs hard seperately for coherence and rule for every sensor
%     at every time point (cluster test) to see if we get something similar
%         - is this the same as amplitude, if we have a range of times at
%           which the value is significantly different?
%         - are coherence clusters happening in the same place as coherence
%           sensitivity?
%         - are rule clusters happening in the same place as rule
%           sensitivity?
%         - how does the timecourse compare in coherence and rule as per
%           the above?

% descriptive -> similar or different
%%
%
% 4. find when the variability in the CPP correlate with
%     the model parameters
%         - sliding window, as per alessandro's paper
%         - but also let's do coherence only and rule only model parameters
%         - so theoretical trial-wise signal e.g. (boundary-bias/2)/(RT-non-decision time)
%         - instead of this, we can correlate the difference in the drift
%           rate between e.g. high/low coherence with the difference in CPP
%           across participants---so a bigger difference in drift rate = a
%           bigger difference in neural signal
%         - do this across timepoints and we can say they correlate at the
%         times as in the timecourse analysis
%         - then we can do space after we see the topology analysis to
%           constrain _where_ we look
%% 5. then rdm
% model fit through time for:
%   - stimulus representation fit through time
%   - categorisation representation fit through time
%   - for high coherence and low coherence for example
%   - and high or low categorisation
%   - we'll do it whole brain decoding (mags+grads) and maybe also eeg
%   - and we can run the same correlations as with the models and CPP and
%   this would be our final analysis (so diff in decoding correlated to
%   diff in models)
% we should also do button presses because we will want time locked to stim
%   and response

%% and yes vary the models for c and r seperately to see if it distinguishes the important parameters


%% so we can do cluster tests like below
%%
%%

%% first calculate neighbours

% meg is standard---can just use FTs version
% megNeighbours = load(fullfile(ftDir,'template','neighbours','neuromag306mag_neighb.mat'));
% or make our own from a subject:
megLayout = fullfile(datadir,'S01','Preprocess','meg_layout.lay');
% for EEG, we have a modified 70channel cap to fit 64 channels, so although
% there is a template for 64 chan eeg at fullfile(ftDir,'template','neighbours','easycap64ch-avg_neighb.mat')
% I think better to make this ourselves from a participant layout
eegLayout = fullfile(datadir,'S01','Preprocess','eeg_layout.lay');
cfg = [];
cfg.channel = 'EEG';
cfg.method = 'distance';
cfg.layout = eegLayout;
eegNeighbours = ft_prepare_neighbours(cfg);

cfg = [];
cfg.channel = 'MEG';
cfg.method = 'distance';
cfg.layout = megLayout;
megNeighbours = ft_prepare_neighbours(cfg);

%% is easy coh different from hard coh

cfg = [];
ecAll = returnStructs(data, 'ec_responseLockedAverage');
ecAve = ft_timelockgrandaverage(cfg, ecAll{:});
cfg = [];
hcAll = returnStructs(data, 'hc_responseLockedAverage');
hcAve = ft_timelockgrandaverage(cfg, hcAll{:});
% and prep this for a plot of our clusters
cfg = [];
cfg.parameter = 'avg';
cfg.operation = 'subtract';
diffEcHc = ft_math(cfg, ecAve, hcAve);

cfg = [];
cfg.channel = {'EEG'};
cfg.latency = 'all';
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.neighbours = eegNeighbours;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 500;

numSubjs  = numel(ecAll);
design = zeros(2, numSubjs*2); % design is 2 x num subjects matrix
design(1,:) = [1:numSubjs 1:numSubjs]; % first row is which subject the data belongs to in the two datasets we are comparing
design(2,:) = [ones(1,numSubjs) ones(1,numSubjs)*2]; % second row is which dataset the data belongs to
cfg.design = design;
cfg.uvar = 1; % uvar is the row that the subject definition is on
cfg.ivar = 2; % ivar is the row that the dataset definition is on

eegStatCoh = ft_timelockstatistics(cfg, ecAll{:}, hcAll{:});
save(fullfile(datadir,'eegCohErpStat.mat'),'eegStatCoh')

% define parameters for plotting
timestep      = 0.05; %(in seconds)
sampling_rate = 1000;
sample_count  = length(eegStatCoh.time);
j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in M/EEG samples
% get relevant values
pos_cluster_pvals = [eegStatCoh.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.025);
pos       = ismember(eegStatCoh.posclusterslabelmat, pos_clust);
% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(diffEcHc.label, eegStatCoh.label);
% plot
for k = 1:16
   cfg.figure     = subplot(4,4,k);
   cfg.xlim       = [j(k) j(k+1)];
   pos_int        = zeros(numel(diffEcHc.label),1);
   pos_int(i1)    = all(pos(i2, m(k):m(k+1)), 2);
   cfg.highlight  = 'on';
   cfg.highlightchannel = find(pos_int);
   cfg.comment    = 'xlim';
   cfg.commentpos = 'title';
   cfg.layout     = eegLayout;
   cfg.figure     = 'gca';
   ft_topoplotER(cfg, diffEcHc);
end
f = gcf; f.Position = [10 10 1600 1600];
print([erpFigDir filesep 'eeg_coh_sig_ERP_clusters.png'], '-dpng');
close all


cfg.channel = {'MEG'};
cfg.neighbours = megNeighbours;
megStatCoh = ft_timelockstatistics(cfg, ecAll{:}, hcAll{:});
save(fullfile(datadir,'megCohErpStat.mat'),'megStatCoh')

% define parameters for plotting
timestep      = 0.05; %(in seconds)
sampling_rate = 1000;
sample_count  = length(megStatCoh.time);
j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in M/EEG samples
% get relevant values
pos_cluster_pvals = [megStatCoh.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.025);
pos       = ismember(megStatCoh.posclusterslabelmat, pos_clust);
% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(diffEcHc.label, megStatCoh.label);
% plot
for k = 1:16
   cfg.figure     = subplot(4,5,k);
   cfg.xlim       = [j(k) j(k+1)];
   pos_int        = zeros(numel(diffEcHc.label),1);
   pos_int(i1)    = all(pos(i2, m(k):m(k+1)), 2);
   cfg.highlight  = 'on';
   cfg.highlightchannel = find(pos_int);
   cfg.comment    = 'xlim';
   cfg.commentpos = 'title';
   cfg.layout     = megLayout;
   cfg.figure     = 'gca';
   ft_topoplotER(cfg, diffEcHc);
end
f = gcf; f.Position = [10 10 1600 1600];
print([erpFigDir filesep 'meg_coh_sig_ERP_clusters.png'], '-dpng');
close all


%% is easy rule different from hard rule?

cfg = [];
erAll = returnStructs(data, 'er_responseLockedAverage');
erAve = ft_timelockgrandaverage(cfg, erAll{:});
cfg = [];
hrAll = returnStructs(data, 'hr_responseLockedAverage');
hrAve = ft_timelockgrandaverage(cfg, hrAll{:});
% and prep this for a plot of our clusters
cfg = [];
cfg.parameter = 'avg';
cfg.operation = 'subtract';
diffErHr = ft_math(cfg, erAve, hrAve);

cfg = [];
cfg.channel = {'EEG'};
cfg.latency = 'all';
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.neighbours = eegNeighbours;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 500;

numSubjs  = numel(erAll);
design = zeros(2, numSubjs*2); % design is 2 x num subjects matrix
design(1,:) = [1:numSubjs 1:numSubjs]; % first row is which subject the data belongs to in the two datasets we are comparing
design(2,:) = [ones(1,numSubjs) ones(1,numSubjs)*2]; % second row is which dataset the data belongs to
cfg.design = design;
cfg.uvar = 1; % uvar is the row that the subject definition is on
cfg.ivar = 2; % ivar is the row that the dataset definition is on

eegStatRule = ft_timelockstatistics(cfg, erAll{:}, hrAll{:});
save(fullfile(datadir,'eegRuleErpStat.mat'),'eegStatRule')

% define parameters for plotting
timestep      = 0.05; %(in seconds)
sampling_rate = 1000;
sample_count  = length(eegStatRule.time);
j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in M/EEG samples
% get relevant values
pos_cluster_pvals = [eegStatRule.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.025);
pos       = ismember(eegStatRule.posclusterslabelmat, pos_clust);
% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(diffErHr.label, eegStatRule.label);
% plot
for k = 1:16
   cfg.figure     = subplot(4,4,k);
   cfg.xlim       = [j(k) j(k+1)];
   pos_int        = zeros(numel(diffErHr.label),1);
   pos_int(i1)    = all(pos(i2, m(k):m(k+1)), 2);
   cfg.highlight  = 'on';
   cfg.highlightchannel = find(pos_int);
   cfg.comment    = 'xlim';
   cfg.commentpos = 'title';
   cfg.layout     = eegLayout;
   cfg.figure     = 'gca';
   ft_topoplotER(cfg, diffErHr);
end
f = gcf; f.Position = [10 10 1600 1600];
print([erpFigDir filesep 'eeg_rule_sig_ERP_clusters.png'], '-dpng');
close all


cfg.channel = {'MEG'};
cfg.neighbours = megNeighbours;
megStatRule = ft_timelockstatistics(cfg, erAll{:}, hrAll{:});
save(fullfile(datadir,'megRuleErpStat.mat'),'megStatRule')


% define parameters for plotting
timestep      = 0.05; %(in seconds)
sampling_rate = 1000;
sample_count  = length(megStatRule.time);
j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in M/EEG samples
% get relevant values
pos_cluster_pvals = [megStatRule.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.025);
pos       = ismember(megStatRule.posclusterslabelmat, pos_clust);
% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(diffErHr.label, megStatRule.label);
% plot
for k = 1:16
   cfg.figure     = subplot(4,4,k);
   cfg.xlim       = [j(k) j(k+1)];
   pos_int        = zeros(numel(diffErHr.label),1);
   pos_int(i1)    = all(pos(i2, m(k):m(k+1)), 2);
   cfg.highlight  = 'on';
   cfg.highlightchannel = find(pos_int);
   cfg.comment    = 'xlim';
   cfg.commentpos = 'title';
   cfg.layout     = megLayout;
   cfg.figure     = 'gca';
   ft_topoplotER(cfg, diffErHr);
end
f = gcf; f.Position = [10 10 1600 1600];
print([erpFigDir filesep 'meg_rule_sig_ERP_clusters.png'], '-dpng');
close all




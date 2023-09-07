%% set up

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

toolbox = fullfile(rootdir,'..','..','Toolboxes','gramm'); addpath(toolbox); clear toolbox
toolbox = fullfile(rootdir,'..','..','Toolboxes','bayesFactor'); addpath(toolbox); clear toolbox

erpFigDir = fullfile(datadir, 'erpFigs');
if ~exist(erpFigDir); mkdir(erpFigDir); end

% we'll loop through subject data dirs, so just get the path from there
inputFileName = ['Preprocess' filesep 'timelocked_averages.mat'];

% some colours
teal = [0.2, 0.6, 0.7];
coral = [0.9, 0.4, 0.3];
lilac = [0.7, 0.5, 0.8];

%% load

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
    
    % so here, load what you care about and play with them!
    disp('loading')
    whichVars = {'ec_responseLockedAverage' 'hc_responseLockedAverage' 'er_responseLockedAverage' 'hr_responseLockedAverage' 'ec_coherenceLockedAverage' 'hc_coherenceLockedAverage' 'er_coherenceLockedAverage' 'hr_coherenceLockedAverage'};
    
    data{subjectNum} = load(thisFile,whichVars{:});
    
end; clear theseFiles thisFile subjectFolders subjectNum

%% get some layouts and calculate neighbours
%

% meg is standard---can just use FTs version
% megNeighbours = load(fullfile(ftDir,'template','neighbours','neuromag306mag_neighb.mat'));
% or make our own from a subject:
megLayout = fullfile(datadir,'S01','Preprocess','meg_layout.lay');
% for EEG, we have a modified 70channel cap to fit 64 channels, so although
% there is a template for 64 chan eeg at fullfile(ftDir,'template','neighbours','easycap64ch-avg_neighb.mat')
% I think better to make this ourselves from a participant layout
eegLayout = fullfile(datadir,'S01','Preprocess','eeg_layout.lay');

% now neighbours

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

% and we'll collect useful sensors

CPP = {'EEG040' 'EEG041' 'EEG042'};

%% grab the data and average it
%

% coherence
cfg = [];
ecRespAll = returnStructs(data, 'ec_responseLockedAverage');
ecRespAve = ft_timelockgrandaverage(cfg, ecRespAll{:});
cfg = [];
ecOnsAll = returnStructs(data, 'ec_coherenceLockedAverage');
ecOnsAve = ft_timelockgrandaverage(cfg, ecOnsAll{:});
cfg = [];
hcRespAll = returnStructs(data, 'hc_responseLockedAverage');
hcRespAve = ft_timelockgrandaverage(cfg, hcRespAll{:});
cfg = [];
hcOnsAll = returnStructs(data, 'hc_coherenceLockedAverage');
hcOnsAve = ft_timelockgrandaverage(cfg, hcOnsAll{:});

%
% categorisation
cfg = [];
erRespAll = returnStructs(data, 'er_responseLockedAverage');
erRespAve = ft_timelockgrandaverage(cfg, erRespAll{:});
cfg = [];
erOnsAll = returnStructs(data, 'er_coherenceLockedAverage');
erOnsAve = ft_timelockgrandaverage(cfg, erOnsAll{:});
cfg = [];
hrRespAll = returnStructs(data, 'hr_responseLockedAverage');
hrRespAve = ft_timelockgrandaverage(cfg, hrRespAll{:});
cfg = [];
hrOnsAll = returnStructs(data, 'hr_coherenceLockedAverage');
hrOnsAve = ft_timelockgrandaverage(cfg, hrOnsAll{:});

%% now get the differences

for subject = 1:numel(ecRespAll)
    fprintf('\nthis is subject %.0f of %.0f\n\n',subject,numel(ecRespAll))
    
    cOnsDiffAll{subject} = getDifference(ecOnsAll{subject}, hcOnsAll{subject});
    cRespDiffAll{subject}= getDifference(ecRespAll{subject}, hcRespAll{subject});
    rOnsDiffAll{subject} = getDifference(erOnsAll{subject}, hrOnsAll{subject});
    rRespDiffAll{subject}= getDifference(erRespAll{subject}, hrRespAll{subject});

end; clear subject

cfg = [];
cOnsDiffAve = ft_timelockgrandaverage(cfg, cOnsDiffAll{:});
cRespDiffAve = ft_timelockgrandaverage(cfg, cRespDiffAll{:});
rOnsDiffAve = ft_timelockgrandaverage(cfg, rOnsDiffAll{:});
rRespDiffAve = ft_timelockgrandaverage(cfg, rRespDiffAll{:});

%% and bayes test the differences

[cOnsDiffAve.bfs cOnsDiffAve.reports cOnsDiffAve.code] = testDiffsAcrossTime(cOnsDiffAll,CPP);
[cRespDiffAve.bfs cRespDiffAve.reports cRespDiffAve.code] = testDiffsAcrossTime(cRespDiffAll,CPP);
[rOnsDiffAve.bfs rOnsDiffAve.reports rOnsDiffAve.code] = testDiffsAcrossTime(rOnsDiffAll,CPP);
[rRespDiffAve.bfs rRespDiffAve.reports rRespDiffAve.code] = testDiffsAcrossTime(rRespDiffAll,CPP);


%% the timecourse of univariate activity
%

% coh onset in eeg
plotTimecourse(eegLayout,ecOnsAve,hcOnsAve,cOnsDiffAve,...
    {'easy coh';'hard coh'},'northwest','Coherence Onset')

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = eegLayout;
% cfg.xlim       = [-0.2 0.6];
% cfg.ylim       = [-1.6e-04 -1.46e-04];
cfg.channel    = CPP;
% cfg.channel    = 'EEG';
cfg.linecolor = [teal; coral];
figure;
ft_singleplotER(cfg, ecOnsAve, hcOnsAve);
line([0 0], get(gca,'YLim'), 'Color', 'k', 'LineStyle', '--')
ylabel('Mean EEG Amplitude (uV)')
xlabel('Time (s)')
% figure;
% ft_singleplotER(cfg, cOnsDiffAve);
yyaxis right;
plot(cOnsDiffAve.time,cOnsDiffAve.code,...
    'LineStyle', 'none', 'Marker', '.', 'Color', lilac)
ylabel('Bayes Evidence for Difference')
ylim([0 10])
yticks([0 1 2])
yticklabels({'weak (<3)' 'moderate (3-10)' 'strong (10+)'})
ax = gca;
ax.YColor = lilac;
clear ax
legend({'easy coh';'hard coh'}, 'Location', 'northwest');
title('CPP (CP1 CPz CP2) in EEG: Coherence Onset', 'FontSize', 14)

% coh response in eeg
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = eegLayout;
% cfg.xlim       = [-0.2 0.6];
% cfg.ylim       = [-1.6e-04 -1.46e-04];
cfg.channel    = CPP;
cfg.linecolor = [teal; coral];
figure;
ft_singleplotER(cfg, ecRespAve, hcRespAve);
line([0 0], get(gca,'YLim'), 'Color', 'k', 'LineStyle', '--')
ylabel('Mean EEG Amplitude (uV)')
xlabel('Time (s)')
% figure;
% ft_singleplotER(cfg, cRespDiffAve);
yyaxis right;
plot(cRespDiffAve.time,cRespDiffAve.code,...
    'LineStyle', 'none', 'Marker', '.', 'Color', lilac)
ylabel('Bayes Evidence for Difference')
ylim([0 10])
yticks([0 1 2])
yticklabels({'weak (<3)' 'moderate (3-10)' 'strong (10+)'})
ax = gca;
ax.YColor = lilac;
clear ax
legend({'easy coh';'hard coh'}, 'Location', 'northwest');
title('CPP (CP1 CPz CP2) in EEG: Coherence Response', 'FontSize', 14)

% rule onset in eeg
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = eegLayout;
% cfg.xlim       = [-0.2 0.6];
% cfg.ylim       = [-1.6e-04 -1.46e-04];
cfg.channel    = CPP;
cfg.linecolor = [teal; coral];
figure;
ft_singleplotER(cfg, erOnsAve, hrOnsAve);
line([0 0], get(gca,'YLim'), 'Color', 'k', 'LineStyle', '--')
ylabel('Mean EEG Amplitude (uV)')
xlabel('Time (s)')
% figure;
% ft_singleplotER(cfg, rOnsDiffAve);
yyaxis right;
plot(rOnsDiffAve.time,rOnsDiffAve.code,...
    'LineStyle', 'none', 'Marker', '.', 'Color', lilac)
ylabel('Bayes Evidence for Difference')
ylim([0 10])
yticks([0 1 2])
yticklabels({'weak (<3)' 'moderate (3-10)' 'strong (10+)'})
ax = gca;
ax.YColor = lilac;
clear ax
legend({'easy rule';'hard rule'}, 'Location', 'southeast');
title('CPP (CP1 CPz CP2) in EEG: Categorisation Onset', 'FontSize', 14)

% rule response in eeg
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = eegLayout;
% cfg.xlim       = [-0.2 0.6];
% cfg.ylim       = [-1.6e-04 -1.46e-04];
cfg.channel    = CPP;
cfg.linecolor = [teal; coral];
figure;
ft_singleplotER(cfg, erRespAve, hrRespAve);
line([0 0], get(gca,'YLim'), 'Color', 'k', 'LineStyle', '--')
ylabel('Mean EEG Amplitude (uV)')
xlabel('Time (s)')
% figure;
% ft_singleplotER(cfg, rRespDiffAve);
yyaxis right;
plot(rOnsDiffAve.time,rOnsDiffAve.code,...
    'LineStyle', 'none', 'Marker', '.', 'Color', lilac)
ylabel('Bayes Evidence for Difference')
ylim([0 10])
yticks([0 1 2])
yticklabels({'weak (<3)' 'moderate (3-10)' 'strong (10+)'})
ax = gca;
ax.YColor = lilac;
clear ax
legend({'easy rule';'hard rule'}, 'Location', 'northwest');
title('CPP (CP1 CPz CP2) in EEG: Categorisation Response', 'FontSize', 14)



%% the topology of univariate activity
%

% where is easy coh different from hard coh

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


% where is easy rule different from hard rule?

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

function difference = getDifference(data1,data2)

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
difference = ft_math(cfg, data1, data2);

return
end

function figHandle = plotSensor(structure,channelRow,time)
% structure should be a fieldtrip data structure
% channelRow should be the row number of the channel you want
% time should be the range of timepoints you want

selected_data = structure.avg(channelRow,time); % MLC24 is the 9th channel, -0.2 to 1.0 is sample 241 to 601
selected_time = structure.time(time);
figure;
figHandle = plot(selected_time, selected_data)

return
end

function [bfs report code] = testDiffsAcrossTime(dataStruct,channels)

for i = 1:numel(dataStruct)
    diffs(i,:) = mean(dataStruct{i}.avg(find(ismember(dataStruct{i}.label,channels)),:));
end; clear i

for i = 1:numel(dataStruct{1}.time)
    bfs(i) = bf.ttest(diffs(:,i));
    if bfs(i) <= 3
        report{i} = 'weak';
        code(i) = 0;
    elseif (bfs(i)>3) && (bfs(i)<=10)
        report{i} = 'moderate';
        code(i) = 1;
    elseif bfs(i)>10
        report{i} = 'strong';
        code(i) = 2;
    end
end; clear i


return
end


function plotTimecourse(layout,data1,data2,data3,leg,legendloc,theTitle)

% coh onset in eeg
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
% cfg.xlim       = [-0.2 0.6];
cfg.ylim       = [-1.6e-04 -1.46e-04];
cfg.channel    = CPP;
% cfg.channel    = 'EEG';
cfg.linecolor = [teal; coral];
figure;
ft_singleplotER(cfg, data1, data2);
line([0 0], get(gca,'YLim'), 'Color', 'k', 'LineStyle', '--')
ylabel('Mean EEG Amplitude (uV)')
xlabel('Time (s)')
% figure;
% ft_singleplotER(cfg, cOnsDiffAve);
yyaxis right;
plot(data3.time,data3.code,...
    'LineStyle', 'none', 'Marker', '.', 'Color', lilac)
ylabel('Bayes Evidence for Difference')
ylim([0 10])
yticks([0 1 2])
yticklabels({'weak (<3)' 'moderate (3-10)' 'strong (10+)'})
ax = gca;
ax.YColor = lilac;
clear ax
legend(leg, 'Location', legendloc);
title(['CPP (CP1 CPz CP2) in EEG: ' theTitle], 'FontSize', 14)

return
end


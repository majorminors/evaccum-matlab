%%%%%%%%%%%%
%% set up %%
%%%%%%%%%%%%

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
peach = [0.9, 0.7, 0.4];


%%%%%%%%%%%%%%%%%%%%%%%%
%% load and prep data %%
%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    whichVars = {...
        'ec_responseLockedAverage' 'hc_responseLockedAverage'...
        'er_responseLockedAverage' 'hr_responseLockedAverage'...
        'ec_coherenceLockedAverage' 'hc_coherenceLockedAverage'...
        'er_coherenceLockedAverage' 'hr_coherenceLockedAverage'...
        };

%     whichVars = {...
%         'ecer_responseLockedAverage' 'ecer_coherenceLockedAverage'...
%         'echr_responseLockedAverage' 'echr_coherenceLockedAverage'...
%         'hcer_responseLockedAverage' 'hcer_coherenceLockedAverage'...
%         'hchr_responseLockedAverage' 'hchr_coherenceLockedAverage'...
%         };

    
    data{subjectNum} = load(thisFile,whichVars{:});
    
    disp('loaded')
    
end; clear theseFiles thisFile subjectFolders subjectNum

disp('loading complete')

% get some layouts and calculate neighbours

disp('prep layouts and neighbours')

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pull coherence data and average it %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('averaging coherence data')

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

disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pull categorisation data and average it %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('averaging categorisation data')

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

disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pull conditionwise data and average it %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('averaging conditionwise data')

cfg = [];
ecerRespAll = returnStructs(data, 'ecer_responseLockedAverage');
ecerRespAve = ft_timelockgrandaverage(cfg, ecerRespAll{:});
cfg = [];
ecerOnsAll = returnStructs(data, 'ecer_coherenceLockedAverage');
ecerOnsAve = ft_timelockgrandaverage(cfg, ecerOnsAll{:});
cfg = [];
echrRespAll = returnStructs(data, 'echr_responseLockedAverage');
echrRespAve = ft_timelockgrandaverage(cfg, echrRespAll{:});
cfg = [];
echrOnsAll = returnStructs(data, 'echr_coherenceLockedAverage');
echrOnsAve = ft_timelockgrandaverage(cfg, echrOnsAll{:});
cfg = [];
hcerRespAll = returnStructs(data, 'hcer_responseLockedAverage');
hcerRespAve = ft_timelockgrandaverage(cfg, hcerRespAll{:});
cfg = [];
hcerOnsAll = returnStructs(data, 'hcer_coherenceLockedAverage');
hcerOnsAve = ft_timelockgrandaverage(cfg, hcerOnsAll{:});
cfg = [];
hchrRespAll = returnStructs(data, 'hchr_responseLockedAverage');
hchrRespAve = ft_timelockgrandaverage(cfg, hchrRespAll{:});
cfg = [];
hchrOnsAll = returnStructs(data, 'hchr_coherenceLockedAverage');
hchrOnsAve = ft_timelockgrandaverage(cfg, hchrOnsAll{:});

disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get the subjectwise differences between manipulations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% and bayes test the differences

[cOnsDiffAve.bfs cOnsDiffAve.reports cOnsDiffAve.code] = testDiffsAcrossTime(cOnsDiffAll,CPP);
[cRespDiffAve.bfs cRespDiffAve.reports cRespDiffAve.code] = testDiffsAcrossTime(cRespDiffAll,CPP);
[rOnsDiffAve.bfs rOnsDiffAve.reports rOnsDiffAve.code] = testDiffsAcrossTime(rOnsDiffAll,CPP);
[rRespDiffAve.bfs rRespDiffAve.reports rRespDiffAve.code] = testDiffsAcrossTime(rRespDiffAll,CPP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the timecourse of univariate activity %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plotTimecourse(layout,dataAv1,dataAv2,dataDiff,dataAll1,dataAll2,xlims,ylims,channels,leg,legendloc,theTitle)

%% coh onset in eeg
plotTimecourse(eegLayout,ecOnsAve,hcOnsAve,cOnsDiffAve,...
    ecOnsAll,hcOnsAll,...
    [-0.2 1.5],[],CPP,...
    {'easy coh';'hard coh'},'northwest','Coherence Onset')
%% coh response in eeg
plotTimecourse(eegLayout,ecRespAve,hcRespAve,cRespDiffAve,...
    ecRespAll,hcRespAll,...
    [],[],CPP,...
    {'easy coh';'hard coh'},'northwest','Coherence Response')
%% rule onset in eeg
plotTimecourse(eegLayout,erOnsAve,hrOnsAve,rOnsDiffAve,...
    erOnsAll,hrOnsAll,...
    [-0.2 1.5],[-3.5e-06 7.5e-06],CPP,...
    {'easy cat';'hard cat'},'southeast','Categorisation Onset')
%% rule response in eeg
plotTimecourse(eegLayout,erRespAve,hrRespAve,rRespDiffAve,...
    erRespAll,hrRespAll,...
    [],[-3.5e-06 7.5e-06],CPP,...
    {'easy cat';'hard cat'},'northwest','Categorisation Response')

%% onset of all conditions in eeg

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = eegLayout;
% cfg.xlim = [];
% cfg.ylim = [];
cfg.channel    = CPP;
cfg.linecolor = [teal; lilac; coral; peach];

h = figure;
theseVars = {ecerOnsAve, echrOnsAve, hcerOnsAve, hchrOnsAve};
theseErrs{1} = getErrorFromStructs(ecerOnsAll,'error'); % or 'deviation';
theseErrs{2} = getErrorFromStructs(echrOnsAll,'error'); % or 'deviation';
theseErrs{3} = getErrorFromStructs(hcerOnsAll,'error'); % or 'deviation';
theseErrs{4} = getErrorFromStructs(hchrOnsAll,'error'); % or 'deviation';

ft_singleplotER(cfg, theseVars{:});
hold on
plotErrorFromStructs(h,theseVars,[teal; lilac; coral; peach-peach*0.2],CPP,theseErrs)
line([0 0], get(gca,'YLim'), 'Color', 'k', 'LineStyle', '--')
hold off;

ylabel('Mean EEG Amplitude (uV)')
xlabel('Time (s)')
legend({'EasyCoh EasyCat';'EasyCoh HardCat';'HardCoh EasyCat';'HardCoh HardCat';},...
    'Location', 'southeast');
title('CPP (CP1 CPz CP2) in EEG: Onset in all Conditions', 'FontSize', 14)

%% response of all conditions in eeg

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = eegLayout;
% cfg.xlim = [];
% cfg.ylim = [];
cfg.channel    = CPP;
cfg.linecolor = [teal; lilac; coral; peach];

h = figure;
theseVars = {ecerRespAve, echrRespAve, hcerRespAve, hchrRespAve};
theseErrs{1} = getErrorFromStructs(ecerRespAll,'error'); % or 'deviation';
theseErrs{2} = getErrorFromStructs(echrRespAll,'error'); % or 'deviation';
theseErrs{3} = getErrorFromStructs(hcerRespAll,'error'); % or 'deviation';
theseErrs{4} = getErrorFromStructs(hchrRespAll,'error'); % or 'deviation';

ft_singleplotER(cfg, theseVars{:});
hold on
plotErrorFromStructs(h,theseVars,[teal; lilac; coral; peach-peach*0.2],CPP,theseErrs)
line([0 0], get(gca,'YLim'), 'Color', 'k', 'LineStyle', '--')
hold off;

ylabel('Mean EEG Amplitude (uV)')
xlabel('Time (s)')
legend({'EasyCoh EasyCat';'EasyCoh HardCat';'HardCoh EasyCat';'HardCoh HardCat';},...
    'Location', 'southeast');
title('CPP (CP1 CPz CP2) in EEG: Response in all Conditions', 'FontSize', 14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the topology of univariate activity %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% where is easy coh different from hard coh

% % we could get the difference between the grand average to plot our clusters
% % on top of
% cfg = [];
% cfg.parameter = 'avg';
% cfg.operation = 'subtract';
% diffEcHc = ft_math(cfg, ecRespAve, hcRespAve);
% but we already have the subjectwise differences, so we can also use them

% first set up the cluster test settings
cfg = [];
cfg.channel = {'EEG'};
cfg.latency = 'all';
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT'; % dependent samples t test (also independent: indepsamplesT)
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05; % alpha level of the sample-specific test statistic that will be used for thresholding
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2; % min number of neighbours required to be considered for clustering
cfg.neighbours = eegNeighbours;
cfg.tail = 0; % one or two sided test? (one would be -1 or 1)
cfg.clustertail = 0;
cfg.alpha = 0.025; % alpha for the permutation
cfg.numrandomization = 500; % how many draws?

% now let's set up the design
numSubjs  = numel(ecRespAll);
design = zeros(2, numSubjs*2); % design is 2 x num subjects matrix
design(1,:) = [1:numSubjs 1:numSubjs]; % first row is which subject the data belongs to in the two datasets we are comparing
design(2,:) = [ones(1,numSubjs) ones(1,numSubjs)*2]; % second row is which dataset the data belongs to
% so we're contrasting two things, and so we repeat subjects twice, then do
% a row of ones and a row of twos
cfg.design = design;
cfg.uvar = 1; % uvar is the row that the subject definition is on
cfg.ivar = 2; % ivar is the row that the dataset (independent variables) definition is on

% get the stats
eegStatCohResp = ft_timelockstatistics(cfg, ecRespAll{:}, hcRespAll{:});
save(fullfile(datadir,'eegCohRespErpStat.mat'),'eegStatCohResp')

% define parameters for plotting
timestep      = 0.05; %(in seconds)
sampling_rate = 1000;
sample_count  = length(eegStatCohResp.time);
j = [0:timestep:0.8];   % in secs, get some timings for each subplot
m = [1:timestep*sampling_rate:sample_count];  % in secs, what timings according to the M/EEG samples
% get relevant values
pos_cluster_pvals = [eegStatCohResp.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.025);
pos       = ismember(eegStatCohResp.posclusterslabelmat, pos_clust);
neg_cluster_pvals = [eegStatCohResp.negclusters(:).prob];
neg_clust = find(neg_cluster_pvals < 0.025);
neg       = ismember(eegStatCohResp.negclusterslabelmat, neg_clust);% ensure the channels to have the same order in the average and in the statistical output
% which might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(cRespDiffAve.label, eegStatCohResp.label);
% plot
for k = 1:16
   cfg.figure     = subplot(4,4,k);
   cfg.xlim       = [j(k) j(k+1)]; % set the xlimits to timepoints
   % now we create an index to the channels that are in clusters to be plotted,
   % by timepoint
   pos_int        = zeros(numel(cRespDiffAve.label),1);
   pos_int(i1)    = all(pos(i2, m(k):m(k+1)), 2);
   neg_int        = zeros(numel(cRespDiffAve.label),1);
   neg_int(i1)    = all(neg(i2, m(k):m(k+1)), 2);
   cfg.highlight  = 'on';
   % now we pull the indices we created earlier
   cfg.highlightchannel = find(pos_int | neg_int);
   cfg.comment    = 'xlim';
   cfg.commentpos = 'title';
   cfg.layout     = eegLayout;
   cfg.interactive = 'no';
   cfg.figure     = 'gca';
%    ft_topoplotER(cfg, cRespDiffAve);
end
f = gcf; f.Position = [10 10 1600 1600];
print([erpFigDir filesep 'eeg_coh_resp_sig_ERP_clusters.png'], '-dpng');
% close all



function difference = getDifference(data1,data2)

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
difference = ft_math(cfg, data1, data2);

return
end

function figHandle = plotSensors(structure,channelRow,time)
% structure should be a fieldtrip data structure
% channelRow should be the row number of the channel you want
% time should be the range of timepoints you want

selected_data = structure.avg(channelRow,time);
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


function plotTimecourse(layout,dataAv1,dataAv2,dataDiff,dataAll1,dataAll2,xlims,ylims,channels,leg,legendloc,theTitle)

teal = [0.2, 0.6, 0.7];
coral = [0.9, 0.4, 0.3];
lilac = [0.7, 0.5, 0.8];

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
% if ~isempty(xlims); cfg.xlim = xlims; end
if ~isempty(ylims); cfg.ylim = ylims; end
cfg.channel    = channels;
% cfg.channel    = 'EEG';
cfg.linecolor = [teal; coral];

theseVars = {dataAv1, dataAv2};
theseErrs{1} = getErrorFromStructs(dataAll1,'error'); % or 'deviation';
theseErrs{2} = getErrorFromStructs(dataAll2,'error'); % or 'deviation';

h = figure;
ft_singleplotER(cfg, theseVars{:});
hold on
plotErrorFromStructs(h,theseVars,[teal; coral],channels,theseErrs)
line([0 0], get(gca,'YLim'), 'Color', 'k', 'LineStyle', '--')
if ~isempty(xlims); xlim(xlims); end
hold off;

ylabel('Mean EEG Amplitude (uV)')
xlabel('Time (s)')
% figure;
% ft_singleplotER(cfg, cOnsDiffAve);
yyaxis right;
plot(dataDiff.time,dataDiff.code,...
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




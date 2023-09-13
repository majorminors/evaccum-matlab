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
% toolbox = fullfile(rootdir,'..','..','Toolboxes','bayesFactor'); addpath(toolbox); clear toolbox
toolbox = fullfile(rootdir,'..','..','Toolboxes','BFF_repo'); addpath(genpath(toolbox)); clear toolbox

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

cOnsDiffAve = getDifference(ecOnsAve,hcOnsAve);
cRespDiffAve = getDifference(ecRespAve,hcRespAve);

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

rOnsDiffAve = getDifference(erOnsAve,hrOnsAve);
rRespDiffAve = getDifference(ecRespAve,hcRespAve);

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
    fprintf('\ngetting subjectwise differences for subject %.0f of %.0f\n\n',subject,numel(ecRespAll))
    
    cOnsDiffAll{subject} = getDifference(ecOnsAll{subject}, hcOnsAll{subject});
    cRespDiffAll{subject}= getDifference(ecRespAll{subject}, hcRespAll{subject});
    rOnsDiffAll{subject} = getDifference(erOnsAll{subject}, hrOnsAll{subject});
    rRespDiffAll{subject}= getDifference(erRespAll{subject}, hrRespAll{subject});

end; clear subject

disp('getting averages')

% cfg = [];
% cOnsDiffAve = ft_timelockgrandaverage(cfg, cOnsDiffAll{:});
% cRespDiffAve = ft_timelockgrandaverage(cfg, cRespDiffAll{:});
% rOnsDiffAve = ft_timelockgrandaverage(cfg, rOnsDiffAll{:});
% rRespDiffAve = ft_timelockgrandaverage(cfg, rRespDiffAll{:});


%% and bayes test the differences

disp('bayes testing differences')

cOnsDiffAve.bfs = testDiffsAcrossTime(cOnsDiffAll,CPP);
cRespDiffAve.bfs = testDiffsAcrossTime(cRespDiffAll,CPP);
rOnsDiffAve.bfs = testDiffsAcrossTime(rOnsDiffAll,CPP);
rRespDiffAve.bfs = testDiffsAcrossTime(rRespDiffAll,CPP);

disp('done bayes testing differences')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the timecourse of univariate activity %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plotTimecourse(layout,dataAv1,dataAv2,dataDiff,dataAll1,dataAll2,xlims,ylims,channels,leg,legendloc,theTitle)

%% coh onset in eeg
plotTimecourse(eegLayout,ecOnsAve,hcOnsAve,cOnsDiffAve,...
    ecOnsAll,hcOnsAll,...
    [-0.2 1.5],[],CPP,...
    {'easy coh';'hard coh'},'northwest','Coherence Onset',...
    [erpFigDir filesep 'eeg_coh_ons_ERP.png'])
%% coh response in eeg
plotTimecourse(eegLayout,ecRespAve,hcRespAve,cRespDiffAve,...
    ecRespAll,hcRespAll,...
    [],[],CPP,...
    {'easy coh';'hard coh'},'northwest','Coherence Response',...
    [erpFigDir filesep 'eeg_coh_resp_ERP.png'])
%% rule onset in eeg
plotTimecourse(eegLayout,erOnsAve,hrOnsAve,rOnsDiffAve,...
    erOnsAll,hrOnsAll,...
    [-0.2 1.5],[-3.5e-06 7.5e-06],CPP,...
    {'easy cat';'hard cat'},'southeast','Categorisation Onset',...
    [erpFigDir filesep 'eeg_cat_ons_ERP.png'])
%% rule response in eeg
plotTimecourse(eegLayout,erRespAve,hrRespAve,rRespDiffAve,...
    erRespAll,hrRespAll,...
    [],[-3.5e-06 7.5e-06],CPP,...
    {'easy cat';'hard cat'},'northwest','Categorisation Response',...
    [erpFigDir filesep 'eeg_cat_resp_ERP.png'])
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
theseStructs = {ecerOnsAve, echrOnsAve, hcerOnsAve, hchrOnsAve};
theseErrs{1} = getErrorFromStructs(ecerOnsAll);
theseErrs{2} = getErrorFromStructs(echrOnsAll);
theseErrs{3} = getErrorFromStructs(hcerOnsAll);
theseErrs{4} = getErrorFromStructs(hchrOnsAll);

ft_singleplotER(cfg, theseStructs{:});
hold on
plotErrorFromStructs(theseStructs,[teal; lilac; coral; peach-peach*0.2],CPP,theseErrs)
line([0 0], get(gca,'YLim'), 'Color', 'k', 'LineStyle', '--')
hold off;

ylabel('Mean EEG Amplitude (uV)')
xlabel('Time (s)')
legend({'EasyCoh EasyCat';'EasyCoh HardCat';'HardCoh EasyCat';'HardCoh HardCat';},...
    'Location', 'southeast');
title('CPP (CP1 CPz CP2) in EEG: Onset in all Conditions', 'FontSize', 14)
f = gcf; f.Position = [10 10 1600 1600];
print([erpFigDir filesep 'eeg_allConds_ons_ERP.png'], '-dpng');
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
theseStructs = {ecerRespAve, echrRespAve, hcerRespAve, hchrRespAve};
theseErrs{1} = getErrorFromStructs(ecerRespAll);
theseErrs{2} = getErrorFromStructs(echrRespAll);
theseErrs{3} = getErrorFromStructs(hcerRespAll);
theseErrs{4} = getErrorFromStructs(hchrRespAll);

ft_singleplotER(cfg, theseStructs{:});
hold on
plotErrorFromStructs(theseStructs,[teal; lilac; coral; peach-peach*0.2],CPP,theseErrs)
line([0 0], get(gca,'YLim'), 'Color', 'k', 'LineStyle', '--')
hold off;

ylabel('Mean EEG Amplitude (uV)')
xlabel('Time (s)')
legend({'EasyCoh EasyCat';'EasyCoh HardCat';'HardCoh EasyCat';'HardCoh HardCat';},...
    'Location', 'southeast');
title('CPP (CP1 CPz CP2) in EEG: Response in all Conditions', 'FontSize', 14)
f = gcf; f.Position = [10 10 1600 1600];
print([erpFigDir filesep 'eeg_allConds_resp_ERP.png'], '-dpng');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the topology of univariate activity %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% where is easy coh different from hard coh

% let's get some zlimits
        
minValue = @(x, pref) min(min(x.avg(startsWith(x.label, pref), :)));
maxValue = @(x, pref) max(max(x.avg(startsWith(x.label, pref), :)));

zlimEegOns = [...
    min(minValue(cOnsDiffAve,'EEG'),minValue(rOnsDiffAve,'EEG')),...
    max(maxValue(cOnsDiffAve,'EEG'),maxValue(rOnsDiffAve,'EEG'))...
    ];
zlimEegResp = [...
    min(minValue(cRespDiffAve,'EEG'),minValue(rRespDiffAve,'EEG')),...
    max(maxValue(cRespDiffAve,'EEG'),maxValue(rRespDiffAve,'EEG'))...
    ];
zlimMegOns = [...
    min(minValue(cOnsDiffAve,'MEG'),minValue(rOnsDiffAve,'MEG')),...
    max(maxValue(cOnsDiffAve,'MEG'),maxValue(rOnsDiffAve,'MEG'))...
    ];
zlimMegResp = [...
    min(minValue(cRespDiffAve,'MEG'),minValue(rRespDiffAve,'MEG')),...
    max(maxValue(cRespDiffAve,'MEG'),maxValue(rRespDiffAve,'MEG'))...
    ];


%% EEG topo of onset across coherence
doTopos({'EEG'},eegNeighbours,ecOnsAll,hcOnsAll,...
    zlimEegOns,0.05,[-0.3,1.5],[6,6],36,...
    cOnsDiffAve,eegLayout,[erpFigDir filesep 'eeg_coh_ons_sig_ERP_clusters.png'])
%% EEG topo of response across coherence
doTopos({'EEG'},eegNeighbours,ecRespAll,hcRespAll,...
    zlimEegResp,0.05,[-0.6,0.2],[4,4],16,...
    cRespDiffAve,eegLayout,[erpFigDir filesep 'eeg_coh_resp_sig_ERP_clusters.png'])
%% EEG topo of onset across categorisation
doTopos({'EEG'},eegNeighbours,erOnsAll,hrOnsAll,...
    zlimEegOns,0.05,[-0.3,1.5],[6,6],36,...
    rOnsDiffAve,eegLayout,[erpFigDir filesep 'eeg_cat_ons_sig_ERP_clusters.png'])
%% EEG topo of response across categorisation
doTopos({'EEG'},eegNeighbours,erRespAll,hrRespAll,...
    zlimEegResp,0.05,[-0.6,0.2],[4,4],16,...
    rRespDiffAve,eegLayout,[erpFigDir filesep 'eeg_cat_resp_sig_ERP_clusters.png'])


%% MEG topo of onset across coherence
doTopos({'MEG'},megNeighbours,ecOnsAll,hcOnsAll,...
    zlimMegOns,0.05,[-0.3,1.5],[6,6],36,...
    cOnsDiffAve,megLayout,[erpFigDir filesep 'meg_coh_ons_sig_ERP_clusters.png'])
%% MEG topo of response across coherence
doTopos({'MEG'},megNeighbours,ecRespAll,hcRespAll,...
    zlimMegResp,0.05,[-0.6,0.2],[4,4],16,...
    cRespDiffAve,megLayout,[erpFigDir filesep 'meg_coh_resp_sig_ERP_clusters.png'])
%% MEG topo of onset across categorisation
doTopos({'MEG'},megNeighbours,erOnsAll,hrOnsAll,...
    zlimMegOns,0.05,[-0.3,1.5],[6,6],36,...
    rOnsDiffAve,megLayout,[erpFigDir filesep 'meg_cat_ons_sig_ERP_clusters.png'])
%% MEG topo of response across categorisation
doTopos({'MEG'},megNeighbours,erRespAll,hrRespAll,...
    zlimMegResp,0.05,[-0.6,0.2],[4,4],16,...
    rRespDiffAve,megLayout,[erpFigDir filesep 'meg_cat_resp_sig_ERP_clusters.png'])


%% subplots

function doTopos(channels,neighbours,dataAll1,dataAll2,zlims,timestep,timebounds,subplots,totalplots,dataDiff,layout,savename)
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
stats = ft_timelockstatistics(cfg, dataAll1{:}, dataAll2{:});

% save(fullfile(datadir,'eegCohRespErpStat.mat'),'stats')

% define parameters for plotting
% timestep      = 0.05; %(in seconds)
sampling_rate = 1000;
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
[i1,i2] = match_str(dataDiff.label, stats.label);
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
   pos_int        = zeros(numel(dataDiff.label),1);
   pos_int(i1)    = all(pos(i2, m(k):m(k+1)), 2);
   neg_int        = zeros(numel(dataDiff.label),1);
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
%    if k == totalplots; cfg.colorbar = 'EastOutside'; end
   ft_topoplotER(cfg, dataDiff);
end

cbh = colorbar('Location', 'eastoutside', 'Position', [0.93 0.1 0.01 0.8]);
caxis([zlims(1),zlims(2)])
cbh.Units = 'normalized';
cbh.Limits = [zlims(1),zlims(2)];
% cbh.Label.String = 'Bayes Factor';
f = gcf; f.Position = [10 10 1600 1600];
print(savename, '-dpng');
% close all

end

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

function bfs = testDiffsAcrossTime(dataStruct,channels,complementary)

for i = 1:numel(dataStruct)
    diffs(:,i) = mean(dataStruct{i}.avg(find(ismember(dataStruct{i}.label,channels)),:));
end; clear i

% add the R module and get the path to Rscript
[status, result] = system('module add R && which Rscript');
if status == 0
    disp('R module added successfully');
    RscriptPath = strtrim(result);
    disp(['Rscript path: ' RscriptPath]);
else
    error('Failed to add R and/or locate Rscript');
end

% now run the rscript version of the bayes analysis
%   we can also get the bf for the complementary interval
%   by specifying complementary = 2. Let's set a default:
if ~exist('complementary','var'); complementary = 1; end
bfs = bayesfactor_R_wrapper(diffs,'Rpath',RscriptPath,'returnindex',complementary,...
    'args','mu=0,rscale="medium",nullInterval=c(0.5,Inf)');

% alternatively they hav implemented it in matlab
% [bfs, bfs_complementary_interval] = bayesfactor(diffs, 'interval',[-Inf Inf]);

return
end


function plotTimecourse(layout,dataAv1,dataAv2,dataDiff,dataAll1,dataAll2,xlims,ylims,channels,leg,legendloc,theTitle,savename)

teal = [0.2, 0.6, 0.7];
coral = [0.9, 0.4, 0.3];
lilac = [0.7, 0.5, 0.8];
colours = [teal;coral;lilac];

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = layout;
% if ~isempty(xlims); cfg.xlim = xlims; end
if ~isempty(ylims); cfg.ylim = ylims; end
cfg.channel    = channels;
% cfg.channel    = 'EEG';
cfg.linecolor = [teal; coral];

theseStructs = {dataAv1, dataAv2};
theseErrs{1} = getErrorFromStructs(dataAll1);
theseErrs{2} = getErrorFromStructs(dataAll2);

subplot(2, 1, 1);

set(gca, 'Position', [0.13 0.57 0.775 0.37]); % set position of the first subplot to make it taller
% ft_singleplotER(cfg, theseStructs{:});
for i = 1:numel(theseStructs)
        plot(theseStructs{i}.time, mean(theseStructs{i}.avg(find(ismember(theseStructs{i}.label,channels)),:)),...
            'Color',colours(i,:))
    hold on
end; clear i

hold on
plotErrorFromStructs(theseStructs,colours,channels,theseErrs)
line([0 0], get(gca,'YLim'), 'Color', 'k', 'LineStyle', '--')
% if ~isempty(xlims); xlim(xlims); end
hold off;
ylabel('Mean EEG Amplitude (uV)')
xlabel('Time (s)')
legend(leg, 'Location', legendloc);
title(['CPP (CP1 CPz CP2) in EEG: ' theTitle], 'FontSize', 14)

subplot(2, 1, 2);
bf = dataDiff.bfs;
load('color_bwr.mat'); % in BFF repo
exponential_minmax=5;
val_col_map = logspace(-exponential_minmax,exponential_minmax,size(co_bwr,1));
for t = 1:length(dataDiff.time)
    [~,idx] = min(abs(val_col_map-bf(t)));  
    st = stem(dataDiff.time(t),bf(t),'Clipping','off','basevalue',1,'Color','k','MarkerFaceColor',co_bwr(idx,1:3),'MarkerEdgeColor','k','LineWidth',1,'MarkerSize',6);
    hold on;

end
ax = gca;
set(ax,'YScale','log','XLim',[dataDiff.time(1),dataDiff.time(end)], ...
        'YLim',[1e-5 1e5],'YTick',10.^(-5:2:5))
xlabel('Time (s)')
ylabel('BF (log scale)')
colormap(co_bwr)
cbh = colorbar;
caxis([-exponential_minmax,exponential_minmax])
cbh.Units = 'normalized';
cbh.Limits = [-exponential_minmax,exponential_minmax];
cbh.Position(1) = 0.92;cbh.Position(3) = 0.01;cbh.Position(4) = ax.Position(4);cbh.Position(2) = ax.Position(2);
cbh.Label.String = 'Bayes Factor';
cbh.TickLabels=arrayfun(@(x) ['10^{' num2str(x) '}'], cbh.Ticks, 'UniformOutput', false);
f = gcf; f.Position = [10 10 1600 1600];
print(savename, '-dpng');


return
end




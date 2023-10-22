%%%%%%%%%%%%
%% set up %%
%%%%%%%%%%%%

clear all %#ok

% set up paths
addpath /hpc-software/matlab/cbu/
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/');

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';
datadir = fullfile(rootdir,'data','meg_pilot_4'); addpath(datadir);
toolsdir = fullfile(rootdir,'tools_analysis');
ftDir = fullfile(toolsdir,'..','..','..','Toolboxes','fieldtrip');
jobdir = fullfile(rootdir,'job_logging');
addpath(ftDir); ft_defaults;
addpath(genpath(fullfile(toolsdir,'lib')))

toolbox = fullfile(rootdir,'..','..','Toolboxes','gramm'); addpath(toolbox); clear toolbox
% toolbox = fullfile(rootdir,'..','..','Toolboxes','bayesFactor'); addpath(toolbox); clear toolbox
toolbox = fullfile(rootdir,'..','..','Toolboxes','BFF_repo'); addpath(genpath(toolbox)); clear toolbox

tfrFigDir = fullfile(datadir, 'tfrFigs');
if ~exist(tfrFigDir,'dir'); mkdir(tfrFigDir); end

% we'll loop through subject data dirs, so just get the path from there
inputFileName = ['Preprocess' filesep 'tfr_hanning.mat'];
% inputFileName = ['Preprocess' filesep 'tfr_multi.mat'];
if contains(inputFileName,'hanning')
    getVar = @(x) sprintf(x,'hann');
elseif contains(inputFileName,'multi')
    getVar = @(x) sprintf(x,'multi');
end

% some colours
teal = [0.2, 0.6, 0.7];
coral = [0.9, 0.4, 0.3];
lilac = [0.7, 0.5, 0.8];
maroon = [0.6579, 0.2821, 0.198];


%%%%%%%%%%%%%%%%%%%%%%%%
%% load and prep data %%
%%%%%%%%%%%%%%%%%%%%%%%%

% since our folders are all e.g. S01, S12 etc, let's load information about
% the directories that have that pattern
tmp = dir([datadir filesep 'S*']);
subjectFolders = tmp([tmp.isdir]); clear tmp
% and loop through them
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    disp('no parallel pool is currently initialised---initialising');
    workers = numel(subjectFolders);
    P = cbupool(workers);
    P.JobStorageLocation = jobdir;
    tempPath = fullfile(jobdir,'tmp');
    if ~exist(tempPath,'dir');mkdir(tempPath);end
    P.JobStorageLocation = tempPath;
    parpool(P,workers);
else
    disp('parallel pool has already been initialized---skipping');
end
parfor subjectNum = 1:numel(subjectFolders)
    % clear thisFile % can't do this in a parfor loop
    
    % full path to the inputfile
    thisFile = fullfile(subjectFolders(subjectNum).folder,subjectFolders(subjectNum).name,inputFileName);
    % skip this loop if we can't find one
    if ~exist(thisFile,'file'); continue; end
    % I seem to be doing a lot of checks in this loop and I can't remember why
    % so let's do another one to
    % account for the fact I can't clear the thisFile variable in a parfor
    % loop. let's check the subject directory part of thisFile matches the
    % subject number
    pathParts = strsplit(thisFile, '/');
    index = find(contains(pathParts, 'S'));
    subjectCode = pathParts{index}; %#ok
    if ~strcmp(subjectCode,subjectFolders(subjectNum).name); error('file doesnt match subject'); end
    
    fprintf('this is subject: %s...',subjectFolders(subjectNum).name)
    
     % grab info about all the meeg data files we care about
    theseFiles = dir([subjectFolders(subjectNum).folder filesep subjectFolders(subjectNum).name filesep inputFileName]);
    if isempty(theseFiles); continue; end
    
    % so here, load what you care about and play with them!
    fprintf('loading...')

    tfrManips = {...
        getVar('ec_responseLockedTFR%s') getVar('hc_responseLockedTFR%s')...
        getVar('er_responseLockedTFR%s') getVar('hr_responseLockedTFR%s')...
        getVar('ec_coherenceLockedTFR%s') getVar('hc_coherenceLockedTFR%s')...
        getVar('er_coherenceLockedTFR%s') getVar('hr_coherenceLockedTFR%s')...
    };
    tfrConds = {...
        getVar('ecer_coherenceLockedTFR%s')...
        getVar('echr_responseLockedTFR%s') getVar('echr_coherenceLockedTFR%s')...
        getVar('hcer_responseLockedTFR%s') getVar('hcer_coherenceLockedTFR%s')...
        getVar('hchr_responseLockedTFR%s') getVar('hchr_coherenceLockedTFR%s')...
    };

    whichVars = {...
            tfrManips{:} tfrConds{:}...
        };
    
    data{subjectNum} = load(thisFile,whichVars{:});
    
    fprintf('loaded\n')
    
end; clear theseFiles thisFile subjectFolders subjectNum subjectCode index pathParts
clear tfrManips tfrConds
response = input('Do you want to delete the parallel pool? (y/n): ','s');
if strcmpi(response,'y') || strcmpi(response,'yes')
    disp('deleting...');
    delete(gcp('nocreate')); clear workers;
    system(['rm -rf ' tempPath]);
    if exist(tempPath,'file'); warning('I couldnt delete the job directory :('); end
else
    disp('not deleting...')
end

fprintf('loading complete with %.0f subjects\n',sum(~cellfun(@isempty,data)))

% get some layouts and calculate neighbours

disp('prep layouts and neighbours')

% meg is standard---can just use FTs version
megLayout = fullfile(ftDir,'template','layout','neuromag306all.lay');
tmp = load(fullfile(ftDir,'template','neighbours','neuromag306mag_neighb.mat'),'neighbours');
megNeighbours = tmp.neighbours;
% or make our own from a subject:
% megLayout = fullfile(datadir,'S01','Preprocess','meg_layout.lay'); % note this layout will squeeze everything into an eeg circle in a topo plot---hard to read
% cfg = [];
% cfg.channel = 'MEG';
% cfg.method = 'distance';
% cfg.layout = megLayout;
% megNeighbours = ft_prepare_neighbours(cfg);

% for EEG, we have a modified 70channel cap to fit 64 channels, so although
% there is a template for 64 chan eeg e.g.:
% eegLayout = fullfile(ftDir,'template','layout','acticap-64ch-standard2.mat');
% eegNeighbours = fullfile(ftDir,'template','neighbours','easycap64ch-avg_neighb.mat');
% I think better to make this ourselves from a participant layout
eegLayout = fullfile(datadir,'S01','Preprocess','eeg_layout.lay');

% now neighbours

cfg = [];
cfg.channel = 'EEG';
cfg.method = 'distance';
cfg.layout = eegLayout;
eegNeighbours = ft_prepare_neighbours(cfg);


% and we'll collect useful sensors

% we can plot these
cfg = [];
cfg.layout = eegLayout;
% % cfg.layout = megLayout;
layout = ft_prepare_layout(cfg);
% ft_plot_layout(layout);
clear layout




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% let's compile some averages and differences %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------%
%-- average coherence data --%
%----------------------------%

disp('averaging coherence data')

cfg = [];
ecRespTfrAll = returnStructs(data, getVar('ec_responseLockedTFR%s'));
ecRespTfrAve = ft_freqgrandaverage(cfg, ecRespTfrAll{:});
cfg = [];
ecOnsTfrAll = returnStructs(data, getVar('ec_coherenceLockedTFR%s'));
ecOnsTfrAve = ft_freqgrandaverage(cfg, ecOnsTfrAll{:});
cfg = [];
hcRespTfrAll = returnStructs(data, getVar('hc_responseLockedTFR%s'));
hcRespTfrAve = ft_freqgrandaverage(cfg, hcRespTfrAll{:});
cfg = [];
hcOnsTfrAll = returnStructs(data, getVar('hc_coherenceLockedTFR%s'));
hcOnsTfrAve = ft_freqgrandaverage(cfg, hcOnsTfrAll{:});

% get an overall average
cfg = [];
cOnsTfrAve = ft_freqgrandaverage(cfg, ecOnsTfrAll{:}, hcOnsTfrAll{:});
cfg = [];
cRespTfrAve = ft_freqgrandaverage(cfg, ecRespTfrAll{:}, hcRespTfrAll{:});

% get diffs
cOnsTfrDiffAve = getTfrDiff(ecOnsTfrAve, hcOnsTfrAve);
cRespTfrDiffAve = getTfrDiff(ecRespTfrAve, hcRespTfrAve);

% % make a dummy structure with the difference of the timefrequencies
% cOnsTfrDiffAve = ecOnsTfrAve;
% cOnsTfrDiffAve.powspctrm = ecOnsTfrAve.powspctrm - hcOnsTfrAve.powspctrm;
% cRespTfrDiffAve = ecRespTfrAve;
% cRespTfrDiffAve.powspctrm = ecRespTfrAve.powspctrm - hcRespTfrAve.powspctrm;

disp('done')

%---------------------------------%
%-- average categorisation data --%
%---------------------------------%

disp('averaging categorisation data')

cfg = [];
erRespTfrAll = returnStructs(data, getVar('er_responseLockedTFR%s'));
erRespTfrAve = ft_freqgrandaverage(cfg, erRespTfrAll{:});
cfg = [];
erOnsTfrAll = returnStructs(data, getVar('er_coherenceLockedTFR%s'));
erOnsTfrAve = ft_freqgrandaverage(cfg, erOnsTfrAll{:});
cfg = [];
hrRespTfrAll = returnStructs(data, getVar('hr_responseLockedTFR%s'));
hrRespTfrAve = ft_freqgrandaverage(cfg, hrRespTfrAll{:});
cfg = [];
hrOnsTfrAll = returnStructs(data, getVar('hr_coherenceLockedTFR%s'));
hrOnsTfrAve = ft_freqgrandaverage(cfg, hrOnsTfrAll{:});

% get an overall average
cfg = [];
rOnsTfrAve = ft_freqgrandaverage(cfg, erOnsTfrAll{:}, hrOnsTfrAll{:});
cfg = [];
rRespTfrAve = ft_freqgrandaverage(cfg, erRespTfrAll{:}, hrRespTfrAll{:});

% get diffs
rOnsTfrDiffAve = getTfrDiff(erOnsTfrAve, hrOnsTfrAve);
rRespTfrDiffAve = getTfrDiff(erRespTfrAve, hrRespTfrAve);

% % make a dummy structure with the difference of the timefrequencies
% rOnsTfrDiffAve = erOnsTfrAve;
% rOnsTfrDiffAve.powspctrm = erOnsTfrAve.powspctrm - hrOnsTfrAve.powspctrm;
% rRespTfrDiffAve = erRespTfrAve;
% rRespTfrDiffAve.powspctrm = erRespTfrAve.powspctrm - hrRespTfrAve.powspctrm;

disp('done')

%--------------------------------%
%-- average conditionwise data --%
%--------------------------------%

disp('averaging conditionwise data')

cfg = [];
ecerRespTfrAll = returnStructs(data, getVar('ecer_responseLockedTFR%s'));
ecerRespTfrAve = ft_freqgrandaverage(cfg, ecerRespTfrAll{:});
cfg = [];
ecerOnsTfrAll = returnStructs(data, getVar('ecer_coherenceLockedTFR%s'));
ecerOnsTfrAve = ft_freqgrandaverage(cfg, ecerOnsTfrAll{:});
cfg = [];
echrRespTfrAll = returnStructs(data, getVar('echr_responseLockedTFR%s'));
echrRespTfrAve = ft_freqgrandaverage(cfg, echrRespTfrAll{:});
cfg = [];
echrOnsTfrAll = returnStructs(data, getVar('echr_coherenceLockedTFR%s'));
echrOnsTfrAve = ft_freqgrandaverage(cfg, echrOnsTfrAll{:});
cfg = [];
hcerRespTfrAll = returnStructs(data, getVar('hcer_responseLockedTFR%s'));
hcerRespTfrAve = ft_freqgrandaverage(cfg, hcerRespTfrAll{:});
cfg = [];
hcerOnsTfrAll = returnStructs(data, getVar('hcer_coherenceLockedTFR%s'));
hcerOnsTfrAve = ft_freqgrandaverage(cfg, hcerOnsTfrAll{:});
cfg = [];
hchrRespTfrAll = returnStructs(data, getVar('hchr_responseLockedTFR%s'));
hchrRespTfrAve = ft_freqgrandaverage(cfg, hchrRespTfrAll{:});
cfg = [];
hchrOnsTfrAll = returnStructs(data, getVar('hchr_coherenceLockedTFR%s'));
hchrOnsTfrAve = ft_freqgrandaverage(cfg, hchrOnsTfrAll{:});

disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get subjectwise differences %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('getting differences')

for subject = 1:numel(ecOnsTfrAll)
    fprintf('\ngetting subjectwise differences for subject %.0f of %.0f\n\n',subject,numel(ecOnsTfrAll))
    
    cOnsTfrDiffAll{subject} = getTfrDiff(ecOnsTfrAll{subject}, hcOnsTfrAll{subject});
    cRespTfrDiffAll{subject} = getTfrDiff(ecRespTfrAll{subject}, hcRespTfrAll{subject});
    rOnsTfrDiffAll{subject} = getTfrDiff(erOnsTfrAll{subject}, hrOnsTfrAll{subject});
    rRespTfrDiffAll{subject} = getTfrDiff(erRespTfrAll{subject}, hrRespTfrAll{subject});

%     % make a dummy structure with the difference of the timefrequencies
%     cOnsTfrDiffAll{subject} = ecOnsTfrAll{subject};
%     cOnsTfrDiffAll{subject}.powspctrm = ecOnsTfrAll{subject}.powspctrm - hcOnsTfrAll{subject}.powspctrm;
%     cRespTfrDiffAll{subject} = ecRespTfrAll{subject};
%     cRespTfrDiffAll{subject}.powspctrm = ecRespTfrAll{subject}.powspctrm - hcRespTfrAll{subject}.powspctrm;
%     rOnsTfrDiffAll{subject} = erOnsTfrAll{subject};
%     rOnsTfrDiffAll{subject}.powspctrm = erOnsTfrAll{subject}.powspctrm - hrOnsTfrAll{subject}.powspctrm;
%     rRespTfrDiffAll{subject} = erRespTfrAll{subject};
%     rRespTfrDiffAll{subject}.powspctrm = erRespTfrAll{subject}.powspctrm - hrRespTfrAll{subject}.powspctrm;

end; clear subject


disp('done')


%% save those differences, if we want

disp('saving')

save(fullfile(datadir,getVar('tfr_differences_%s.mat')),...
    'cOnsTfrDiffAll', 'cRespTfrDiffAll', 'rOnsTfrDiffAll', 'rRespTfrDiffAll', 'cOnsTfrDiffAve', 'cRespTfrDiffAve', 'rOnsTfrDiffAve', 'rRespTfrDiffAve'...
    )

disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the timecourse of univariate activity %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%% coh onset
cfg = [];
% cfg.baseline     = [-0.5 -0.2]; % if we haven't normalised by difference
% cfg.baselinetype = 'absolute';
% cfg.zlim         = [-0.23 0.21];
cfg.showlabels   = 'yes';
cfg.layout       = megLayout;
cfg.maskstyle = 'saturation';
% cfg.xlim=[-0.5 1.5];
% whichChannels = 'parietal'; cfg.channel = getMegLabels(whichChannels);

if isfield(cfg,'channel')
    savename = [tfrFigDir filesep whichChannels '_' getVar('easy_coherence_onsetLockedTFR%s') '_%s.png'];
else
    savename = [tfrFigDir filesep getVar('easy_coherence_onsetLockedTFR%s') '_%s.png'];
end
ecOnsTfrAve.bfs = testDiffsAcrossTime(ecOnsTfrAll,cfg);
plotTfr(ecOnsTfrAve, cfg, ecOnsTfrAve.bfs, sprintf(savename, 'tfr'))
if contains(savename,'hann'); freqs = '(1-30Hz)'; elseif contains(savename,'multi'); freqs = '(30-150Hz)'; end
suptitle(['TFR Difference in Onset Locked Coherence Easy Difficulty ' freqs])

if isfield(cfg,'channel')
    savename = [tfrFigDir filesep whichChannels '_' getVar('hard_coherence_onsetLockedTFR%s') '_%s.png'];
else
    savename = [tfrFigDir filesep getVar('hard_coherence_onsetLockedTFR%s') '_%s.png'];
end
hcOnsTfrAve.bfs = testDiffsAcrossTime(hcOnsTfrAll,cfg);
plotTfr(hcOnsTfrAve, cfg, hcOnsTfrAve.bfs, sprintf(savename, 'tfr'))
if contains(savename,'hann'); freqs = '(1-30Hz)'; elseif contains(savename,'multi'); freqs = '(30-150Hz)'; end
suptitle(['TFR Difference in Onset Locked Coherence Hard Difficulty ' freqs])

%% coh onset diff
cfg = [];
% cfg.baseline     = [-0.5 -0.2]; % if we haven't normalised by difference
% cfg.baselinetype = 'absolute';
% cfg.zlim         = [-0.23 0.21];
cfg.showlabels   = 'yes';
cfg.layout       = megLayout;
cfg.maskstyle = 'saturation';
% cfg.xlim=[-0.5 1.5];
% whichChannels = 'parietal'; cfg.channel = getMegLabels(whichChannels);

if isfield(cfg,'channel')
    savename = [tfrFigDir filesep whichChannels '_' getVar('coherence_onsetLockedTFR%s') '_%s.png'];
else
    savename = [tfrFigDir filesep getVar('coherence_onsetLockedTFR%s') '_%s.png'];
end
cOnsTfrDiffAve.bfs = testDiffsAcrossTime(cOnsTfrDiffAll,cfg);
plotTfr(cOnsTfrDiffAve, cfg, cOnsTfrDiffAve.bfs, sprintf(savename, 'tfr'))
if contains(savename,'hann'); freqs = '(1-30Hz)'; elseif contains(savename,'multi'); freqs = '(30-150Hz)'; end
suptitle(['TFR Difference in Onset Locked Coherence Difficulty ' freqs])

% ft_multiplotTFR(cfg, cOnsTfrDiffAve);

if contains(savename,'hann')
    
    % delta positivity
    cfg.ylim = [1 5]; % freqs for topo
    cfg.xlim = [0.25 1.05];
    ft_topoplotTFR(cfg, cOnsTfrDiffAve);
    title(['TFR Difference in Onset Locked Coherence Difficulty '...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))])
    print(sprintf(savename, ['topo_'...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))]), '-dpng');
    
    % beta positivity
    cfg.ylim = [15 30]; % freqs for topo
    cfg.xlim = [1.0 1.5];
    ft_topoplotTFR(cfg, cOnsTfrDiffAve);
    title(['TFR Difference in Onset Locked Coherence Difficulty '...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))])
    print(sprintf(savename, ['topo_'...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))]), '-dpng');
    
    % beta negativity
    cfg.ylim = [10 30]; % freqs for topo
    cfg.xlim = [0.25 0.65];
    ft_topoplotTFR(cfg, cOnsTfrDiffAve);
    title(['TFR Difference in Onset Locked Coherence Difficulty '...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))])
    print(sprintf(savename, ['topo_'...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))]), '-dpng');
    
elseif contains(savename,'multi')
    
    % low gamma negativity
    cfg.ylim = [30 38]; % freqs for topo
    cfg.xlim = [0.3 0.65];
    ft_topoplotTFR(cfg, cOnsTfrDiffAve);
    title(['TFR Difference in Onset Locked Coherence Difficulty '...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))])
    print(sprintf(savename, ['topo_'...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))]), '-dpng');
    
    % low gamma positivity
    cfg.ylim = [30 41]; % freqs for topo
    cfg.xlim = [0.9 1.5];
    ft_topoplotTFR(cfg, cOnsTfrDiffAve);
    title(['TFR Difference in Onset Locked Coherence Difficulty '...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))])
    print(sprintf(savename, ['topo_'...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))]), '-dpng');
end

%% coh response
cfg = [];
% cfg.baseline     = [-0.5 -0.2]; % if we haven't normalised by difference
% cfg.baselinetype = 'absolute';
% cfg.zlim         = [-0.23 0.21];
cfg.showlabels   = 'yes';
cfg.layout       = megLayout;
cfg.maskstyle = 'saturation';
% cfg.xlim=[-0.5 1.5];
% whichChannels = 'parietal'; cfg.channel = getMegLabels(whichChannels);

if isfield(cfg,'channel')
    savename = [tfrFigDir filesep whichChannels '_' getVar('easy_coherence_responseLockedTFR%s') '_%s.png'];
else
    savename = [tfrFigDir filesep getVar('easy_coherence_responseLockedTFR%s') '_%s.png'];
end
ecRespTfrAve.bfs = testDiffsAcrossTime(ecRespTfrAll,cfg);
plotTfr(ecRespTfrAve, cfg, ecRespTfrAve.bfs, sprintf(savename, 'tfr'))
if contains(savename,'hann'); freqs = '(1-30Hz)'; elseif contains(savename,'multi'); freqs = '(30-150Hz)'; end
suptitle(['TFR Difference in Response Locked Coherence Easy Difficulty ' freqs])

if isfield(cfg,'channel')
    savename = [tfrFigDir filesep whichChannels '_' getVar('hard_coherence_responseLockedTFR%s') '_%s.png'];
else
    savename = [tfrFigDir filesep getVar('hard_coherence_responseLockedTFR%s') '_%s.png'];
end
hcRespTfrAve.bfs = testDiffsAcrossTime(hcRespTfrAll,cfg);
plotTfr(hcRespTfrAve, cfg, hcRespTfrAve.bfs, sprintf(savename, 'tfr'))
if contains(savename,'hann'); freqs = '(1-30Hz)'; elseif contains(savename,'multi'); freqs = '(30-150Hz)'; end
suptitle(['TFR Difference in Response Locked Coherence Hard Difficulty ' freqs])

%% coh response diff
cfg = [];
% cfg.baseline     = [-1.0 -0.5];
% cfg.baselinetype = 'absolute';
cfg.showlabels   = 'yes';
cfg.layout       = megLayout;
cfg.maskstyle = 'saturation';
% cfg.zlim         = [-0.12 0.22];
% cfg.xlim=[-0.5 1.5];
% whichChannels = 'parietal'; cfg.channel = getMegLabels(whichChannels);

if isfield(cfg,'channel')
    savename = [tfrFigDir filesep whichChannels '_' getVar('coherence_responseLockedTFR%s') '_%s.png'];
else
    savename = [tfrFigDir filesep getVar('coherence_responseLockedTFR%s') '_%s.png'];
end
cRespTfrDiffAve.bfs = testDiffsAcrossTime(cRespTfrDiffAll,cfg);
plotTfr(cRespTfrDiffAve, cfg, cRespTfrDiffAve.bfs, sprintf(savename, 'tfr'))
if contains(savename,'hann'); freqs = '(1-30Hz)'; elseif contains(savename,'multi'); freqs = '(30-150Hz)'; end
suptitle(['TFR Difference in Response Locked Coherence Difficulty ' freqs])

% ft_multiplotTFR(cfg, cRespTfrDiffAve);

if contains(savename,'hann')
    
    % beta negativity to zero
    cfg.ylim = [10 21]; % freqs for topo
    cfg.xlim = [-0.25 0];
    ft_topoplotTFR(cfg, cRespTfrDiffAve);
    title(['TFR Difference in Response Locked Coherence Difficulty '...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))])
    print(sprintf(savename, ['topo_'...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))]), '-dpng');
    
    % delta positivity
    cfg.ylim = [1 5]; % freqs for topo
    cfg.xlim = [-0.6 0.2];
    ft_topoplotTFR(cfg, cRespTfrDiffAve);
    title(['TFR Difference in Response Locked Coherence Difficulty '...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))])
    print(sprintf(savename, ['topo_'...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))]), '-dpng');
    
elseif contains(savename,'multi')
    
    % low gamma negativity
    cfg.ylim = [30 38]; % freqs for topo
    cfg.xlim = [-0.3 -0.05];
    ft_topoplotTFR(cfg, cOnsTfrDiffAve);
    title(['TFR Difference in Onset Locked Coherence Difficulty '...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))])
    print(sprintf(savename, ['topo_'...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))]), '-dpng');
    
    % low gamma positivity
    cfg.ylim = [30 38]; % freqs for topo
    cfg.xlim = [0.15 0.2];
    ft_topoplotTFR(cfg, cOnsTfrDiffAve);
    title(['TFR Difference in Onset Locked Coherence Difficulty '...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))])
    print(sprintf(savename, ['topo_'...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))]), '-dpng');
end

%% cat onset
cfg = [];
% cfg.baseline     = [-0.5 -0.2]; % if we haven't normalised by difference
% cfg.baselinetype = 'absolute';
% cfg.zlim         = [-0.23 0.21];
cfg.showlabels   = 'yes';
cfg.layout       = megLayout;
cfg.maskstyle = 'saturation';
% cfg.xlim=[-0.5 1.5];
% whichChannels = 'parietal'; cfg.channel = getMegLabels(whichChannels);

if isfield(cfg,'channel')
    savename = [tfrFigDir filesep whichChannels '_' getVar('easy_categorisation_onsetLockedTFR%s') '_%s.png'];
else
    savename = [tfrFigDir filesep getVar('easy_categorisation_onsetLockedTFR%s') '_%s.png'];
end
% erOnsTfrAve.bfs = testDiffsAcrossTime(erOnsTfrAll,cfg);
plotTfr(erOnsTfrAve, cfg, erOnsTfrAve.bfs, sprintf(savename, 'tfr'))
if contains(savename,'hann'); freqs = '(1-30Hz)'; elseif contains(savename,'multi'); freqs = '(30-150Hz)'; end
suptitle(['TFR Difference in Onset Locked Categorisation Easy Difficulty ' freqs])

if isfield(cfg,'channel')
    savename = [tfrFigDir filesep whichChannels '_' getVar('hard_categorisation_onsetLockedTFR%s') '_%s.png'];
else
    savename = [tfrFigDir filesep getVar('hard_categorisation_onsetLockedTFR%s') '_%s.png'];
end
hrOnsTfrAve.bfs = testDiffsAcrossTime(hrOnsTfrAll,cfg);
plotTfr(hrOnsTfrAve, cfg, hrOnsTfrAve.bfs, sprintf(savename, 'tfr'))
if contains(savename,'hann'); freqs = '(1-30Hz)'; elseif contains(savename,'multi'); freqs = '(30-150Hz)'; end
suptitle(['TFR Difference in Onset Locked Categorisation Hard Difficulty ' freqs])

%% cat onset diff
cfg = [];
% cfg.baseline     = [-0.5 -0.2];
% cfg.baselinetype = 'absolute';
% cfg.zlim         = [-0.23 0.21];
cfg.showlabels   = 'yes';
cfg.layout       = megLayout;
cfg.maskstyle = 'saturation';
% cfg.xlim=[-0.5 1.5];
% whichChannels = 'parietal'; cfg.channel = getMegLabels(whichChannels);

if isfield(cfg,'channel')
    savename = [tfrFigDir filesep whichChannels '_' getVar('categorisation_onsetLockedTFR%s') '_%s.png'];
else
    savename = [tfrFigDir filesep getVar('categorisation_onsetLockedTFR%s') '_%s.png'];
end
rOnsTfrDiffAve.bfs = testDiffsAcrossTime(rOnsTfrDiffAll,cfg);
plotTfr(rOnsTfrDiffAve, cfg, rOnsTfrDiffAve.bfs, sprintf(savename, 'tfr'))
if contains(savename,'hann'); freqs = '(1-30Hz)'; elseif contains(savename,'multi'); freqs = '(30-150Hz)'; end
suptitle(['TFR Difference in Onset Locked Categorisation Difficulty ' freqs])

% ft_multiplotTFR(cfg, rOnsTfrDiffAve);

if contains(savename,'hann')
    
    % beta positivity
    cfg.ylim = [15 21]; % freqs for topo
    cfg.xlim = [0.55 1.15];
    ft_topoplotTFR(cfg, cRespTfrDiffAve);
    title(['TFR Difference in Onset Locked Categorisation Difficulty '...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))])
    print(sprintf(savename, ['topo_'...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))]), '-dpng');
    
elseif contains(savename,'multi')
end

%% cat response
cfg = [];
% cfg.baseline     = [-0.5 -0.2]; % if we haven't normalised by difference
% cfg.baselinetype = 'absolute';
% cfg.zlim         = [-0.23 0.21];
cfg.showlabels   = 'yes';
cfg.layout       = megLayout;
cfg.maskstyle = 'saturation';
% cfg.xlim=[-0.5 1.5];
% whichChannels = 'parietal'; cfg.channel = getMegLabels(whichChannels);

if isfield(cfg,'channel')
    savename = [tfrFigDir filesep whichChannels '_' getVar('easy_categorisation_responseLockedTFR%s') '_%s.png'];
else
    savename = [tfrFigDir filesep getVar('easy_categorisation_responseLockedTFR%s') '_%s.png'];
end
erRespTfrAve.bfs = testDiffsAcrossTime(erRespTfrAll,cfg);
plotTfr(erRespTfrAve, cfg, erRespTfrAve.bfs, sprintf(savename, 'tfr'))
if contains(savename,'hann'); freqs = '(1-30Hz)'; elseif contains(savename,'multi'); freqs = '(30-150Hz)'; end
suptitle(['TFR Difference in Response Locked Categorisation Easy Difficulty ' freqs])

if isfield(cfg,'channel')
    savename = [tfrFigDir filesep whichChannels '_' getVar('hard_categorisation_responseLockedTFR%s') '_%s.png'];
else
    savename = [tfrFigDir filesep getVar('hard_categorisation_responseLockedTFR%s') '_%s.png'];
end
hrRespTfrAve.bfs = testDiffsAcrossTime(hrRespTfrAll,cfg);
plotTfr(hrRespTfrAve, cfg, hrRespTfrAve.bfs, sprintf(savename, 'tfr'))
if contains(savename,'hann'); freqs = '(1-30Hz)'; elseif contains(savename,'multi'); freqs = '(30-150Hz)'; end
suptitle(['TFR Difference in Response Locked Categorisation Hard Difficulty ' freqs])

%% cat response diff
cfg = [];
% cfg.baseline     = [-1.0 -0.5];
% cfg.baselinetype = 'absolute';
cfg.showlabels   = 'yes';
cfg.layout       = megLayout;
cfg.maskstyle = 'saturation';
% cfg.zlim         = [-0.12 0.22];
% cfg.xlim=[-0.5 1.5];
% whichChannels = 'parietal'; cfg.channel = getMegLabels(whichChannels);

if isfield(cfg,'channel')
    savename = [tfrFigDir filesep whichChannels '_' getVar('categorisation_responseLockedTFR%s') '_%s.png'];
else
    savename = [tfrFigDir filesep getVar('categorisation_responseLockedTFR%s') '_%s.png'];
end
rRespTfrDiffAve.bfs = testDiffsAcrossTime(rRespTfrDiffAll,cfg);
plotTfr(rRespTfrDiffAve, cfg, rRespTfrDiffAve.bfs, sprintf(savename, 'tfr'))
if contains(savename,'hann'); freqs = '(1-30Hz)'; elseif contains(savename,'multi'); freqs = '(30-150Hz)'; end
suptitle(['TFR Difference in Response Locked Categorisation Difficulty ' freqs])

% ft_multiplotTFR(cfg, rRespTfrDiffAve);

if contains(savename,'hann')
    
    cfg.ylim = [18 22]; % freqs for topo
    cfg.xlim = [-0.6 -0.25];
    ft_topoplotTFR(cfg, cRespTfrDiffAve);
    title(['TFR Difference in Onset Locked Categorisation Difficulty '...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))])
    print(sprintf(savename, ['topo_'...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))]), '-dpng');
    
elseif contains(savename,'multi')
    
        cfg.ylim = [30 32]; % freqs for topo
    cfg.xlim = [-0.6 -0.20];
    ft_topoplotTFR(cfg, cRespTfrDiffAve);
    title(['TFR Difference in Onset Locked Categorisation Difficulty '...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))])
    print(sprintf(savename, ['topo_'...
        num2str(cfg.ylim(1)) '-' num2str(cfg.ylim(2)) '@'...
        num2str(cfg.xlim(1)) '-' num2str(cfg.xlim(2))]), '-dpng');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the topology of univariate activity %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% where is easy coh different from hard coh

% let's get some zlimits
        

%% subfunctions

function diffStruct = getTfrDiff(struct1, struct2)
% with two conditions of similar baseline, we can compare by subtracting 
% and then dividing the two powerspectra by their sum: normalising common acty

cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = '(x1-x2)/(x1+x2)';
diffStruct = ft_math(cfg, struct1, struct2);

return
end

function bfs = testDiffsAcrossTime(dataStruct,cfg,complementary)

numSubjs = numel(dataStruct);
numFreqs = numel(dataStruct{1}.freq);
numTimes = numel(dataStruct{1}.time);
for subject = 1:numSubjs
    for frequency = 1:numFreqs
        for timepoint = 1:numTimes

            if isfield(cfg,'channel')
                chanIdx = find(ismember(dataStruct{subject}.label,cfg.channel));
            else
                chanIdx = (1:numel(dataStruct{subject}.label))';
            end
            % rows = tests
            % cols = observations (subjects)
            diffs(timepoint,frequency,subject) = mean(dataStruct{subject}.powspctrm(chanIdx,frequency,timepoint));
            
        end
    end; clear frequency
end; clear subj
reshapedDiffs = reshape(diffs, numTimes*numFreqs, numSubjs);
    

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
bfs = bayesfactor_R_wrapper(reshapedDiffs,'Rpath',RscriptPath,'returnindex',complementary,...
    'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');

bfs = reshape(bfs, numTimes, numFreqs); % comes out rotated
bfs = flip(bfs.'); % so rotate it back


% alternatively they have implemented it in matlab
% [bfs, bfs_complementary_interval] = bayesfactor(diffs, 'interval',[-Inf Inf]);

return
end

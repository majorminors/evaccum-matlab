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

erpFigDir = fullfile(datadir, 'erpFigs');
if ~exist(erpFigDir,'dir'); mkdir(erpFigDir); end

% we'll loop through subject data dirs, so just get the path from there
inputFileName = ['Preprocess' filesep 'tfr_hanning.mat'];

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
    
    disp('this is subject:')
    disp(subjectFolders(subjectNum).name);
    
     % grab info about all the meeg data files we care about
    theseFiles = dir([subjectFolders(subjectNum).folder filesep subjectFolders(subjectNum).name filesep inputFileName]);
    if isempty(theseFiles); continue; end
    
    % so here, load what you care about and play with them!
    disp('loading')

    tfrManips = {...
        'ec_responseLockedTFRhann' 'hc_responseLockedTFRhann'...
        'er_responseLockedTFRhann' 'hr_responseLockedTFRhann'...
        'ec_coherenceLockedTFRhann' 'hc_coherenceLockedTFRhann'...
        'er_coherenceLockedTFRhann' 'hr_coherenceLockedTFRhann'...
    };
    tfrConds = {...
        'ecer_responseLockedTFRhann' 'ecer_coherenceLockedTFRhann'...
        'echr_responseLockedTFRhann' 'echr_coherenceLockedTFRhann'...
        'hcer_responseLockedTFRhann' 'hcer_coherenceLockedTFRhann'...
        'hchr_responseLockedTFRhann' 'hchr_coherenceLockedTFRhann'...
    };

    whichVars = {...
            tfrManips{:} tfrConds{:}...
        };
    
    data{subjectNum} = load(thisFile,whichVars{:});
    
    disp('loaded')
    
end; clear theseFiles thisFile subjectFolders subjectNum subjectCode index pathParts
clear tfrManips tfrConds
delete(gcp('nocreate')); clear workers;
rmdir(tempPath,'s');

disp('loading complete')

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
% clear layout


CPP = {'EEG040' 'EEG041' 'EEG042'};
frontal = {'EEG004' 'EEG002' 'EEG008'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% let's compile some averages and differences %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------%
%-- average coherence data --%
%----------------------------%

disp('averaging coherence data')

cfg = [];
ecRespTfrAll = returnStructs(data, 'ec_responseLockedTFRhann');
ecRespTfrAve = ft_freqgrandaverage(cfg, ecRespTfrAll{:});
cfg = [];
ecOnsTfrAll = returnStructs(data, 'ec_coherenceLockedTFRhann');
ecOnsTfrAve = ft_freqgrandaverage(cfg, ecOnsTfrAll{:});
cfg = [];
hcRespTfrAll = returnStructs(data, 'hc_responseLockedTFRhann');
hcRespTfrAve = ft_freqgrandaverage(cfg, hcRespTfrAll{:});
cfg = [];
hcOnsTfrAll = returnStructs(data, 'hc_coherenceLockedTFRhann');
hcOnsTfrAve = ft_freqgrandaverage(cfg, hcOnsTfrAll{:});

% get an overall average
cfg = [];
cOnsTfrAve = ft_freqgrandaverage(cfg, ecOnsTfrAll{:}, hcOnsTfrAll{:});
cfg = [];
cRespTfrAve = ft_freqgrandaverage(cfg, ecRespTfrAll{:}, hcRespTfrAll{:});

% get diffs
% make a dummy structure with the difference of the timefrequencies
cOnsTfrDiffAve = ecOnsTfrAve;
cOnsTfrDiffAve.powspctrm = ecOnsTfrAve.powspctrm - hcOnsTfrAve.powspctrm;
cRespTfrDiffAve = ecRespTfrAve;
cRespTfrDiffAve.powspctrm = ecRespTfrAve.powspctrm - hcRespTfrAve.powspctrm;

disp('done')

%---------------------------------%
%-- average categorisation data --%
%---------------------------------%

disp('averaging categorisation data')

cfg = [];
erRespTfrAll = returnStructs(data, 'er_responseLockedTFRhann');
erRespTfrAve = ft_freqgrandaverage(cfg, erRespTfrAll{:});
cfg = [];
erOnsTfrAll = returnStructs(data, 'er_coherenceLockedTFRhann');
erOnsTfrAve = ft_freqgrandaverage(cfg, erOnsTfrAll{:});
cfg = [];
hrRespTfrAll = returnStructs(data, 'hr_responseLockedTFRhann');
hrRespTfrAve = ft_freqgrandaverage(cfg, hrRespTfrAll{:});
cfg = [];
hrOnsTfrAll = returnStructs(data, 'hr_coherenceLockedTFRhann');
hrOnsTfrAve = ft_freqgrandaverage(cfg, hrOnsTfrAll{:});

% get an overall average
cfg = [];
rOnsTfrAve = ft_freqgrandaverage(cfg, erOnsTfrAll{:}, hrOnsTfrAll{:});
cfg = [];
rRespTfrAve = ft_freqgrandaverage(cfg, erRespTfrAll{:}, hrRespTfrAll{:});

% get diffs
% make a dummy structure with the difference of the timefrequencies
rOnsTfrDiffAve = erOnsTfrAve;
rOnsTfrDiffAve.powspctrm = erOnsTfrAve.powspctrm - hrOnsTfrAve.powspctrm;
rRespTfrDiffAve = erRespTfrAve;
rRespTfrDiffAve.powspctrm = erRespTfrAve.powspctrm - hrRespTfrAve.powspctrm;

disp('done')

%--------------------------------%
%-- average conditionwise data --%
%--------------------------------%

disp('averaging conditionwise data')

cfg = [];
ecerRespTfrAll = returnStructs(data, 'ecer_responseLockedTFRhann');
ecerRespTfrAve = ft_freqgrandaverage(cfg, ecerRespTfrAll{:});
cfg = [];
ecerOnsTfrAll = returnStructs(data, 'ecer_coherenceLockedTFRhann');
ecerOnsTfrAve = ft_freqgrandaverage(cfg, ecerOnsTfrAll{:});
cfg = [];
echrRespTfrAll = returnStructs(data, 'echr_responseLockedTFRhann');
echrRespTfrAve = ft_freqgrandaverage(cfg, echrRespTfrAll{:});
cfg = [];
echrOnsTfrAll = returnStructs(data, 'echr_coherenceLockedTFRhann');
echrOnsTfrAve = ft_freqgrandaverage(cfg, echrOnsTfrAll{:});
cfg = [];
hcerRespTfrAll = returnStructs(data, 'hcer_responseLockedTFRhann');
hcerRespTfrAve = ft_freqgrandaverage(cfg, hcerRespTfrAll{:});
cfg = [];
hcerOnsTfrAll = returnStructs(data, 'hcer_coherenceLockedTFRhann');
hcerOnsTfrAve = ft_freqgrandaverage(cfg, hcerOnsTfrAll{:});
cfg = [];
hchrRespTfrAll = returnStructs(data, 'hchr_responseLockedTFRhann');
hchrRespTfrAve = ft_freqgrandaverage(cfg, hchrRespTfrAll{:});
cfg = [];
hchrOnsTfrAll = returnStructs(data, 'hchr_coherenceLockedTFRhann');
hchrOnsTfrAve = ft_freqgrandaverage(cfg, hchrOnsTfrAll{:});

disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get subjectwise differences %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('getting differences')

for subject = 1:numel(ecOnsTfrAll)
    fprintf('\ngetting subjectwise differences for subject %.0f of %.0f\n\n',subject,numel(ecOnsTfrAll))
    
    % make a dummy structure with the difference of the timefrequencies
    cOnsTfrDiffAll{subject} = ecOnsTfrAll{subject};
    cOnsTfrDiffAll{subject}.powspctrm = ecOnsTfrAll{subject}.powspctrm - hcOnsTfrAll{subject}.powspctrm;
    cRespTfrDiffAll{subject} = ecRespTfrAll{subject};
    cRespTfrDiffAll{subject}.powspctrm = ecRespTfrAll{subject}.powspctrm - hcRespTfrAll{subject}.powspctrm;
    rOnsTfrDiffAll{subject} = erOnsTfrAll{subject};
    rOnsTfrDiffAll{subject}.powspctrm = erOnsTfrAll{subject}.powspctrm - hrOnsTfrAll{subject}.powspctrm;
    rRespTfrDiffAll{subject} = erRespTfrAll{subject};
    rRespTfrDiffAll{subject}.powspctrm = erRespTfrAll{subject}.powspctrm - hrRespTfrAll{subject}.powspctrm;

end; clear subject


disp('done')

%% save those differences, if we want

disp('saving')

save(fullfile(datadir,'tfr_differences.mat'),...
    'cOnsTfrDiffAll', 'cRespTfrDiffAll', 'rOnsTfrDiffAll', 'rRespTfrDiffAll', 'cOnsTfrDiffAve', 'cRespTfrDiffAve', 'rOnsTfrDiffAve', 'rRespTfrDiffAve'...
    )

disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the timecourse of univariate activity %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%


cfg = [];
% cfg.baseline     = [-.5 0];
% cfg.baselinetype = 'absolute';
% cfg.zlim         = [-2.5e-24 2.5e-24];
cfg.showlabels   = 'yes';
cfg.layout       = megLayout;
cfg.xlim=[-0.5 1.5];
cfg.channel = getMegLabels('parietal');
ft_multiplotTFR(cfg, rRespTfrAve);
ft_singleplotTFR(cfg, rOnsTfrAve);
ft_topoplotTFR(cfg, rRespTfrDiffAve);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the topology of univariate activity %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% where is easy coh different from hard coh

% let's get some zlimits
        

%% subfunctions


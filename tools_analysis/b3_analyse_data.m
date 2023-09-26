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
inputFileName = ['Preprocess' filesep 'timelocked_averages.mat'];

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

    erpManips = {...
        'ec_responseLockedAverage' 'hc_responseLockedAverage'...
        'er_responseLockedAverage' 'hr_responseLockedAverage'...
        'ec_coherenceLockedAverage' 'hc_coherenceLockedAverage'...
        'er_coherenceLockedAverage' 'hr_coherenceLockedAverage'...
    };
    erpConds = {...
        'ecer_responseLockedAverage' 'ecer_coherenceLockedAverage'...
        'echr_responseLockedAverage' 'echr_coherenceLockedAverage'...
        'hcer_responseLockedAverage' 'hcer_coherenceLockedAverage'...
        'hchr_responseLockedAverage' 'hchr_coherenceLockedAverage'...
    };
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

    % erpManips{:} erpConds{:} tfrManips{:} tfrConds{:}

    whichVars = {...
            [erpManips erpConds tfrManips tfrConds]...
        };
    
    data{subjectNum} = load(thisFile,whichVars{:});
    
    disp('loaded')
    
end; clear theseFiles thisFile subjectFolders subjectNum subjectCode index pathParts
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
ft_plot_layout(layout);
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
ecRespAll = returnStructs(data, 'ec_responseLockedAverage');
ecRespAve = ft_timelockgrandaverage(cfg, ecRespAll{:});
ecRespTfrAll = returnStructs(data, 'ec_responseLockedTFRhann');
ecRespTfrAve = ft_freqgrandaverage(cfg, ecRespTfrAll{:});
cfg = [];
ecOnsAll = returnStructs(data, 'ec_coherenceLockedAverage');
ecOnsAve = ft_timelockgrandaverage(cfg, ecOnsAll{:});
ecOnsTfrAll = returnStructs(data, 'ec_coherenceLockedTFRhann');
ecOnsTfrAve = ft_freqgrandaverage(cfg, ecOnsTfrAll{:});
cfg = [];
hcRespAll = returnStructs(data, 'hc_responseLockedAverage');
hcRespAve = ft_timelockgrandaverage(cfg, hcRespAll{:});
hcRespTfrAll = returnStructs(data, 'hc_responseLockedTFRhann');
hcRespTfrAve = ft_freqgrandaverage(cfg, hcRespTfrAll{:});
cfg = [];
hcOnsAll = returnStructs(data, 'hc_coherenceLockedAverage');
hcOnsAve = ft_timelockgrandaverage(cfg, hcOnsAll{:});
hcOnsTfrAll = returnStructs(data, 'hc_coherenceLockedTFRhann');
hcOnsTfrAve = ft_freqgrandaverage(cfg, hcOnsTfrAll{:});

% get an overall average
cfg = [];
cOnsAve = ft_timelockgrandaverage(cfg, ecOnsAll{:}, hcOnsAll{:});
cOnsTfrAve = ft_freqgrandaverage(cfg, ecOnsTfrAll{:}, hcOnsTfrAll{:});
cfg = [];
cRespAve = ft_timelockgrandaverage(cfg, ecRespAll{:}, hcRespAll{:});
cRespTfrAve = ft_freqgrandaverage(cfg, ecRespTfrAll{:}, hcRespTfrAll{:});

% get diffs
cOnsDiffAve = getDifference(ecOnsAve,hcOnsAve);
cRespDiffAve = getDifference(ecRespAve,hcRespAve);
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
erRespAll = returnStructs(data, 'er_responseLockedAverage');
erRespAve = ft_timelockgrandaverage(cfg, erRespAll{:});
erRespTfrAll = returnStructs(data, 'er_responseLockedTFRhann');
erRespTfrAve = ft_freqgrandaverage(cfg, erRespTfrAll{:});
cfg = [];
erOnsAll = returnStructs(data, 'er_coherenceLockedAverage');
erOnsAve = ft_timelockgrandaverage(cfg, erOnsAll{:});
erOnsTfrAll = returnStructs(data, 'er_coherenceLockedTFRhann');
erOnsTfrAve = ft_freqgrandaverage(cfg, erOnsTfrAll{:});
cfg = [];
hrRespAll = returnStructs(data, 'hr_responseLockedAverage');
hrRespAve = ft_timelockgrandaverage(cfg, hrRespAll{:});
hrRespTfrAll = returnStructs(data, 'hr_responseLockedTFRhann');
hrRespTfrAve = ft_freqgrandaverage(cfg, hrRespTfrAll{:});
cfg = [];
hrOnsAll = returnStructs(data, 'hr_coherenceLockedAverage');
hrOnsAve = ft_timelockgrandaverage(cfg, hrOnsAll{:});
hrOnsTfrAll = returnStructs(data, 'hr_coherenceLockedTFRhann');
hrOnsTfrAve = ft_freqgrandaverage(cfg, hrOnsTfrAll{:});

% get an overall average
cfg = [];
rOnsAve = ft_timelockgrandaverage(cfg, erOnsAll{:}, hrOnsAll{:});
rOnsTfrAve = ft_freqgrandaverage(cfg, erOnsTfrAll{:}, hrOnsTfrAll{:});
cfg = [];
rRespAve = ft_timelockgrandaverage(cfg, erRespAll{:}, hrRespAll{:});
rRespTfrAve = ft_freqgrandaverage(cfg, erRespTfrAll{:}, hrRespTfrAll{:});

% get diffs
rOnsDiffAve = getDifference(erOnsAve,hrOnsAve);
rRespDiffAve = getDifference(ecRespAve,hcRespAve);
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
ecerRespAll = returnStructs(data, 'ecer_responseLockedAverage');
ecerRespAve = ft_timelockgrandaverage(cfg, ecerRespAll{:});
ecerRespTfrAll = returnStructs(data, 'ecer_responseLockedTFRhann');
ecerRespTfrAve = ft_freqgrandaverage(cfg, ecerRespTfrAll{:});
cfg = [];
ecerOnsAll = returnStructs(data, 'ecer_coherenceLockedAverage');
ecerOnsAve = ft_timelockgrandaverage(cfg, ecerOnsAll{:});
ecerOnsTfrAll = returnStructs(data, 'ecer_coherenceLockedTFRhann');
ecerOnsTfrAve = ft_freqgrandaverage(cfg, ecerOnsTfrAll{:});
cfg = [];
echrRespAll = returnStructs(data, 'echr_responseLockedAverage');
echrRespAve = ft_timelockgrandaverage(cfg, echrRespAll{:});
echrRespTfrAll = returnStructs(data, 'echr_responseLockedTFRhann');
echrRespTfrAve = ft_freqgrandaverage(cfg, echrRespTfrAll{:});
cfg = [];
echrOnsAll = returnStructs(data, 'echr_coherenceLockedAverage');
echrOnsAve = ft_timelockgrandaverage(cfg, echrOnsAll{:});
echrOnsTfrAll = returnStructs(data, 'echr_coherenceLockedTFRhann');
echrOnsTfrAve = ft_freqgrandaverage(cfg, echrOnsTfrAll{:});
cfg = [];
hcerRespAll = returnStructs(data, 'hcer_responseLockedAverage');
hcerRespAve = ft_timelockgrandaverage(cfg, hcerRespAll{:});
hcerRespTfrAll = returnStructs(data, 'hcer_responseLockedTFRhann');
hcerRespTfrAve = ft_freqgrandaverage(cfg, hcerRespTfrAll{:});
cfg = [];
hcerOnsAll = returnStructs(data, 'hcer_coherenceLockedAverage');
hcerOnsAve = ft_timelockgrandaverage(cfg, hcerOnsAll{:});
hcerOnsTfrAll = returnStructs(data, 'hcer_coherenceLockedTFRhann');
hcerOnsTfrAve = ft_freqgrandaverage(cfg, hcerOnsTfrAll{:});
cfg = [];
hchrRespAll = returnStructs(data, 'hchr_responseLockedAverage');
hchrRespAve = ft_timelockgrandaverage(cfg, hchrRespAll{:});
hchrRespTfrAll = returnStructs(data, 'hchr_responseLockedTFRhann');
hchrRespTfrAve = ft_freqgrandaverage(cfg, hchrRespTfrAll{:});
cfg = [];
hchrOnsAll = returnStructs(data, 'hchr_coherenceLockedAverage');
hchrOnsAve = ft_timelockgrandaverage(cfg, hchrOnsAll{:});
hchrOnsTfrAll = returnStructs(data, 'hchr_coherenceLockedTFRhann');
hchrOnsTfrAve = ft_freqgrandaverage(cfg, hchrOnsTfrAll{:});

disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get subjectwise differences %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('getting differences')

for subject = 1:numel(ecRespAll)
    fprintf('\ngetting subjectwise differences for subject %.0f of %.0f\n\n',subject,numel(ecRespAll))
    
    cOnsDiffAll{subject} = getDifference(ecOnsAll{subject}, hcOnsAll{subject});
    cRespDiffAll{subject}= getDifference(ecRespAll{subject}, hcRespAll{subject});
    rOnsDiffAll{subject} = getDifference(erOnsAll{subject}, hrOnsAll{subject});
    rRespDiffAll{subject}= getDifference(erRespAll{subject}, hrRespAll{subject});

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

% disp('getting averages')
% 
% cfg = [];
% consdiffave = ft_timelockgrandaverage(cfg, consdiffall{:});
% crespdiffave = ft_timelockgrandaverage(cfg, crespdiffall{:});
% ronsdiffave = ft_timelockgrandaverage(cfg, ronsdiffall{:});
% rrespdiffave = ft_timelockgrandaverage(cfg, rrespdiffall{:});

disp('done')

%% save those differences, if we want

disp('saving')

save(fullfile(datadir,'differences.mat'),...
    'cOnsDiffAll', 'cRespDiffAll', 'rOnsDiffAll', 'rRespDiffAll', 'cOnsDiffAve', 'cRespDiffAve', 'rOnsDiffAve', 'rRespDiffAve',...
    'cOnsTfrDiffAll', 'cRespTfrDiffAll', 'rOnsTfrDiffAll', 'rRespTfrDiffAll', 'cOnsTfrDiffAve', 'cRespTfrDiffAve', 'rOnsTfrDiffAve', 'rRespTfrDiffAve'...
    )

disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the timecourse of univariate activity %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%% coh onset in eeg
cOnsDiffAve.bfs = testDiffsAcrossTime(cOnsDiffAll,CPP);
plotTimecourse(eegLayout,[teal;coral],... % layout and colours
    {ecOnsAve,hcOnsAve},cOnsDiffAve,... % {average,data},averageDifference (with bayes factors)
    {ecOnsAll,hcOnsAll},... % {subjectwise,data}
    [-0.2 1.5],[],CPP,... % xlims,ylims,channels
    {'EEG','uV'},... % sensor type and units for ylablel
    {'easy coh';'hard coh'},'northwest','CPP (CP1 CPz CP2) in EEG: Coherence Onset',... % legend, legend location, plot title
    [erpFigDir filesep 'eeg_coh_ons_ERP.png']) % save loc
%% coh response in eeg
cRespDiffAve.bfs = testDiffsAcrossTime(cRespDiffAll,CPP);
plotTimecourse(eegLayout,[teal;coral],...
    {ecRespAve,hcRespAve},cRespDiffAve,...
    {ecRespAll,hcRespAll},...
    [],[],CPP,...
    {'EEG','uV'},...
    {'easy coh';'hard coh'},'northwest','CPP (CP1 CPz CP2) in EEG: Coherence Response',...
    [erpFigDir filesep 'eeg_coh_resp_ERP.png'])
%% rule onset in eeg
rOnsDiffAve.bfs = testDiffsAcrossTime(rOnsDiffAll,CPP);
plotTimecourse(eegLayout,[teal;coral],...
    {erOnsAve,hrOnsAve},rOnsDiffAve,...
    {erOnsAll,hrOnsAll},...
    [-0.2 1.5],[-3.5e-06 7.5e-06],CPP,...
    {'EEG','uV'},...
    {'easy cat';'hard cat'},'southeast','CPP (CP1 CPz CP2) in EEG: Categorisation Onset',...
    [erpFigDir filesep 'eeg_cat_ons_ERP.png'])
%% rule response in eeg
rRespDiffAve.bfs = testDiffsAcrossTime(rRespDiffAll,CPP);
plotTimecourse(eegLayout,[teal;coral],...
    {erRespAve,hrRespAve},rRespDiffAve,...
    {erRespAll,hrRespAll},...
    [],[-3.5e-06 7.5e-06],CPP,...
    {'EEG','uV'},...
    {'easy cat';'hard cat'},'northwest','CPP (CP1 CPz CP2) in EEG: Categorisation Response',...
    [erpFigDir filesep 'eeg_cat_resp_ERP.png'])
%% onset of all conditions in eeg
plotTimecourse(eegLayout,[lilac;teal;coral;maroon],... % layout and colours
    {ecerOnsAve, echrOnsAve, hcerOnsAve, hchrOnsAve},[],... % {average,data},averageDifference (with bayes factors)
    {ecerOnsAll,echrOnsAll,hcerOnsAll,hchrOnsAll},... % {subjectwise,data}
    [],[],CPP,... % xlims,ylims,channels
    {'EEG','uV'},... % sensor type and units for ylablel
    {'EasyCoh EasyCat';'EasyCoh HardCat';'HardCoh EasyCat';'HardCoh HardCat'},'southeast','CPP (CP1 CPz CP2) in EEG: Onset in all Conditions',... % legend, legend location, plot title
    [erpFigDir filesep 'eeg_allConds_ons_ERP.png']) % save loc
%% response of all conditions in eeg
plotTimecourse(eegLayout,[lilac;teal;coral;maroon],... % layout and colours
    {ecerRespAve, echrRespAve, hcerRespAve, hchrRespAve},[],... % {average,data},averageDifference (with bayes factors)
    {ecerRespAll,echrRespAll,hcerRespAll,hchrRespAll},... % {subjectwise,data}
    [],[],CPP,... % xlims,ylims,channels
    {'EEG','uV'},... % sensor type and units for ylablel
    {'EasyCoh EasyCat';'EasyCoh HardCat';'HardCoh EasyCat';'HardCoh HardCat'},'southeast','CPP (CP1 CPz CP2) in EEG: Response in all Conditions',... % legend, legend location, plot title
    [erpFigDir filesep 'eeg_allConds_resp_ERP.png']) % save loc

%%

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
        
minValue = @(x, pref) min(min(x.avg(startsWith(x.label, pref), :)));
maxValue = @(x, pref) max(max(x.avg(startsWith(x.label, pref), :)));

zlimEegOns = [...
    min([minValue(cOnsDiffAve,'EEG'),minValue(rOnsDiffAve,'EEG'),minValue(cOnsAve,'EEG'),minValue(rOnsAve,'EEG')]),...
    max([maxValue(cOnsDiffAve,'EEG'),maxValue(rOnsDiffAve,'EEG'),maxValue(cOnsAve,'EEG'),maxValue(rOnsAve,'EEG')])...
    ];
zlimEegResp = [...
    min([minValue(cRespDiffAve,'EEG'),minValue(rRespDiffAve,'EEG'),minValue(cRespAve,'EEG'),minValue(rRespAve,'EEG')]),...
    max([maxValue(cRespDiffAve,'EEG'),maxValue(rRespDiffAve,'EEG'),maxValue(cRespAve,'EEG'),maxValue(rRespAve,'EEG')])...
    ];
zlimMegOns = [...
    min([minValue(cOnsDiffAve,'MEG'),minValue(rOnsDiffAve,'MEG'),minValue(cOnsAve,'MEG'),minValue(rOnsAve,'MEG')]),...
    max([maxValue(cOnsDiffAve,'MEG'),maxValue(rOnsDiffAve,'MEG'),maxValue(cOnsAve,'MEG'),maxValue(rOnsAve,'MEG')])...
    ];
zlimMegResp = [...
    min([minValue(cRespDiffAve,'MEG'),minValue(rRespDiffAve,'MEG'),minValue(cRespAve,'MEG'),minValue(rRespAve,'MEG')]),...
    max([maxValue(cRespDiffAve,'MEG'),maxValue(rRespDiffAve,'MEG'),maxValue(cRespAve,'MEG'),maxValue(rRespAve,'MEG')])...
    ];


%% EEG topo of onset across coherence
doTopos({'EEG'},eegNeighbours,ecOnsAll,hcOnsAll,...
    zlimEegOns,0.05,[-0.3,1.5],[6,6],36,...
    cOnsDiffAve,eegLayout,'eeg_coh_ons_on_diff_sig_ERP_clusters','eeg_coh_ons_sig_ERP_clusters',datadir,erpFigDir)
%% EEG topo of response across coherence
doTopos({'EEG'},eegNeighbours,ecRespAll,hcRespAll,...
    zlimEegResp,0.05,[-0.6,0.2],[4,4],16,...
    cRespDiffAve,eegLayout,'eeg_coh_resp_on_diff_sig_ERP_clusters','eeg_coh_resp_sig_ERP_clusters',datadir,erpFigDir)
%% EEG topo of onset across categorisation
doTopos({'EEG'},eegNeighbours,erOnsAll,hrOnsAll,...
    zlimEegOns,0.05,[-0.3,1.5],[6,6],36,...
    rOnsDiffAve,eegLayout,'eeg_cat_ons_on_diff_sig_ERP_clusters','eeg_cat_ons_sig_ERP_clusters',datadir,erpFigDir)
%% EEG topo of response across categorisation
doTopos({'EEG'},eegNeighbours,erRespAll,hrRespAll,...
    zlimEegResp,0.05,[-0.6,0.2],[4,4],16,...
    rRespDiffAve,eegLayout,'eeg_cat_resp_on_diff_sig_ERP_clusters','eeg_cat_resp_sig_ERP_clusters',datadir,erpFigDir)
%% EEG topo of onset across coherence ON AVERAGE, NOT DIFF
doTopos({'EEG'},eegNeighbours,ecOnsAll,hcOnsAll,...
    zlimEegOns,0.05,[-0.3,1.5],[6,6],36,...
    cOnsAve,eegLayout,'eeg_coh_ons_on_ave_sig_ERP_clusters','eeg_coh_ons_sig_ERP_clusters',datadir,erpFigDir)
%% EEG topo of response across coherence ON AVERAGE, NOT DIFF
doTopos({'EEG'},eegNeighbours,ecRespAll,hcRespAll,...
    zlimEegResp,0.05,[-0.6,0.2],[4,4],16,...
    cRespAve,eegLayout,'eeg_coh_resp_on_ave_sig_ERP_clusters','eeg_resp_ons_sig_ERP_clusters',datadir,erpFigDir)
%% let's plot that frontal stuff: coherence
cOnsDiffAve.bfs = testDiffsAcrossTime(cOnsDiffAll,frontal);
plotTimecourse(eegLayout,[teal;coral],... % layout and colours
    {ecOnsAve,hcOnsAve},cOnsDiffAve,... % {average,data},averageDifference (with bayes factors)
    {ecOnsAll,hcOnsAll},... % {subjectwise,data}
    [],[],frontal,... % xlims,ylims,channels
    {'EEG','uV'},... % sensor type and units for ylablel
    {'easy coh';'hard coh'},'northwest','Frontal sensors in EEG: Coherence Onset',... % legend, legend location, plot title
    [erpFigDir filesep 'eeg_coh_ons_frontal_ERP.png']) % save loc
cRespDiffAve.bfs = testDiffsAcrossTime(cRespDiffAll,frontal);
plotTimecourse(eegLayout,[teal;coral],... % layout and colours
    {ecRespAve,hcRespAve},cRespDiffAve,... % {average,data},averageDifference (with bayes factors)
    {ecRespAll,hcRespAll},... % {subjectwise,data}
    [],[],frontal,... % xlims,ylims,channels
    {'EEG','uV'},... % sensor type and units for ylablel
    {'easy coh';'hard coh'},'northwest','Frontal sensors in EEG: Coherence Response',... % legend, legend location, plot title
    [erpFigDir filesep 'eeg_coh_resp_frontal_ERP.png']) % save loc
%% let's plot that frontal stuff: categorisation
rOnsDiffAve.bfs = testDiffsAcrossTime(rOnsDiffAll,frontal);
plotTimecourse(eegLayout,[teal;coral],... % layout and colours
    {erOnsAve,hrOnsAve},rOnsDiffAve,... % {average,data},averageDifference (with bayes factors)
    {erOnsAll,hrOnsAll},... % {subjectwise,data}
    [],[],frontal,... % xlims,ylims,channels
    {'EEG','uV'},... % sensor type and units for ylablel
    {'easy cat';'hard cat'},'northwest','Frontal sensors in EEG: Categorisation Onset',... % legend, legend location, plot title
    [erpFigDir filesep 'eeg_cat_ons_frontal_ERP.png']) % save loc
rRespDiffAve.bfs = testDiffsAcrossTime(rRespDiffAll,frontal);
plotTimecourse(eegLayout,[teal;coral],... % layout and colours
    {erRespAve,hrRespAve},rRespDiffAve,... % {average,data},averageDifference (with bayes factors)
    {erRespAll,hrRespAll},... % {subjectwise,data}
    [],[],frontal,... % xlims,ylims,channels
    {'EEG','uV'},... % sensor type and units for ylablel
    {'easy cat';'hard cat'},'northwest','Frontal sensors in EEG: Categorisation Response',... % legend, legend location, plot title
    [erpFigDir filesep 'eeg_cat_resp_frontal_ERP.png']) % save loc
%% let's plot a suss sensor: coherence
cOnsDiffAve.bfs = testDiffsAcrossTime(cOnsDiffAll,'EEG061');
plotTimecourse(eegLayout,[teal;coral],... % layout and colours
    {ecOnsAve,hcOnsAve},cOnsDiffAve,... % {average,data},averageDifference (with bayes factors)
    {ecOnsAll,hcOnsAll},... % {subjectwise,data}
    [-0.2 1.5],[],'EEG061',... % xlims,ylims,channels
    {'EEG','uV'},... % sensor type and units for ylablel
    {'easy coh';'hard coh'},'northwest','Iz in EEG: Coherence Onset',... % legend, legend location, plot title
    [erpFigDir filesep 'eeg_coh_ons_Suss.png']) % save loc

%% MEG topo of onset across coherence
doTopos({'MEG'},megNeighbours,ecOnsAll,hcOnsAll,...
    zlimMegOns,0.1,[-0.1,1.5],[4,4],16,...
    cOnsDiffAve,megLayout,'meg_coh_ons_on_diff_sig_ERP_clusters','meg_coh_ons_sig_ERP_clusters',datadir,erpFigDir)
%% MEG topo of response across coherence
doTopos({'MEG'},megNeighbours,ecRespAll,hcRespAll,...
    zlimMegResp,0.05,[-0.6,0.2],[4,4],16,...
    cRespDiffAve,megLayout,'meg_coh_resp_on_diff_sig_ERP_clusters','meg_coh_resp_sig_ERP_clusters',datadir,erpFigDir)
%% MEG topo of onset across categorisation
doTopos({'MEG'},megNeighbours,erOnsAll,hrOnsAll,...
    zlimMegOns,0.1,[-0.1,1.5],[4,4],16,...
    rOnsDiffAve,megLayout,'meg_cat_ons_on_diff_sig_ERP_clusters','meg_cat_ons_sig_ERP_clusters',datadir,erpFigDir)
%% MEG topo of response across categorisation
doTopos({'MEG'},megNeighbours,erRespAll,hrRespAll,...
    zlimMegResp,0.05,[-0.6,0.2],[4,4],16,...
    rRespDiffAve,megLayout,'meg_cat_resp_on_diff_sig_ERP_clusters','meg_cat_resp_sig_ERP_clusters',datadir,erpFigDir)


%% subfunctions

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
    % if multiple channels, get the average, else get the channel
    if size(dataStruct{i}.avg(find(ismember(dataStruct{i}.label,channels)),:),1) > 1
        diffs(:,i) = mean(dataStruct{i}.avg(find(ismember(dataStruct{i}.label,channels)),:));
    else
        diffs(:,i) = dataStruct{i}.avg(find(ismember(dataStruct{i}.label,channels)),:);
    end
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




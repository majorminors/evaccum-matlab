function b2_aggregate_data(thisSubject,datadir,toolsdir,overwrite)
% here we will:
% 1. produce a file for each subject with the timelocked averages we care
% about for analysis
% 2. produce a file for each subject with the TFRs

if ~exist('overwrite','var'); overwrite = 0; end

% set up paths
addpath /hpc-software/matlab/cbu/
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/');
ftDir = fullfile(toolsdir,'..','..','..','Toolboxes','fieldtrip'); % I pass this to my saving function too
addpath(ftDir); ft_defaults;
addpath(genpath(fullfile(toolsdir,'lib')))
rootDir = fullfile(datadir,thisSubject.id);
preProcDir = fullfile(rootDir,'Preprocess');
addpath(rootDir);

% file definition
% sprintf to enter missing data ('%s')
inputFilePattern = fullfile(preProcDir,'run*_if_transmid.fif');
erpOutputFilename = fullfile(preProcDir,'timelocked_averages.mat');
tfrOutputFilename = fullfile(preProcDir,'tfr_hanning.mat');
trlfilePattern = fullfile(preProcDir, '%s_trl.mat'); % this contains our trial related info for epoching - should exist (generated with megTriggers script)
behavfile = fullfile(datadir,'behavioural',[thisSubject.id '_MEGRTs.mat']); % grab our behavioural data + output from a3_megtriggers


%% print subject for logging

disp('>>> this is subject:')
disp(thisSubject.id);

if exist(erpOutputFilename,'file') && ~overwrite
    warning('>>> save file exists---skipping')
    return
elseif exist(erpOutputFilename,'file') && overwrite
    warning('>>> file exists and overwrite: ill delete is so it wont get loaded into your next analysis')
    system(['rm -f ' erpOutputFilename])
end

theseFiles = dir(inputFilePattern);
% skip this participant if we don't have any files to work with
if isempty(theseFiles)
    warning('>>> no files match this inputFilePattern---skipping this participant');
    return
else
    fprintf('>>> found %.0f files\n',numel(theseFiles));
end

% loop through those files
for fileNum = 1:numel(theseFiles)
    % let's make the current data file easy to reference
    thisFile = [theseFiles(fileNum).folder filesep theseFiles(fileNum).name];
    % let's find what the run label is (e.g. 'run1' etc) (looks for
    % 'run' followed by any amount of numerical digits, and returns
    % {'run[numbers]'}, so we pop it out of the cell array. This would
    % break if you had more than one 'run[numbers]' in your filename
    % because it would return more than one item in the cell array
    runLabel = regexp(thisFile,'run\d+','match'); runLabel = runLabel{:};
    % and get the run number in a similar way, since I later found out I need this too
    thisRun = regexp(thisFile,'run(\d*)_','tokens'); thisRun = thisRun{1}{1};
    % now we can find the trial information for this run
    thisTrlFile = sprintf(trlfilePattern,runLabel);
    
    fprintf('>>> loading data %.0f of %.0f\n',fileNum,numel(theseFiles))
    % ok, now let's load in the data file
    cfg = [];
    cfg.continuous = 'yes'; % this data is not epoched, it's continuous
    if endsWith(thisFile,'.fif')
        cfg.dataset = thisFile;
        rawData = convertFif(cfg,ftDir); % see function: fieldtrip clears the mne toolbox after first use? so let's force this in a distinct environment, in case this is something they did on purpose
    elseif endsWith(thisFile,'.mat')
        rawData = ft_preprocessing(cfg,loadFtData(thisFile)); % load it
    end
    
    % add subject number to the thing
    rawData.subj = thisSubject.num;
    % see: https://www.fieldtriptoolbox.org/development/datastructure/
    % but unfortunately this does not carry through to future stuff!
    
    % let's do any additional filtering
    cfg = [];
    cfg.continuous = 'yes'; % this data is not epoched, it's continuous
%     cfg.hpfilter        = 'yes';
%     cfg.hpfreq          = 0.5;
%     cfg.lpfilter        = 'yes';
%     cfg.lpfreq          = 200;
    rawData = ft_preprocessing(cfg,rawData);
    
    % let's make a layout using the data
    disp('>>> prepping layouts')
    % make an eeg layout
    cfg = [];
    cfg.elec = rawData.elec;
    cfg.output = fullfile(preProcDir,'eeg_layout.lay');
    eegLayout = ft_prepare_layout(cfg);
    % and an meg layout
    cfg = [];
    cfg.grad = rawData.grad;
    cfg.output = fullfile(preProcDir,'meg_layout.lay');
    megLayout = ft_prepare_layout(cfg);
    % combine them
    cfg = [];
    combinedLayout = ft_appendlayout(cfg, eegLayout, megLayout);
    % won't bother saving this one

    %% now onset of coherent dots
    disp('>>> compiling dataset based on onset')
    cfg = [];
    cfg.trlfile = thisTrlFile;
    cfg.behavfile = behavfile;
    cfg.pre = -500;
    cfg.run = thisRun;
    cfg.trl = trlFromFile(cfg);
    cohOnsetsData{fileNum} = ft_redefinetrial(cfg,rawData);    
    % baseline correction
    cfg = [];
    cfg.demean          = 'yes';
    cfg.baselinewindow  = [-0.2 0]; % assumes we're demeaning coherence-locked stimuli
    cohOnsetsData{fileNum} = ft_preprocessing(cfg,cohOnsetsData{fileNum});
    % reject artefacts
    cohOnsetsData{fileNum} = rejectArtefacts(cohOnsetsData{fileNum}, combinedLayout,ftDir);
    
    %% get response locked data
    disp('>>> compiling dataset based on responses')
    cfg = [];
    cfg.trlfile = thisTrlFile;
    cfg.behavfile = behavfile;
    cfg.run = thisRun;
    cfg.lockTo = 'response';
    cfg.pre = -600;
    cfg.post = 200;
    cfg.trl = trlFromFile(cfg);
    % if we want to look at the raw data
%     respLockedData{fileNum} = ft_redefinetrial(cfg,rawData);
    % if we want to look at the demeaned data
    respLockedData{fileNum} = ft_redefinetrial(cfg,cohOnsetsData{fileNum});
    % reject artefacts
    respLockedData{fileNum} = rejectArtefacts(respLockedData{fileNum}, combinedLayout,ftDir);

    
    
end

clear rawData

%% append these now

disp('>>> appending datasets')

% append them all
cfg = [];
cfg.keepsampleinfo  = 'yes';
coherenceOnsetData = ft_appenddata(cfg,cohOnsetsData{:}); clear cohOnsetsData
responseLockedData = ft_appenddata(cfg,respLockedData{:}); clear respLockedData

% % we can compile the overall averages if we want like so:
% disp('>>> compiling overall averages')
% 
% cfg = [];
% cfg.trials = find(coherenceOnsetData.trialinfo(:,2) == 1); % grab just the good trials
% coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetData);
% 
% cfg = [];
% cfg.trials = find(responseLockedData.trialinfo(:,2) == 1); % grab just the good trials
% responseLockedAverage = ft_timelockanalysis(cfg, responseLockedData);


%% compile averages broken into our conditions - first response-locked

disp('>>> compiling condition-wise averages')

cfg = [];
cfg.trials = find(responseLockedData.trialinfo(:,2) == 1 & responseLockedData.trialinfo(:,1) == 1); % good trials, 1st condition
ecer_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedData);
ecer_responseLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials), responseLockedData);

% cfg = [];
% % cfg.baseline     = [-0.5 -0.4];
% cfg.baseline     = [-0.2 0];
% cfg.baselinetype = 'absolute';
% % cfg.zlim         = [-2.5e-24 2.5e-24];
% cfg.showlabels   = 'yes';
% cfg.layout       = megLayout;
% cfg.maskstyle = 'saturation';
% % cfg.xlim=[-0.5 1.5];
% % cfg.channel = getMegLabels('parietal');
% ft_multiplotTFR(cfg, ecer_responseLockedTFRhann);
% ft_singleplotTFR(cfg, ecer_coherenceLockedTFRhann);
% ft_topoplotTFR(cfg, ecer_responseLockedTFRhann);
% % ecer_coherenceLockedTFRhann
% % ecer_responseLockedTFRhann

cfg = [];
cfg.trials = find(responseLockedData.trialinfo(:,2) == 1 & responseLockedData.trialinfo(:,1) == 2); % good trials, 2nd condition
echr_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedData);
echr_responseLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials), responseLockedData);

cfg = [];
cfg.trials = find(responseLockedData.trialinfo(:,2) == 1 & responseLockedData.trialinfo(:,1) == 3); % good trials, 3rd condition
hcer_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedData);
hcer_responseLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials), responseLockedData);

cfg = [];
cfg.trials = find(responseLockedData.trialinfo(:,2) == 1 & responseLockedData.trialinfo(:,1) == 4); % good trials, 4th condition
hchr_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedData);
hchr_responseLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials), responseLockedData);

%% and averaged across manipulations - response-locked

cfg = [];
cfg.trials = find(responseLockedData.trialinfo(:,2) == 1 & (responseLockedData.trialinfo(:,1) == 1 | responseLockedData.trialinfo(:,1) == 2)); % good trials, easy coherenc3
ec_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedData);
ec_responseLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials), responseLockedData);

cfg = [];
cfg.trials = find(responseLockedData.trialinfo(:,2) == 1 & (responseLockedData.trialinfo(:,1) == 3 | responseLockedData.trialinfo(:,1) == 4)); % good trials, hard coherence
hc_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedData);
hc_responseLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials), responseLockedData);

cfg = [];
cfg.trials = find(responseLockedData.trialinfo(:,2) == 1 & (responseLockedData.trialinfo(:,1) == 1 | responseLockedData.trialinfo(:,1) == 3)); % good trials, easy rule
er_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedData);
er_responseLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials), responseLockedData);

cfg = [];
cfg.trials = find(responseLockedData.trialinfo(:,2) == 1 & (responseLockedData.trialinfo(:,1) == 2 | responseLockedData.trialinfo(:,1) == 4)); % good trials, hard rule
hr_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedData);
hr_responseLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials), responseLockedData);

%% same again for conditions but onset-locked

cfg = [];
cfg.trials = find(coherenceOnsetData.trialinfo(:,2) == 1 & coherenceOnsetData.trialinfo(:,1) == 1); % good trials, 1st condition
ecer_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetData);
ecer_coherenceLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials), coherenceOnsetData);

cfg = [];
cfg.trials = find(coherenceOnsetData.trialinfo(:,2) == 1 & coherenceOnsetData.trialinfo(:,1) == 2); % good trials, 2nd condition
echr_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetData);
echr_coherenceLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials), coherenceOnsetData);

cfg = [];
cfg.trials = find(coherenceOnsetData.trialinfo(:,2) == 1 & coherenceOnsetData.trialinfo(:,1) == 3); % good trials, 3rd condition
hcer_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetData);
hcer_coherenceLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials), coherenceOnsetData);

cfg = [];
cfg.trials = find(coherenceOnsetData.trialinfo(:,2) == 1 & coherenceOnsetData.trialinfo(:,1) == 4); % good trials, 4th condition
hchr_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetData);
hchr_coherenceLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials), coherenceOnsetData);

%% and again averaged across manipulations - onset-locked

cfg = [];
cfg.trials = find(coherenceOnsetData.trialinfo(:,2) == 1 & (coherenceOnsetData.trialinfo(:,1) == 1 | coherenceOnsetData.trialinfo(:,1) == 2)); % good trials, easy coherence
ec_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetData);
ec_coherenceLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials), coherenceOnsetData);

cfg = [];
cfg.trials = find(coherenceOnsetData.trialinfo(:,2) == 1 & (coherenceOnsetData.trialinfo(:,1) == 3 | coherenceOnsetData.trialinfo(:,1) == 4)); % good trials, hard coherence
hc_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetData);
hc_coherenceLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials), coherenceOnsetData);

cfg.trials = find(coherenceOnsetData.trialinfo(:,2) == 1 & (coherenceOnsetData.trialinfo(:,1) == 1 | coherenceOnsetData.trialinfo(:,1) == 3)); % good trials, easy rule
er_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetData);
er_coherenceLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials), coherenceOnsetData);

cfg = [];
cfg.trials = find(coherenceOnsetData.trialinfo(:,2) == 1 & (coherenceOnsetData.trialinfo(:,1) == 2 | coherenceOnsetData.trialinfo(:,1) == 4)); % good trials, hard rule
hr_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetData);
hr_coherenceLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials), coherenceOnsetData);

%% save!

disp('>>> saving')
% save(erpOutputFilename,'coherenceLockedAverage','responseLockedAverage',...
save(erpOutputFilename,...
    'ecer_responseLockedAverage','echr_responseLockedAverage','hcer_responseLockedAverage','hchr_responseLockedAverage',...
    'ec_responseLockedAverage','hc_responseLockedAverage','er_responseLockedAverage','hr_responseLockedAverage',...
    'ecer_coherenceLockedAverage','echr_coherenceLockedAverage','hcer_coherenceLockedAverage','hchr_coherenceLockedAverage',...
    'ec_coherenceLockedAverage','hc_coherenceLockedAverage','er_coherenceLockedAverage','hr_coherenceLockedAverage',...
    '-v7.3');
save(tfrOutputFilename,...
    'ecer_responseLockedTFRhann','echr_responseLockedTFRhann','hcer_responseLockedTFRhann','hchr_responseLockedTFRhann',...
    'ec_responseLockedTFRhann','hc_responseLockedTFRhann','er_responseLockedTFRhann','hr_responseLockedTFRhann',...
    'ecer_coherenceLockedTFRhann','echr_coherenceLockedTFRhann','hcer_coherenceLockedTFRhann','hchr_coherenceLockedTFRhann',...
    'ec_coherenceLockedTFRhann','hc_coherenceLockedTFRhann','er_coherenceLockedTFRhann','hr_coherenceLockedTFRhann',...
    '-v7.3');
disp('>>> saved')

disp('>>> done')
    


return
end

function cleanedData = rejectArtefacts(inputData, layout,ftDir)
addpath(ftDir); ft_defaults;

disp('>>> doing auto-artefact rejection')

cfg = [];
% here's an example of badchannels (NO CHANNELS SPECIFIED!)
% but we won't worry about that for now
%         cfg.layout = layout;
%         cfg.method = 'distance';
%         cfg.neighbours = ft_prepare_neighbours(cfg, inputData);
%         cfg = ft_badchannel(cfg, inputData);
%         badMegChans = cfg.badchannel;
cfg.artfctdef.reject = 'complete'; % remove 'complete' trials
cfg.metric = 'zvalue';
cfg.threshold = 4;
cfg.channel = 'MEG';
cfg = ft_badsegment(cfg, inputData);
cleanedData = ft_rejectartifact(cfg, inputData);
cfg.channel = 'EEG';
cfg = ft_badsegment(cfg, cleanedData);
disp('dont worry about warnings about nans---this wont spread, its contained to the trial youre naning')
cleanedData = ft_rejectartifact(cfg, cleanedData);

return
end

function convertedData = convertFif(cfg, ftDir)
% fieldtrip clears the mne toolbox after first use? so let's force this in
% a distinct environment in case fieldtrip does this for a good reason

addpath(ftDir); ft_defaults; addpath(fullfile(ftDir, 'external','mne'));
convertedData = ft_preprocessing(cfg); % load it
            
return
end



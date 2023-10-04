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
epochedOutputFilename = fullfile(preProcDir,'epoched_data.mat');
erpOutputFilename = fullfile(preProcDir,'timelocked_averages.mat');
tfrLowOutputFilename = fullfile(preProcDir,'tfr_hanning.mat');
tfrHighOutputFilename = fullfile(preProcDir,'tfr_multi.mat');
trlfilePattern = fullfile(preProcDir, '%s_trl.mat'); % this contains our trial related info for epoching - should exist (generated with megTriggers script)
behavfile = fullfile(datadir,'behavioural',[thisSubject.id '_MEGRTs.mat']); % grab our behavioural data + output from a3_megtriggers


%% print subject for logging

disp('>>> this is subject:')
disp(thisSubject.id);

if exist(epochedOutputFilename,'file') && ~overwrite
    warning('>>> epoched data save file exists and overwrite off---skipping erp')
    doEpoch = 0;
elseif ~exist(epochedOutputFilename,'file') || overwrite
    if exist(epochedOutputFilename,'file')
        warning('>>> epoched data exists and overwrite is on: ill delete is so it wont get loaded into your next analysis')
        system(['rm -f ' epochedOutputFilename])
    end
    doEpoch = 1;
end
if exist(erpOutputFilename,'file') && ~overwrite
    warning('>>> erp save file exists and overwrite off---skipping erp')
    doErp = 0;
elseif ~exist(erpOutputFilename,'file') || overwrite
    if exist(erpOutputFilename,'file')
        warning('>>> erp file exists and overwrite is on: ill delete is so it wont get loaded into your next analysis')
        system(['rm -f ' erpOutputFilename])
    end
    doErp = 1;
end
if exist(tfrLowOutputFilename,'file') && ~overwrite
    warning('>>> erp save file exists and overwrite off---skipping trf low')
    doTfrLow = 0;
elseif ~exist(tfrLowOutputFilename,'file') || overwrite
    if exist(tfrLowOutputFilename,'file')
        warning('>>> tfr low file exists and overwrite is on: ill delete is so it wont get loaded into your next analysis')
        system(['rm -f ' tfrLowOutputFilename])
    end
    doTfrLow = 1;
end
if exist(tfrHighOutputFilename,'file') && ~overwrite
    warning('>>> erp save file exists and overwrite off---skipping trf high')
    doTfrHigh = 0;
elseif ~exist(tfrHighOutputFilename,'file') || overwrite
    if exist(tfrHighOutputFilename,'file')
        warning('>>> tfr high file exists and overwrite is on: ill delete is so it wont get loaded into your next analysis')
        system(['rm -f ' tfrHighOutputFilename])
    end
    doTfrHigh = 0;
end

if ~doEpoch && ~doErp && ~doTfrLow && ~doTfrHigh
    disp('nothing to do for this participant: skipping')
    return
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
    cfg.lpfilter        = 'yes';
    cfg.lpfreq          = 200; % cut off everything above high gamma
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

    %% now let's epoch our data around the onset and the response for our erp analysis
    
    % first to onset
    disp('>>> compiling erp dataset locked to onset')
    cfg = [];
    cfg.trlfile = thisTrlFile;
    cfg.behavfile = behavfile;
    cfg.pre = -500; % how much pre onset
    cfg.run = thisRun;
    cfg.trl = trlFromFile(cfg); % my trl is by default locked to the trial onset to end, so we're just looking at the trial plus 500ms prior
    cohOnsetErpData{fileNum} = ft_redefinetrial(cfg,rawData);    
    % baseline correction
    cfg = [];
    cfg.demean          = 'yes';
    cfg.baselinewindow  = [-0.2 0]; % assumes we're demeaning coherence-locked stimuli
    cohOnsetErpData{fileNum} = ft_preprocessing(cfg,cohOnsetErpData{fileNum});
    % reject artefacts
    cohOnsetErpData{fileNum} = rejectArtefacts(cohOnsetErpData{fileNum}, combinedLayout,ftDir);
    
    % now to response
    disp('>>> compiling erp dataset locked to response')
    cfg = [];
    cfg.trlfile = thisTrlFile;
    cfg.behavfile = behavfile;
    cfg.run = thisRun;
    cfg.lockTo = 'response'; % my trl function will use this to adjust the trl to the response, and you need to specify how much time pre and post response you want
    cfg.pre = -600; % pre response
    cfg.post = 200; % post response (when locked to response)
    cfg.trl = trlFromFile(cfg);
    % if we want to look at the raw data
%     respLockedErpData{fileNum} = ft_redefinetrial(cfg,rawData);
    % if we want to look at the demeaned data
    respLockedErpData{fileNum} = ft_redefinetrial(cfg,cohOnsetErpData{fileNum});
    % reject artefacts
    respLockedErpData{fileNum} = rejectArtefacts(respLockedErpData{fileNum}, combinedLayout,ftDir);

    %% now let's do the same thing for tfr analysis: onset and response locked epochs
    % we want more padding for this
    
    % first to onset
    disp('>>> compiling tfr dataset locked to onset')
    cfg = [];
    cfg.trlfile = thisTrlFile;
    cfg.behavfile = behavfile;
    % so now we want a little more of the trial, for purposes of padding
    cfg.pre = -1000; % so 1 sec prior to onset
    cfg.post = 500; % and 500 ms post trial end (leaving us with 3 sec epochs)
    cfg.run = thisRun;
    cfg.trl = trlFromFile(cfg);
    cohOnsetTfrData{fileNum} = ft_redefinetrial(cfg,rawData);    
    % baseline correction
    cfg = [];
    cfg.demean          = 'yes'; % now, as in the ft tutorial, we demean based on the whole trial
    cohOnsetTfrData{fileNum} = ft_preprocessing(cfg,cohOnsetTfrData{fileNum});
    % reject artefacts
    cohOnsetTfrData{fileNum} = rejectArtefacts(cohOnsetTfrData{fileNum}, combinedLayout,ftDir);
    
    % then response
    disp('>>> compiling tfr dataset locked to response')
    cfg = [];
    cfg.trlfile = thisTrlFile;
    cfg.behavfile = behavfile;
    cfg.run = thisRun;
    cfg.lockTo = 'response'; % my trl function will use this to adjust the trl to the response, and you need to specify how much time pre and post response you want
    cfg.pre = -1500; % 1500 ms pre response
    cfg.post = 1500; % 1500 ms post response, so now we have a 3 second epoch
    cfg.trl = trlFromFile(cfg);
    % if we want to look at the raw data
    respLockedTfrData{fileNum} = ft_redefinetrial(cfg,rawData);
    % again as in the tutorial, we'll demean based on the whole trial
    cfg = [];
    cfg.demean          = 'yes'; 
    respLockedTfrData{fileNum} = ft_preprocessing(cfg,respLockedTfrData{fileNum});
    % reject artefacts
    respLockedTfrData{fileNum} = rejectArtefacts(respLockedTfrData{fileNum}, combinedLayout,ftDir);
    
    
    
end

clear rawData

%% append these now

disp('>>> appending datasets')

% append them all
cfg = [];
cfg.keepsampleinfo  = 'yes';
coherenceOnsetErpData = ft_appenddata(cfg,cohOnsetErpData{:}); clear cohOnsetErpData
responseLockedErpData = ft_appenddata(cfg,respLockedErpData{:}); clear respLockedErpData
coherenceOnsetTfrData = ft_appenddata(cfg,cohOnsetTfrData{:}); clear cohOnsetTfrData
responseLockedTfrData = ft_appenddata(cfg,respLockedTfrData{:}); clear respLockedTfrData

% saving this produces huge save files
% if doEpoch
%     disp('>>> saving epoched data')
%     save(epochedOutputFilename,...
%         'coherenceOnsetErpData','responseLockedErpData','coherenceOnsetTfrData','responseLockedTfrData',...
%         '-v7.3');
% end

% % we can compile the overall averages if we want like so:
% disp('>>> compiling overall averages')
% 
% cfg = [];
% cfg.trials = find(coherenceOnsetErpData.trialinfo(:,2) == 1); % grab just the good trials
% coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetErpData);
% 
% cfg = [];
% cfg.trials = find(responseLockedErpData.trialinfo(:,2) == 1); % grab just the good trials
% responseLockedAverage = ft_timelockanalysis(cfg, responseLockedErpData);


%% compile averages broken into our conditions

% a quick indexing function
getTrialIndicesConditions = @(x,condCode) find(x.trialinfo(:,2) == 1 & x.trialinfo(:,1) == condCode);
% trialinfo(:,2), from trlFromFile() is a code for good vs bad trials based on the behavioural data coding
% trialinfo(:,1), is the condition codes 1:4---again from trlFromFile

% first response-locked

disp('>>> compiling condition-wise averages')

if doErp
    cfg = [];
    cfg.trials = getTrialIndicesConditions(responseLockedErpData,1); % good trials, 1st condition
    ecer_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedErpData);
    cfg.trials = getTrialIndicesConditions(responseLockedErpData,2); % good trials, 2nd condition
    echr_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedErpData);
    cfg.trials = getTrialIndicesConditions(responseLockedErpData,3); % good trials, 3rd condition
    hcer_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedErpData);
    cfg.trials = getTrialIndicesConditions(responseLockedErpData,4); % good trials, 4th condition
    hchr_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedErpData);
end

if doTfrLow
    cfg = [];
    cfg.trials = getTrialIndicesConditions(responseLockedTfrData,1); % good trials, 1st condition
    ecer_responseLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials, 'low'), responseLockedTfrData);
    cfg.trials = getTrialIndicesConditions(responseLockedTfrData,2); % good trials, 2nd condition
    echr_responseLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials, 'low'), responseLockedTfrData);
    cfg.trials = getTrialIndicesConditions(responseLockedTfrData,3); % good trials, 3rd condition
    hcer_responseLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials, 'low'), responseLockedTfrData);
    cfg.trials = getTrialIndicesConditions(responseLockedTfrData,4); % good trials, 4th condition
    hchr_responseLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials, 'low'), responseLockedTfrData);
end

if doTfrHigh
    cfg = [];
    cfg.trials = getTrialIndicesConditions(responseLockedTfrData,1); % good trials, 1st condition
    ecer_responseLockedTFRmulti = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials, 'high'), responseLockedTfrData);
    cfg.trials = getTrialIndicesConditions(responseLockedTfrData,2); % good trials, 2nd condition
    echr_responseLockedTFRmulti = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials, 'high'), responseLockedTfrData);
    cfg.trials = getTrialIndicesConditions(responseLockedTfrData,3); % good trials, 3rd condition
    hcer_responseLockedTFRmulti = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials, 'high'), responseLockedTfrData);
    cfg.trials = getTrialIndicesConditions(responseLockedTfrData,4); % good trials, 4th condition
    hchr_responseLockedTFRmulti = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials, 'high'), responseLockedTfrData);
end

% then onset locked

if doErp
    cfg = [];
    cfg.trials = getTrialIndicesConditions(coherenceOnsetErpData,1); % good trials, 1st condition
    ecer_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetErpData);
    cfg.trials = getTrialIndicesConditions(coherenceOnsetErpData,2); % good trials, 2nd condition
    echr_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetErpData);
    cfg.trials = getTrialIndicesConditions(coherenceOnsetErpData,3); % good trials, 3rd condition
    hcer_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetErpData);
    cfg.trials = getTrialIndicesConditions(coherenceOnsetErpData,4); % good trials, 4th condition
    hchr_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetErpData);
end

if doTfrLow
    cfg = [];
    cfg.trials = getTrialIndicesConditions(coherenceOnsetTfrData, 1); % good trials, 1st condition
    ecer_coherenceLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials, 'low'), coherenceOnsetTfrData);
    cfg.trials = getTrialIndicesConditions(coherenceOnsetTfrData, 2); % good trials, 2nd condition
    echr_coherenceLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials, 'low'), coherenceOnsetTfrData);
    cfg.trials = getTrialIndicesConditions(coherenceOnsetTfrData, 3); % good trials, 3rd condition
    hcer_coherenceLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials, 'low'), coherenceOnsetTfrData);
    cfg.trials = getTrialIndicesConditions(coherenceOnsetTfrData, 4); % good trials, 4th condition
    hchr_coherenceLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials, 'low'), coherenceOnsetTfrData);
end

if doTfrHigh
    cfg = [];
    cfg.trials = getTrialIndicesConditions(coherenceOnsetTfrData, 1); % good trials, 1st condition
    ecer_coherenceLockedTFRmulti = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials, 'high'), coherenceOnsetTfrData);
    cfg.trials = getTrialIndicesConditions(coherenceOnsetTfrData, 2); % good trials, 2nd condition
    echr_coherenceLockedTFRmulti = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials, 'high'), coherenceOnsetTfrData);
    cfg.trials = getTrialIndicesConditions(coherenceOnsetTfrData, 3); % good trials, 3rd condition
    hcer_coherenceLockedTFRmulti = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials, 'high'), coherenceOnsetTfrData);
    cfg.trials = getTrialIndicesConditions(coherenceOnsetTfrData, 4); % good trials, 4th condition
    hchr_coherenceLockedTFRmulti = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials, 'high'), coherenceOnsetTfrData);
end

%% now averaged across manipulations

% first response-locked

% another quick indexing function
getTrialIndicesManipulations = @(x,condCode) find(x.trialinfo(:,2) == 1 & (x.trialinfo(:,1) == condCode(1) | x.trialinfo(:,1) == condCode(2)));
% trialinfo(:,2), from trlFromFile() is a code for good vs bad trials based on the behavioural data coding
% trialinfo(:,1), is the condition codes 1:4---again from trlFromFile, so we'll pass a vector with two numbers to filter

if doErp
    cfg = [];
    cfg.trials = getTrialIndicesManipulations(responseLockedErpData,[1,2]); % good trials, easy coherence
    ec_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedErpData);
    cfg.trials = getTrialIndicesManipulations(responseLockedErpData,[3,4]); % good trials, hard coherence
    hc_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedErpData);
    cfg.trials = getTrialIndicesManipulations(responseLockedErpData,[1,3]); % good trials, easy rule
    er_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedErpData);
    cfg.trials = getTrialIndicesManipulations(responseLockedErpData,[2,4]); % good trials, hard rule
    hr_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedErpData);
end

if doTfrLow
    cfg = [];
    cfg.trials = getTrialIndicesManipulations(responseLockedTfrData,[1,2]); % good trials, easy coherence
    ec_responseLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials, 'low'), responseLockedTfrData);
    cfg.trials = getTrialIndicesManipulations(responseLockedTfrData,[3,4]); % good trials, hard coherence
    hc_responseLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials, 'low'), responseLockedTfrData);
    cfg.trials = getTrialIndicesManipulations(responseLockedTfrData,[1,3]); % good trials, easy rule
    er_responseLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials, 'low'), responseLockedTfrData);
    cfg.trials = getTrialIndicesManipulations(responseLockedTfrData,[2,4]); % good trials, hard rule
    hr_responseLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials, 'low'), responseLockedTfrData);
end

if doTfrHigh
    cfg = [];
    cfg.trials = getTrialIndicesManipulations(responseLockedTfrData,[1,2]); % good trials, easy coherence
    ec_responseLockedTFRmulti = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials, 'high'), responseLockedTfrData);
    cfg.trials = getTrialIndicesManipulations(responseLockedTfrData,[3,4]); % good trials, hard coherence
    hc_responseLockedTFRmulti = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials, 'high'), responseLockedTfrData);
    cfg.trials = getTrialIndicesManipulations(responseLockedTfrData,[1,3]); % good trials, easy rule
    er_responseLockedTFRmulti = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials, 'high'), responseLockedTfrData);
    cfg.trials = getTrialIndicesManipulations(responseLockedTfrData,[2,4]); % good trials, hard rule
    hr_responseLockedTFRmulti = ft_freqanalysis(ftTfrAnalysisConfig('response',cfg.trials, 'high'), responseLockedTfrData);
end


% and onset-locked

if doErp
    cfg = [];
    cfg.trials = getTrialIndicesManipulations(coherenceOnsetErpData,[1,2]); % good trials, easy coherence
    ec_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetErpData);
    cfg.trials = getTrialIndicesManipulations(coherenceOnsetErpData,[3,4]); % good trials, hard coherence
    hc_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetErpData);
    cfg.trials = getTrialIndicesManipulations(coherenceOnsetErpData,[1,3]); % good trials, easy rule
    er_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetErpData);
    cfg.trials = getTrialIndicesManipulations(coherenceOnsetErpData,[2,4]); % good trials, hard rule
    hr_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetErpData);
end

if doTfrLow
    cfg = [];
    cfg.trials = getTrialIndicesManipulations(coherenceOnsetTfrData,[1,2]); % good trials, easy coherence
    ec_coherenceLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials, 'low'), coherenceOnsetTfrData);
    cfg.trials = getTrialIndicesManipulations(coherenceOnsetTfrData,[3,4]); % good trials, hard coherence
    hc_coherenceLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials, 'low'), coherenceOnsetTfrData);
    cfg.trials = getTrialIndicesManipulations(coherenceOnsetTfrData,[1,3]); % good trials, easy rule
    er_coherenceLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials, 'low'), coherenceOnsetTfrData);
    cfg.trials = getTrialIndicesManipulations(coherenceOnsetTfrData,[2,4]); % good trials, hard rule
    hr_coherenceLockedTFRhann = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials, 'low'), coherenceOnsetTfrData);
end

if doTfrHigh
    cfg = [];
    cfg.trials = getTrialIndicesManipulations(coherenceOnsetTfrData,[1,2]); % good trials, easy coherence
    ec_coherenceLockedTFRmulti = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials, 'high'), coherenceOnsetTfrData);
    cfg.trials = getTrialIndicesManipulations(coherenceOnsetTfrData,[3,4]); % good trials, hard coherence
    hc_coherenceLockedTFRmulti = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials, 'high'), coherenceOnsetTfrData);
    cfg.trials = getTrialIndicesManipulations(coherenceOnsetTfrData,[1,3]); % good trials, easy rule
    er_coherenceLockedTFRmulti = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials, 'high'), coherenceOnsetTfrData);
    cfg.trials = getTrialIndicesManipulations(coherenceOnsetTfrData,[2,4]); % good trials, hard rule
    hr_coherenceLockedTFRmulti = ft_freqanalysis(ftTfrAnalysisConfig('onset',cfg.trials, 'high'), coherenceOnsetTfrData);
end

%% save!

disp('>>> saving')
% save(erpOutputFilename,'coherenceLockedAverage','responseLockedAverage',...
if doErp
    disp('>>> saving erp data')
    save(erpOutputFilename,...
        'ecer_responseLockedAverage','echr_responseLockedAverage','hcer_responseLockedAverage','hchr_responseLockedAverage',...
        'ec_responseLockedAverage','hc_responseLockedAverage','er_responseLockedAverage','hr_responseLockedAverage',...
        'ecer_coherenceLockedAverage','echr_coherenceLockedAverage','hcer_coherenceLockedAverage','hchr_coherenceLockedAverage',...
        'ec_coherenceLockedAverage','hc_coherenceLockedAverage','er_coherenceLockedAverage','hr_coherenceLockedAverage',...
        '-v7.3');
    disp('>>> done saving erp data')
end
if doTfrLow
    disp('>>> saving tfr low data')
    save(tfrLowOutputFilename,...
        'ecer_responseLockedTFRhann','echr_responseLockedTFRhann','hcer_responseLockedTFRhann','hchr_responseLockedTFRhann',...
        'ec_responseLockedTFRhann','hc_responseLockedTFRhann','er_responseLockedTFRhann','hr_responseLockedTFRhann',...
        'ecer_coherenceLockedTFRhann','echr_coherenceLockedTFRhann','hcer_coherenceLockedTFRhann','hchr_coherenceLockedTFRhann',...
        'ec_coherenceLockedTFRhann','hc_coherenceLockedTFRhann','er_coherenceLockedTFRhann','hr_coherenceLockedTFRhann',...
        '-v7.3');
    disp('>>> done saving tfr low data')
end
if doTfrHigh
    disp('>>> saving tfr high data')
    save(tfrHighOutputFilename,...
        'ecer_responseLockedTFRmulti','echr_responseLockedTFRmulti','hcer_responseLockedTFRmulti','hchr_responseLockedTFRmulti',...
        'ec_responseLockedTFRmulti','hc_responseLockedTFRmulti','er_responseLockedTFRmulti','hr_responseLockedTFRmulti',...
        'ecer_coherenceLockedTFRmulti','echr_coherenceLockedTFRmulti','hcer_coherenceLockedTFRmulti','hchr_coherenceLockedTFRmulti',...
        'ec_coherenceLockedTFRmulti','hc_coherenceLockedTFRmulti','er_coherenceLockedTFRmulti','hr_coherenceLockedTFRmulti',...
        '-v7.3');
    disp('>>> done saving tfr high data')
end
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



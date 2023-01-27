function a4_preProc01_atypical_artefacts(thisSubject,datadir,toolsdir,overwrite)
% here we will:
% 1. view the data if we want (just uncomment and run locally)
% 2. do some basic artefact detection prior to ica


if ~exist('overwrite','var'); overwrite = 0; end

% set up paths
addpath /hpc-software/matlab/cbu/
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/');
addpath(fullfile(toolsdir,'..','..','..','Toolboxes','fieldtrip')); ft_defaults;
addpath(genpath(fullfile(toolsdir,'lib')))
rootDir = fullfile(datadir,thisSubject.id);
addpath(rootDir);

% file definition
% here I'm using the transmid output of my maxfiltering---this is the
% run transformed to the coordinates of the middle run, so all my runs
% are comparable
inputFileName = fullfile(rootDir,'MaxfilterOutput','%s_transmid.fif'); % I'll use thisSubject.meg_labs to cycle through with sprintf
outputFileName = fullfile(rootDir,'Preprocess','run%s_C_transmid.mat'); % run number to full this out
trlfile  = fullfile(rootDir,'Preprocess', 'run%s_trl.mat'); % this contains our trial related info for epoching - should exist (generated with megTriggers script)
behavfile = fullfile(datadir,'behavioural',[thisSubject.id '_MEGRTs.mat']); % grab our behavioural data + output from a3_megtriggers

%% print subject for logging

disp('this is subject:')
disp(thisSubject.id);

warning('you need to run this locally!')

% first we'll load and epoch it all, so we can do other stuff while we wait
goodRuns = [];
for runNum = 1:numel(thisSubject.meg_runs)
    
    thisFile{runNum} = sprintf(outputFileName,runNum);
    
    if ~exist(thisFile{runNum},'file') || overwrite
        
        fprintf('loading and epoching %.0f of %.0f run files\n',runNum,numel(thisSubject.meg_runs))
        
        %% set up a fieldtrip config
        cfg = [];
        cfg.continuous = 'yes'; % this data is not epoched, it's continuous
        cfg.dataset = sprintf(inputFileName,thisSubject.meg_labs{runNum});
        
        %% we can look at the data
        % e.g. to see if we need to filter harmonics of line noise (see filtering section below)
        % or do manual artefact detection at the outset if you know something
        % weird happened
        %     cfg.channel = {'MEG'}; cfg.ylim = [-5e-11 5e-11];% grab just these, so we aren't plagued by all the STI channels and triggers etc
        % %     cfg.channel = {'EEG'}; cfg.ylim = [-0.0005 0.0005];
        %     cfg.blocksize = 10; % how many seconds of data to show at once
        %     % choose one of
        %     cfg.viewmode = 'butterfly';
        % %     cfg.viewmode = 'vertical';
        %     cfg = ft_databrowser(cfg);
        
        
        %% NOTE: if you looked at your data, be sure to re-run the previous ft config section before moving on
        
        %% now we want to define the trials
        % I have already produced a trial definition in spm/osl trl file
        % format, so I will not do it as per fieldtrip tutorials, but instead
        % import that file in the fieldtrip format (which is very similar)
        cfg.trialfun = 'trlFromFile'; % this is my own trial function that simply reads in from the spm/osl trl file I made in a3_megTriggers
        % I just arbitrarily add these to cfg struct so the function has something to read to find the files it needs
        cfg.trlfile = sprintf(trlfile,num2str(runNum));
        cfg.behavfile = behavfile;
        cfg = ft_definetrial(cfg);
        
        %% now we can do some preliminary artefact removal
        % we want to remove atypical artefacts, prior to any ICA
        % this is (a) because we don't want to waste components and
        % (b) because ICA can introduce artefacts (see e.g.
        % https://ohba-analysis.github.io/osl-docs/matlab/osl_example_africa.html#12)
        % since my trials are quite short (1500ms), and I care about the
        % whole timecourse of the trial,
        % I'm just going to remove any weird trials, and not worry about trying to remove
        % artefactual parts of trials
        rawData{runNum} = ft_preprocessing(cfg); % load it
        goodRuns = [goodRuns, runNum];
        
    else
        disp([thisFile{runNum} ' exists, skipping'])
    end
end

% now loop through whatever we have left and visually reject
count = 0;
for goodIdx = goodRuns
    count = count+1;
    fprintf('artefact checking %.0f of %.0f run files\n',count,numel(goodRuns))
    
    % make an eeg layout
    cfg = [];
    cfg.elec = rawData{goodIdx}.elec;
    eeglayout = ft_prepare_layout(cfg);
    % and an meg layout
    cfg = [];
    cfg.grad = rawData{goodIdx}.grad;
    meglayout = ft_prepare_layout(cfg);
    % combine them
    cfg = [];
    combinedLayout = ft_appendlayout(cfg, eeglayout, meglayout);
    clear meglayout eeglayout
    
    cfg = [];
    % if we want to do it manually
    cfg.method = 'summary';
    % ICA is better on continuous data, so lets fill these with nans,
    % rather than removing them
    cfg.keepchannel = 'nan';
    cfg.keeptrial = 'nan';
    cfg.metric = 'zvalue';
    cfg.layout = combinedLayout;
    
    % if we remove channels, and don't fill them with nans, then this would
    % change the channel numbers (e.g. removing EEG004 would make
    % EEG005 number 4 in later analyses) so keep that in mind
    % also, if you subselect (just taking MEG and EEG for example) this
    % will automatically reject the ones you didn't select
    % so careful on the channels 1:3 because these are eog and ecg
    % channels---don't fill them with NaNs or ica won't be able to use them to remove eye and heart artefacts!
    % also, anything not selected here will be removed automatically
    cfg.channel = {'MEG' 'EEG' 'EOG' 'ECG'};
    cleanedData{goodIdx} = ft_rejectvisual(cfg, rawData{goodIdx});
    
    % you can apparently do this automatically, but honestly this has
    % been a mess to figure out---in particular, my recent update of the
    % fieldtrip package to jan 23 has made the visual selection better
    % (don't need to do eeg and meg seperately and append them later)
    % but removes mne after ft_badsegment/badchannel so I can't loop through
    % runs. Perhaps someone else can work out these, or
    % ft_artefact_zvalue variations (muscle, jump, etc)
    
end

%% save these
% now loop through whatever we have left and save
count = 0;
for goodIdx = goodRuns
    count = count+1;
    fprintf('saving %.0f of %.0f run files\n',count,numel(goodRuns))
    thisData = cleanedData{goodIdx};
    save(thisFile{goodIdx},'thisData','-v7.3');
    clear thisData
    
end

end
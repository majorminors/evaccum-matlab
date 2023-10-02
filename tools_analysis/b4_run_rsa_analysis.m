function b4_run_rsa_analysis(thisSubject,datadir,toolsdir,overwrite)

if ~exist('overwrite','var'); overwrite = 0; end

%%%%%%%%%%%%
%% set up %%
%%%%%%%%%%%%

% set up paths
addpath /hpc-software/matlab/cbu/
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/');

addpath(datadir);
addpath(genpath(fullfile(toolsdir,'lib')))

modeldir = fullfile(toolsdir,'rdms');
ftDir = fullfile(toolsdir,'..','..','..','Toolboxes','fieldtrip');
addpath(ftDir); ft_defaults;
addpath(genpath(fullfile(toolsdir,'lib')))
toolbox = fullfile(toolsdir,'..','..','..','Toolboxes','BFF_repo'); addpath(genpath(toolbox)); clear toolbox
tmp = pwd;
cosmoDir = fullfile(toolsdir,'..','..','..','Toolboxes','CoSMoMVPA');
cd(fullfile(cosmoDir,'mvpa'));
cosmo_set_path();
cd(tmp); clear tmp
rootDir = fullfile(datadir,thisSubject.id);addpath(rootDir);
preProcDir = fullfile(rootDir,'Preprocess');
%alternatively, add the mvpa and external directories, and their subdirectories to the path, using addpath

behavfile = fullfile(datadir,'behavioural',[thisSubject.id '_MEGRTs.mat']); % grab our behavioural data + output from a3_megtriggers
trlfilePattern = fullfile(preProcDir, '%s_trl.mat'); % this contains our trial related info for epoching - should exist (generated with megTriggers script)
inputFilePattern = fullfile(preProcDir,'run*_if_transmid.fif');

%%%%%%%%%%%%%%%%%%%%%%%%
%% load and prep data %%
%%%%%%%%%%%%%%%%%%%%%%%%


disp('this is subject:')
disp(thisSubject.id);

% so here, load what you care about and play with them!
disp('loading')

theseFiles = dir(sprintf(inputFilePattern,thisSubject.id));
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
disp('>>> appending datasets')
% append them all
cfg = [];
cfg.keepsampleinfo  = 'yes';
coherenceOnsetErpData = ft_appenddata(cfg,cohOnsetErpData{:}); clear cohOnsetErpData
responseLockedErpData = ft_appenddata(cfg,respLockedErpData{:}); clear respLockedErpData
coherenceOnsetTfrData = ft_appenddata(cfg,cohOnsetTfrData{:}); clear cohOnsetTfrData
responseLockedTfrData = ft_appenddata(cfg,respLockedTfrData{:}); clear respLockedTfrData

% now we'll get the ft datasets we want
cfg = [];
cfg.trials = getTrialIndicesConditions(responseLockedErpData,1); % good trials, 1st condition
ecer_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedErpData);

% prep for cosmo
% https://www.cosmomvpa.org/contents_demo.html#demo-meeg-timelock-searchlight

% first the conditions of interest
% get a template of zeroes for the number of trials
trialTemplate = zeros(size(coherenceOnsetErpData.trialinfo,1),0);

% anon function to get trial indices where trial is good and
% condition==condCode
ds.sa.conditions = trialTemplate;
getTrialIndicesConditions = @(x,condCode) find(x.trialinfo(:,2) == 1 & x.trialinfo(:,1) == condCode);
for i = 1:4 % since there are four conditions, loop though 1:4
    % use that to grab the trial indices
    x = getTrialIndicesConditions(coherenceOnsetErpData,i);
    % create a trial index of the condition codes
    ds.sa.conditions(x) = i;
end; clear x i

% anon function to get trial indices where trial is good and
% [coherence difficulty, categorisation difficulty] selected
getTrialIndicesManipulations = @(x,condCode) find(x.trialinfo(:,2) == 1 & (x.trialinfo(:,1) == condCode(1) | x.trialinfo(:,1) == condCode(2)));
ds.sa.coherence = trialTemplate;
ecoh = getTrialIndicesManipulations(coherenceOnsetErpData,[1,2]);
hcoh = getTrialIndicesManipulations(coherenceOnsetErpData,[3,4]);
ds.sa.coherence(ecoh) = 1;
ds.sa.coherence(hcoh) = 2;
ds.sa.categorisation = trialTemplate;
ecat = getTrialIndicesManipulations(coherenceOnsetErpData,[1,3]);
hcat = getTrialIndicesManipulations(coherenceOnsetErpData,[2,4]);
ds.sa.categorisation(ecat) = 1;
ds.sa.categorisation(hcat) = 2;
clear ecoh hcoh ecat hcat

ds = cosmo_meeg_dataset(coherenceOnsetErpData);

%     % lina did a pca
%     ds_pca=ds;
%     ds_pca.samples=[];
%     ds_pca.fa.time=[];
%     ds_pca.fa.chan=[];
%     disp('computing pca')
%     for t=1:length(ds.a.fdim.values{2})
%         if ~mod(t,3);fprintf('.');end;
%         maskidx = ds.fa.time==t;
%         dat = ds.samples(:,maskidx);
%         [~,datpca,~,~,explained] = pca(dat);
%         datretain = datpca(:,cumsum(explained)<=99);
%         ds_pca.samples = cat(2,ds_pca.samples,datretain);
%         ds_pca.fa.time = [ds_pca.fa.time t*ones(1,size(datretain,2))];
%         ds_pca.fa.chan = [ds_pca.fa.chan 1:size(datretain,2)];
%     end
%     ds_nopca=ds;

% make meg dsms and run correlation between models and meg dsms
%     dsa = ds_pca;
ds=cosmo_average_samples(ds, 'repeats',1, 'split_by', {'conditions'});
%     ds=cosmo_average_samples(ds, 'repeats',1, 'split_by', {'coherence'});
%     ds=cosmo_average_samples(ds, 'repeats',1, 'split_by', {'categorisation'});


% is this what we want? or do we want to do the conditions separately?
ds.sa.targets = (1:length(ds.sa.conditions))'; % treats each trial as one 'condition'
ds.sa.chunks = (1:length(ds.sa.conditions))'; % treats each trial as being independent

nbrhood=cosmo_interval_neighborhood(ds,'time','radius',1); % to calculate the correlations it looks at each timepoint separately

measure=@cosmo_target_dsm_corr_measure;
measure_args=struct();
measure_args.target_dsm=object_model;
measure_args.type='Spearman'; %correlation type between target and MEG dsms
measure_args.metric='Spearman'; %metric to use to compute MEG dsm
measure_args.center_data=true; %removes the mean pattern before correlating
%run searchlight
ds_rsm_binary=cosmo_searchlight(ds,nbrhood,measure,measure_args);






% matrix = csvread('output.csv');

% Stimulus (motion) representation fit through time
%   for easy coh
matrix = csvread(fullfile(modeldir,'rdm_stim_ec.csv'));
%   for hard coh
matrix = csvread(fullfile(modeldir,'rdm_stim_hc.csv'));
%   for easy cat
matrix = csvread(fullfile(modeldir,'rdm_stim_ed.csv'));
%   for hard cat
matrix = csvread(fullfile(modeldir,'rdm_stim_hd.csv'));

% Categorisation (decision boundary) representation fit through time
%   for easy coh
matrix = csvread(fullfile(modeldir,'rdm_decbdry_ec.csv'));
%   for hard coh
matrix = csvread(fullfile(modeldir,'rdm_decbdry_hc.csv'));
%   for easy cat
matrix = csvread(fullfile(modeldir,'rdm_decbdry_ed.csv'));
%   for hard cat
matrix = csvread(fullfile(modeldir,'rdm_decbdry_hd.csv'));

% Categorisation (?) (button press) representation fit through time
%   for easy coh
matrix = csvread(fullfile(modeldir,'rdm_resp_ec.csv'));
%   for hard coh
matrix = csvread(fullfile(modeldir,'rdm_resp_hc.csv'));
%   for easy cat
matrix = csvread(fullfile(modeldir,'rdm_resp_ed.csv'));
%   for hard cat
matrix = csvread(fullfile(modeldir,'rdm_resp_hd.csv'));

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

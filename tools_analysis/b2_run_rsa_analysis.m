function b2_run_rsa_analysis(thisSubject,datadir,toolsdir,overwrite)
%%%%%%%%%%%%
%% set up %%
%%%%%%%%%%%%

clear all %#ok

% set up paths
addpath /hpc-software/matlab/cbu/
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/');

rootDir = fullfile(datadir,thisSubject.id);
preProcDir = fullfile(rootDir,'Preprocess');
addpath(genpath(fullfile(toolsdir,'lib')))
toolboxdir = fullfile(toolsdir,'..','..','..','Toolboxes');
ftDir = fullfile(toolboxdir,'fieldtrip');
addpath(ftDir); ft_defaults;
toolbox = fullfile(toolboxdir,'BFF_repo'); addpath(genpath(toolbox)); clear toolbox
cosmoDir = fullfile(toolboxdir,'Toolboxes','CoSMoMVPA');cd(fullfile(cosmoDir,'mvpa'));cosmo_set_path();cd(rootdir)
%alternatively, add the mvpa and external directories, and their subdirectories to the path, using addpath

rdm = collectRdms(fullfile(toolsdir,'rdms'));

behavfile = fullfile(datadir,'behavioural',[thisSubject.id '_MEGRTs.mat']); % grab our behavioural data + output from a3_megtriggers
trlfilePattern = fullfile(preProcDir, '%s_trl.mat'); % this contains our trial related info for epoching - should exist (generated with megTriggers script)

% here we get the filtered (but not ICAed) data
inputFilePattern = fullfile(preProcDir,'run*_f_transmid.fif');

%%%%%%%%%%%%%%%
%% prep data %%
%%%%%%%%%%%%%%%

%% print subject for logging
disp('>>> this is subject:')
disp(thisSubject.id);

% set up some trial indexing functions we'll use for averaging and cosmo

% for conditions, where x is ft data structure and condCode is 1:4 for each condition
getTrialIndicesConditions = @(x,condCode) find(x.trialinfo(:,2) == 1 & x.trialinfo(:,1) == condCode);
% trialinfo(:,2), from trlFromFile() is a code for good vs bad trials based on the behavioural data coding
% trialinfo(:,1), is the condition codes 1:4---again from trlFromFile
% for manipulations were x is ft data and cond code is [coherence difficulty, categorisation difficulty]
getTrialIndicesManipulations = @(x,condCode) find(x.trialinfo(:,2) == 1 & (x.trialinfo(:,1) == condCode(1) | x.trialinfo(:,1) == condCode(2)));
% trialinfo(:,2), from trlFromFile() is a code for good vs bad trials based on the behavioural data coding
% trialinfo(:,1), is the condition codes 1:4---again from trlFromFile, so we'll pass a vector with two numbers to filter

% now we need to epoch the data
% this is more or less the same as b2_aggregate_data, but we want to do some slightly different stuff here (including keeping the trials)
theseFiles = dir(inputFilePattern);
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

    %% ok, now we need to epoch our data
    
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
    
end

clear rawData

%% append these now

disp('>>> appending datasets')

% append them all
cfg = [];
cfg.keepsampleinfo  = 'yes';
coherenceOnsetErpData = ft_appenddata(cfg,cohOnsetErpData{:}); clear cohOnsetErpData
responseLockedErpData = ft_appenddata(cfg,respLockedErpData{:}); clear respLockedErpData

%% compile averages broken into our conditions

% first response-locked

disp('>>> compiling condition-wise averages')

cfg = [];
cfg.keeptrials = 1; % we want to keep trials for cosmo
cfg.trials = getTrialIndicesConditions(responseLockedErpData,1); % good trials, 1st condition
data.ecer_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedErpData);
cfg.trials = getTrialIndicesConditions(responseLockedErpData,2); % good trials, 2nd condition
data.echr_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedErpData);
cfg.trials = getTrialIndicesConditions(responseLockedErpData,3); % good trials, 3rd condition
data.hcer_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedErpData);
cfg.trials = getTrialIndicesConditions(responseLockedErpData,4); % good trials, 4th condition
data.hchr_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedErpData);

% now onset locked

cfg = [];
cfg.keeptrials = 1; % we want to keep trials for cosmo
cfg.trials = getTrialIndicesConditions(coherenceOnsetErpData,1); % good trials, 1st condition
data.ecer_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetErpData);
cfg.trials = getTrialIndicesConditions(coherenceOnsetErpData,2); % good trials, 2nd condition
data.echr_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetErpData);
cfg.trials = getTrialIndicesConditions(coherenceOnsetErpData,3); % good trials, 3rd condition
data.hcer_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetErpData);
cfg.trials = getTrialIndicesConditions(coherenceOnsetErpData,4); % good trials, 4th condition
data.hchr_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetErpData);

%% now averaged across manipulations

% first response-locked


cfg = [];
cfg.keeptrials = 1; % we want to keep trials for cosmo
cfg.trials = getTrialIndicesManipulations(responseLockedErpData,[1,2]); % good trials, easy coherence
data.ec_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedErpData);
cfg.trials = getTrialIndicesManipulations(responseLockedErpData,[3,4]); % good trials, hard coherence
data.hc_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedErpData);
cfg.trials = getTrialIndicesManipulations(responseLockedErpData,[1,3]); % good trials, easy rule
data.er_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedErpData);
cfg.trials = getTrialIndicesManipulations(responseLockedErpData,[2,4]); % good trials, hard rule
data.hr_responseLockedAverage = ft_timelockanalysis(cfg, responseLockedErpData);

% and onset-locked

cfg = [];
cfg.trials = getTrialIndicesManipulations(coherenceOnsetErpData,[1,2]); % good trials, easy coherence
data.ec_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetErpData);
cfg.trials = getTrialIndicesManipulations(coherenceOnsetErpData,[3,4]); % good trials, hard coherence
data.hc_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetErpData);
cfg.trials = getTrialIndicesManipulations(coherenceOnsetErpData,[1,3]); % good trials, easy rule
data.er_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetErpData);
cfg.trials = getTrialIndicesManipulations(coherenceOnsetErpData,[2,4]); % good trials, hard rule
data.hr_coherenceLockedAverage = ft_timelockanalysis(cfg, coherenceOnsetErpData);

clear coherenceOnsetErpData
clear responseLockedErpData

disp('>>> done averaging')

disp('>>> now cosmo')

%% prep for cosmo
% https://www.cosmomvpa.org/contents_demo.html#demo-meeg-timelock-searchlight

for manipulation = {'ec' 'hc' 'ed' 'hd'}
    for timeLock = {'%s_responseLockedAverage' '%s_coherenceLockedAverage'}
        thisManipulation = manipulation{:};
        thisTimelock = sprintf(timeLock{:},thisManipulation);
        disp('begin rsm loop')
        thisStruct = data.(thisTimelock);
        
        % first the conditions of interest
        % get a template of zeroes for the number of trials
        ds.sa.condition = zeros(size(thisStruct.trialinfo,1),0);

        % I don't think I want to do this. I think maybe my timelock is already only the trials I care about?
        % if not, I don't think this would do what I wanted---I'd be *comparing* things, not selecting things
        % Instead, I can slice the structure to remove the trials I don't need with cosmo_slice
        switch thisManipulation
           case 'ec'
               ds.sa.condition(getTrialIndicesManipulations(thisStruct,[1,2])) = 1;
           case 'hc'
               ds.sa.condition(getTrialIndicesManipulations(thisStruct,[3,4])) = 1;
           case 'ed'
               ds.sa.condition(getTrialIndicesManipulations(thisStruct,[1,3])) = 1;
           case 'hd'
               ds.sa.condition(getTrialIndicesManipulations(thisStruct,[2,4])) = 1;
        end
        
        ds = cosmo_meeg_dataset(thisStruct);
             
        % % lina did a pca
        % ds_pca=ds;
        % ds_pca.samples=[];
        % ds_pca.fa.time=[];
        % ds_pca.fa.chan=[];
        % disp('computing pca')
        % for t=1:length(ds.a.fdim.values{2})
        %     if ~mod(t,3);fprintf('.');end;
        %     maskidx = ds.fa.time==t;
        %     dat = ds.samples(:,maskidx);
        %     [~,datpca,~,~,explained] = pca(dat);
        %     datretain = datpca(:,cumsum(explained)<=99);
        %     ds_pca.samples = cat(2,ds_pca.samples,datretain);
        %     ds_pca.fa.time = [ds_pca.fa.time t*ones(1,size(datretain,2))];
        %     ds_pca.fa.chan = [ds_pca.fa.chan 1:size(datretain,2)];
        % end
        % ds_nopca=ds;

        % ok, now we want to average across samples, to create pseudotrials. this reduces the noise.
        % Cat has some simulations that help you work out what kind of averaging is worth doing
        % so we should do this first
        ds = cosmo_average_samples(ds, 'repeats',1, 'split_by', {'condition'});

        % then we pick our targets
        ds.sa.targets = ds.sa.condition;
        % so we could do chunks as runs, if we kept those, but
        % different number for different participants
        % we can instead use the cosmo_chunkize to balance chunks
        % across targets
        ds.sa.chunks = cosmo_chunkize(ds,2); % 2nd input arg = num chunks

        % to calculate the correlations it looks at each timepoint separately 
        nbrhood=cosmo_interval_neighborhood(ds,'time','radius',1); 

        % select our measurement
        % I think we also want to normalise channels, for MEG. this is a measure_args property
        measure=@cosmo_target_dsm_corr_measure;
        measure_args=struct();
        measure_args.type='Spearman'; %correlation type between target and MEG dsms
        measure_args.metric='Spearman'; %metric to use to compute MEG dsm
        measure_args.center_data=true; %removes the mean pattern before correlating
        % run searchlight between our RDM and our trial data
        for model = fieldnames(rdm)'
            thisModel = model(:);
            measure_args.target_dsm = rdm.(thisModel).(thisManipulation);
            rsa.(thisTimelock).(thisModel) = cosmo_searchlight(ds,nbrhood,measure,measure_args);
        end % model loop
    end %timelock loop
end %manipulation loop
    
disp('>>> cosmo complete')

return
end % end function

%%%%%%%%%%%%%%%%%%
%% subfunctions %%
%%%%%%%%%%%%%%%%%%

function rdm = collectRdms(modeldir)
% Stimulus (motion) representation fit through time
%   for easy coh
rdm.stim.ec = csvread(fullfile(modeldir,'rdm_stim_ec.csv'));
%   for hard coh
rdm.stim.hc = csvread(fullfile(modeldir,'rdm_stim_hc.csv'));
%   for easy cat
rdm.stim.ed = csvread(fullfile(modeldir,'rdm_stim_ed.csv'));
%   for hard cat
rdm.stim.hd = csvread(fullfile(modeldir,'rdm_stim_hd.csv'));

% Categorisation (decision boundary) representation fit through time
%   for easy coh
rdm.decbdry.ec = csvread(fullfile(modeldir,'rdm_decbdry_ec.csv'));
%   for hard coh
rdm.decbdry.hc = csvread(fullfile(modeldir,'rdm_decbdry_hc.csv'));
%   for easy cat
rdm.decbdry.ed = csvread(fullfile(modeldir,'rdm_decbdry_ed.csv'));
%   for hard cat
rdm.decbdry.hd = csvread(fullfile(modeldir,'rdm_decbdry_hd.csv'));

% Categorisation (?) (button press) representation fit through time
%   for easy coh
rdm.resp.ec = csvread(fullfile(modeldir,'rdm_resp_ec.csv'));
%   for hard coh
rdm.resp.hc = csvread(fullfile(modeldir,'rdm_resp_hc.csv'));
%   for easy cat
rdm.resp.ed = csvread(fullfile(modeldir,'rdm_resp_ed.csv'));
%   for hard cat
rdm.resp.hd = csvread(fullfile(modeldir,'rdm_resp_hd.csv'));
end

function convertedData = convertFif(cfg, ftDir)
% fieldtrip clears the mne toolbox after first use? so let's force this in
% a distinct environment in case fieldtrip does this for a good reason

addpath(ftDir); ft_defaults; addpath(fullfile(ftDir, 'external','mne'));
convertedData = ft_preprocessing(cfg); % load it
            
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


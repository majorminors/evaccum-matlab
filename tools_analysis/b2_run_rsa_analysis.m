function b2_run_rsa_analysis(thisSubject,datadir,toolsdir,overwrite)
%%%%%%%%%%%%
%% set up %%
%%%%%%%%%%%%

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
cosmoDir = fullfile(toolboxdir,'CoSMoMVPA');cd(fullfile(cosmoDir,'mvpa'));cosmo_set_path();cd(toolsdir)
%alternatively, add the mvpa and external directories, and their subdirectories to the path, using addpath

behavfile = fullfile(datadir,'behavioural',[thisSubject.id '_MEGRTs.mat']); % grab our behavioural data + output from a3_megtriggers
trlfilePattern = fullfile(preProcDir, '%s_trl.mat'); % this contains our trial related info for epoching - should exist (generated with megTriggers script)

saveFileName = [rootDir filesep 'rsa.mat'];

% here we get the filtered (but not ICAed) data
inputFilePattern = fullfile(preProcDir,'run*_f_transmid.fif');

templateRdms = {'stim' 'decbdry' 'dec_simple' 'resp'}; % 'dec_detail_pred' 'dec_detail_null' 

%%%%%%%%%%%%%%%
%% prep data %%
%%%%%%%%%%%%%%%

%% print subject for logging
disp('>>> this is subject:')
disp(thisSubject.id);

if exist(saveFileName,'file')
    if ~overwrite
        warning('file exists & overwrite off: skipping participant')
        return
    elseif overwrite
        warning('file exists & overwrite on: removing file')
        system(['rm -rf ' saveFileName])
    end
end

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
%     % then baseline correction
%     baselineWindow = [-0.2 0];
%     cfg = [];
%     cfg.latency = baselineWindow;
%     baselineData = ft_selectdata(cfg,cohOnsetErpData{fileNum}); % hold onto this data for response locked
%     cfg = [];
%     cfg.demean          = 'yes';
%     cfg.baselinewindow  = baselineWindow; % assumes we're demeaning coherence-locked stimuli
%     cohOnsetErpData{fileNum} = ft_preprocessing(cfg,cohOnsetErpData{fileNum});
    % reject artefacts
    cohOnsetErpData{fileNum} = rejectArtefacts(cohOnsetErpData{fileNum}, combinedLayout,ftDir);
%     any(cellfun(@(x) any(isnan(x(:))), cohOnsetErpData{1}.trial))

    % downsample
    cfg = [];
    cfg.method = 'downsample';
    cfg.resamplefs = 250; % frequency at which the data will be resampled (default = 256 Hz) must be at least double your highest frequency
    cohOnsetErpData{fileNum} = ft_resampledata(cfg, cohOnsetErpData{fileNum});

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
    respLockedErpData{fileNum} = ft_redefinetrial(cfg,rawData);
    % we can calculate our own baseline correction here
    % so first get the mean of the baseline window data
%     baselineData = cellfun(@(x) mean(x,2), baselineData.trial,'UniformOutput',false);
%     for tidx = 1:length(baselineData)
%         % then subtract each trial baseline mean from the response-locked
%         % epoch data
%         respLockedErpData{fileNum}.trial{tidx} = respLockedErpData{fileNum}.trial{tidx}-baselineData{tidx};
%     end
    % or we could try to carry over the baseline correction from the onset
    % locked data, but this produces a bunch of nans for reasons I couldn't
    % figure out in a cursory investigation
%     respLockedErpData{fileNum} = ft_redefinetrial(cfg,cohOnsetErpData{fileNum});
    % reject artefacts
    respLockedErpData{fileNum} = rejectArtefacts(respLockedErpData{fileNum}, combinedLayout,ftDir);
%     any(cellfun(@(x) any(isnan(x(:))), respLockedErpData{1}.trial))

    % downsample
    cfg = [];
    cfg.method = 'downsample';
    cfg.resamplefs = 250; % frequency at which the data will be resampled (default = 256 Hz) must be at least double your highest frequency
    respLockedErpData{fileNum} = ft_resampledata(cfg, respLockedErpData{fileNum});
    
    clear baselineWindow baselineData
end

clear rawData

%% append these now

disp('>>> appending datasets')

% append them all
cfg = [];
cfg.keepsampleinfo='no';
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
cfg.keeptrials = 1; % we want to keep trials for cosmo
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

for manipulation = {'ec' 'hc' 'er' 'hr' 'ecer' 'echr' 'hcer' 'hchr'}
    for timeLock = {'%s_responseLockedAverage' '%s_coherenceLockedAverage'}
        thisManipulation = manipulation{:};
        thisTimelock = sprintf(timeLock{:},thisManipulation);
        fprintf('doing rsa for %s\n', thisTimelock)
        thisStruct = data.(thisTimelock);
        
        % get trialwise nans, because cosmo can't deal with these: 
        %   so get logical array of nans
        %   then see where all across dim 2 (channels) are nans
        %   and where all across dim 3 (timepoints) are nans
        %   so you have now a vector where each trial is all nans across
        %   channels and timepoints
%         trials_with_nans_ind = find(all(all(isnan(thisStruct.trial), 2), 3));
        trials_with_nans_log = all(all(isnan(thisStruct.trial), 2), 3);
%         % if you want to check channelwise nans:
        channels_with_nans = find(all(all(isnan(thisStruct.trial), 1), 3));
%         % if you want to check timepointwise nans:
        timepoints_with_nans = find(all(all(isnan(thisStruct.trial), 1), 2));

        % make a cosmo dataset
        ds = cosmo_meeg_dataset(thisStruct);
        % slice out the nan-trials
        ds = cosmo_slice(ds,~trials_with_nans_log);
        
        % grab the trial ids to construct our RDMs
        theseTrialIds = ds.sa.trialinfo(:,3);
             
        % % lina did a pca, but only once because it can reduce noise but
        % she thinks it's not necessary---Tijl apparently has a tutorial
        % somewhere that shows in LDA that it makes no difference if you
        % apply PCA or not
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

        % ok, now if we were decoding we'd want to average across samples, to create pseudotrials. this reduces the noise.
        % Cat has some simulations that help you work out what kind of averaging is worth doing
        % howeverm Catriona and Lina didn't do any trial averaging before
        % RSA. for Lina, in the RDM were averages of e.g. decoding
        % accuracies of each pair, but she didn't do any trial averaging
        % before running any decoding
%         ds = cosmo_average_samples(ds, 'repeats',1, 'split_by', {'condition'});

%         % then we pick our targets
%         ds.sa.targets = ds.sa.condition;
%         % so we could do chunks as runs, if we kept those, but
%         % different number for different participants
%         % we can instead use the cosmo_chunkize to balance chunks
%         % across targets
%         ds.sa.chunks = cosmo_chunkize(ds,2); % 2nd input arg = num chunks
        % or we could do what Lina did and treat each trial independently
        ds.sa.targets = (1:size(ds.samples,1))'; % each trial is one condition
        ds.sa.chunks = (1:size(ds.samples,1))'; % each trial is independent
        % but if we do this, we need to work out how to reformat our model
        % dsm to correspond to the pairings
        % i think we need something like a lookup table? use the
        % trlFromFile function to pull all the relevant info, then on the
        % models have something that allows us to index into the pairings?
        % and construct an rsm from that
        
        % to calculate the correlations it looks at each timepoint separately 
        nbrhood=cosmo_interval_neighborhood(ds,'time','radius',1); 

        % select our measurement
        % I think we also want to normalise channels, for MEG. this is a measure_args property
        measure=@cosmo_target_dsm_corr_measure;
        measure_args=struct();
        measure_args.type='Spearman'; % correlation type between target and MEG dsms
        measure_args.metric='Spearman'; % metric to use to compute MEG dsm
        measure_args.center_data=true; % removes the mean pattern before correlating
        % run searchlight between our RDM and our trial data to get the raw correlations
        if numel(templateRdms) > 1;comparisonRdms = {};end % init this, in case we want to check partials
        for thisModel = templateRdms
            fprintf('getting raw correlations for model %s\n', thisModel{:})
            thisRdm = makeRdm(thisModel{:},theseTrialIds);
            measure_args.target_dsm = thisRdm;
            rsa.(thisTimelock).(thisModel{:}) = cosmo_searchlight(ds,nbrhood,measure,measure_args);
            % fisher transform the correlation values so they are more
            % normally distributed for group analysis
            rsa.(thisTimelock).(thisModel{:}).fisher_transformed_samples=atanh(rsa.(thisTimelock).(thisModel{:}).samples);
            rsa.(thisTimelock).(thisModel{:}).rdm = thisRdm;
            if numel(templateRdms) > 1; comparisonRdms = [comparisonRdms {thisRdm}];end % collate these in case we want to check partials
            clear thisRdm
        end; clear thisModel % model loop
        if numel(templateRdms) > 1
            % ok, now we want to do the same thing, but to get partial correlations out
            % we could, in theory, do this in the same loop as before, but for some reason I was getting very different results when I did it, so rather than work out why I will just do it seperately
            measure_args=struct();
            measure_args.type='Spearman'; % correlation type between target and MEG dsms
            measure_args.metric='Spearman'; % metric to use to compute MEG dsm
            measure_args.center_data=true; % removes the mean pattern before correlating
            for thisModel = templateRdms
                fprintf('getting partial correlations for model %s\n', thisModel{:})
                thisRdm = makeRdm(thisModel{:},theseTrialIds);
                measure_args.target_dsm = thisRdm;
                compareFunc = @(x) ~isequal(x, thisRdm);
                toCompare = comparisonRdms(cellfun(compareFunc, comparisonRdms));
                measure_args.regress_dsm = toCompare; clear toCompare
                tmp = cosmo_searchlight(ds,nbrhood,measure,measure_args);
                rsa.(thisTimelock).(thisModel{:}).partial_samples = tmp.samples;
                rsa.(thisTimelock).(thisModel{:}).fisher_transformed_partial_samples = atanh(tmp.samples);
                clear tmp compareFunc toCompare thisRdm
            end % model loop
        end % if multiple models
    end %timelock loop
end %manipulation loop
    
disp('>>> cosmo complete')

disp('saving')
save(saveFileName,'rsa')
disp('done saving')

return
end % end function

%%%%%%%%%%%%%%%%%%
%% subfunctions %%
%%%%%%%%%%%%%%%%%%


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


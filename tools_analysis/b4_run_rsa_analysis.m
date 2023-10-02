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
addpath(genpath(fullfile(toolsdir,'lib')))
modeldir = fullfile(toolsdir,'rdms');
ftDir = fullfile(toolsdir,'..','..','..','Toolboxes','fieldtrip');
addpath(ftDir); ft_defaults;
addpath(genpath(fullfile(toolsdir,'lib')))
jobdir = fullfile(rootdir,'job_logging');
toolbox = fullfile(rootdir,'..','..','Toolboxes','BFF_repo'); addpath(genpath(toolbox)); clear toolbox
cosmoDir = fullfile(toolsdir,'..','..','..','Toolboxes','CoSMoMVPA');cd(fullfile(cosmoDir,'mvpa'));cosmo_set_path();cd(rootdir)
%alternatively, add the mvpa and external directories, and their subdirectories to the path, using addpath

inputFileName = ['Preprocess' filesep 'epoched_data.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%
%% load and prep data %%
%%%%%%%%%%%%%%%%%%%%%%%%

% since our folders are all e.g. S01, S12 etc, let's load information about
% the directories that have that pattern
tmp = dir([datadir filesep 'S*']);
subjectFolders = tmp([tmp.isdir]); clear tmp
for subjectNum = 1:numel(subjectFolders)
    % full path to the inputfile
    thisFile = fullfile(subjectFolders(subjectNum).folder,subjectFolders(subjectNum).name,inputFileName);
    % skip this loop if we can't find one
    if ~exist(thisFile,'file'); continue; end
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
    
    whichVars = {...
        % what ft preprocessed data do you want?
        'coherenceOnsetErpData','responseLockedErpData','coherenceOnsetTfrData','responseLockedTfrData',...
        };
    
    data{subjectNum} = load(thisFile,whichVars{:});

    % prep for cosmo
    % https://www.cosmomvpa.org/contents_demo.html#demo-meeg-timelock-searchlight
    
    % first the conditions of interest
    % get a template of zeroes for the number of trials
    trialTemplate = zeros(size(data{subjectNum}.coherenceOnsetErpData.trialinfo,1),0);
    
    % anon function to get trial indices where trial is good and
    % condition==condCode
    ds.sa.conditions = trialTemplate;
    getTrialIndicesConditions = @(x,condCode) find(x.trialinfo(:,2) == 1 & x.trialinfo(:,1) == condCode);
    for i = 1:4 % since there are four conditions, loop though 1:4
        % use that to grab the trial indices
        x = getTrialIndicesConditions(data{subjectNum}.coherenceOnsetErpData,i);
        % create a trial index of the condition codes
        ds.sa.conditions(x) = i;
    end; clear x i
    
    % anon function to get trial indices where trial is good and
    % [coherence difficulty, categorisation difficulty] selected
    getTrialIndicesManipulations = @(x,condCode) find(x.trialinfo(:,2) == 1 & (x.trialinfo(:,1) == condCode(1) | x.trialinfo(:,1) == condCode(2)));
    ds.sa.coherence = trialTemplate;
    ecoh = getTrialIndicesManipulations(data{subjectNum}.coherenceOnsetErpData,[1,2]);
    hcoh = getTrialIndicesManipulations(data{subjectNum}.coherenceOnsetErpData,[3,4]);
    ds.sa.coherence(ecoh) = 1;
    ds.sa.coherence(hcoh) = 2;
    ds.sa.categorisation = trialTemplate;
    ecat = getTrialIndicesManipulations(data{subjectNum}.coherenceOnsetErpData,[1,3]);
    hcat = getTrialIndicesManipulations(data{subjectNum}.coherenceOnsetErpData,[2,4]);
    ds.sa.categorisation(ecat) = 1;
    ds.sa.categorisation(hcat) = 2;
    clear ecoh hcoh ecat hcat
    
    ds = cosmo_meeg_dataset(data{subjectNum}.coherenceOnsetErpData);

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

end; clear theseFiles thisFile subjectFolders subjectNum subjectCode index pathParts

disp('loading complete')




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

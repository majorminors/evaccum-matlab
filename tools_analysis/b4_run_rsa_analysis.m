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
ftDir = fullfile(toolsdir,'..','..','..','Toolboxes','fieldtrip');
addpath(ftDir); ft_defaults;
addpath(genpath(fullfile(toolsdir,'lib')))
jobdir = fullfile(rootdir,'job_logging');
toolbox = fullfile(rootdir,'..','..','Toolboxes','BFF_repo'); addpath(genpath(toolbox)); clear toolbox
cosmoDir = fullfile(toolsdir,'..','..','..','Toolboxes','CoSMoMVPA');cd(fullfile(cosmoDir,'mvpa'));cosmo_set_path();cd(rootdir)
%alternatively, add the mvpa and external directories, and their subdirectories to the path, using addpath

loadData = 1;
inputFileName = ['Preprocess' filesep 'timelocked_averages.mat'];
rdm = collectRdms(fullfile(toolsdir,'rdms'));


%%%%%%%%%%%%%%%%%%%%%%%%
%% load and prep data %%
%%%%%%%%%%%%%%%%%%%%%%%%

% since our folders are all e.g. S01, S12 etc, let's load information about
% the directories that have that pattern
tmp = dir([datadir filesep 'S*']);
subjectFolders = tmp([tmp.isdir]); clear tmp
for subjectNum = 1:numel(subjectFolders)
    if loadData
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
        
        erpManips = {...
            'ec_responseLockedAverage' 'hc_responseLockedAverage'...
            'er_responseLockedAverage' 'hr_responseLockedAverage'...
            'ec_coherenceLockedAverage' 'hc_coherenceLockedAverage'...
            'er_coherenceLockedAverage' 'hr_coherenceLockedAverage'...
            };
        
        
        whichVars = {...
            erpManips{:}...
            };
        
        data{subjectNum} = load(thisFile,whichVars{:});
    elseif ~loadData
%         [coherenceOnsetErpData,responseLockedErpData,coherenceOnsetTfrData,responseLockedTfrData]...
%             = epochData(thisSubject,datadir,toolsdir);
    end
    
    % prep for cosmo
    % https://www.cosmomvpa.org/contents_demo.html#demo-meeg-timelock-searchlight
    
    for manipulation = {'ec' 'hc' 'ed' 'hd'}
        for timeLock = {'%s_responseLockedAverage' '%s_coherenceLockedAverage'}
            thisManipulation = manipulation{:};
            thisTimelock = sprintf(timeLock{:},thisManipulation);
            disp('begin rsm loop')
            thisStruct = data{subjectNum}.(thisTimelock);
            % hc onset
            
            % first the conditions of interest
            % get a template of zeroes for the number of trials
            ds.sa.condition = zeros(size(thisStruct.trialinfo,1),0);
%             
%             % anon function to get trial indices where trial is good and
%             % condition==condCode
%             ds.sa.conditions = trialTemplate;
%             getTrialIndicesConditions = @(x,condCode) find(x.trialinfo(:,2) == 1 & x.trialinfo(:,1) == condCode);
%             for i = 1:4 % since there are four conditions, loop though 1:4
%                 % use that to grab the trial indices
%                 x = getTrialIndicesConditions(coherenceOnsetErpData,i);
%                 % create a trial index of the condition codes
%                 ds.sa.conditions(x) = i;
%             end; clear x i
%             
            % anon function to get trial indices where trial is good and
            % [coherence difficulty, categorisation difficulty] selected
            getTrialIndicesManipulations = @(x,condCode) find(x.trialinfo(:,2) == 1 & (x.trialinfo(:,1) == condCode(1) | x.trialinfo(:,1) == condCode(2)));
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
%             
%             %     % lina did a pca
%             %     ds_pca=ds;
%             %     ds_pca.samples=[];
%             %     ds_pca.fa.time=[];
%             %     ds_pca.fa.chan=[];
%             %     disp('computing pca')
%             %     for t=1:length(ds.a.fdim.values{2})
%             %         if ~mod(t,3);fprintf('.');end;
%             %         maskidx = ds.fa.time==t;
%             %         dat = ds.samples(:,maskidx);
%             %         [~,datpca,~,~,explained] = pca(dat);
%             %         datretain = datpca(:,cumsum(explained)<=99);
%             %         ds_pca.samples = cat(2,ds_pca.samples,datretain);
%             %         ds_pca.fa.time = [ds_pca.fa.time t*ones(1,size(datretain,2))];
%             %         ds_pca.fa.chan = [ds_pca.fa.chan 1:size(datretain,2)];
%             %     end
%             %     ds_nopca=ds;
%             
%             % make meg dsms and run correlation between models and meg dsms
%             %     dsa = ds_pca;
%
%             do we want to average across samples?
%             ds = cosmo_average_samples(ds, 'repeats',1, 'split_by', {'condition'});
%             
%             

% questions:
% do I want to average across samples? why/whynot?
% how do I reduce the thing I give to cosmo to the trials I care about?
% does it make sense to just run the rsa on the model (modified to include
% only trials of interes) and the trials of interest?

            ds.sa.targets = ds.sa.condition;
            % so we could do chunks as runs, if we kept those, but
            % different number for different participants
            % we can instead use the cosmo_chunkize to balance chunks
            % across targets
            ds.sa.chunks = cosmo_chunkize(ds,2); % 2nd input arg = num chunks
            
            nbrhood=cosmo_interval_neighborhood(ds,'time','radius',1); % to calculate the correlations it looks at each timepoint separately          
            measure=@cosmo_target_dsm_corr_measure;
            measure_args=struct();
            measure_args.type='Spearman'; %correlation type between target and MEG dsms
            measure_args.metric='Spearman'; %metric to use to compute MEG dsm
            measure_args.center_data=true; %removes the mean pattern before correlating
            % run searchlight
            for model = fieldnames(rdm)'
                thisModel = model(:);
                measure_args.target_dsm = rdm.(thisModel).(thisManipulation);
                rsa.(thisTimelock).(thisModel) = cosmo_searchlight(ds,nbrhood,measure,measure_args);
            end % model loop
        end %timelock loop
    end %manipulation loop
    
end; clear theseFiles thisFile subjectFolders subjectNum subjectCode index pathParts

disp('loading complete')




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

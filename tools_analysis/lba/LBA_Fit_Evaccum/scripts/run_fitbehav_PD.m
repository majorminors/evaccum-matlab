%% Fitting LBA to behavioural data
% Adapted from Alessandro Tomassini's RDK task
% DM Last Edit: Feb 2020
%
% path to saved file can be found in p.save_path

% to run go to matlab opened terminal and run:
% squeue -u username (dm06)

% do 100 optimisations

%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);
p = struct(); % for parameters
t = struct(); % for temp vars

%% set up variables

% required
rootdir = 'C:\Users\doria\Nextcloud\desiderata\desiderata\04 Research\05 Evidence Accumulation\01 EvAccum Code';%'\\cbsu\data\Group\Woolgar-Lab\projects\Dorian\EvAccum'; % root directory - used to inform directory mappings

% only required if not testing
datadir = fullfile(rootdir,'data\behav_pilot_2');
lbadatadir = fullfile(datadir,'lba_fit'); % expects to find your data here and will save results in a sub-folder here
p.data_name = 'prepped_data.mat'; % data file name

p.save_name = 'Model_%s.mat';
p.rng_seed = 17; % the rng seed number - fixed for reproducibility
t.local = 1; % run locally? Or 0 will use cbu scheduler
p.testing = 1; % if you want to use testing data, then switch to 1 and add the data folder to the path, else to 0. will save to pwd/test_results/
t.subject = 1; % if testing, which subject do you want to run?
t.test_data_name = 'lba_test_data.mat'; % name of your test data

% these are all the model variants we want to test - different combinations of free parameters
p.design_space={[1,3],[1,4],[1,3,4],[1,3,4,5],[1,2],[1,2,3],[1,2,4],[1,2,3,4],[1,2,3,4,5],[1,5],[1,3,5],[1,4,5],[1,2,5],[1,2,3,5],[1,2,4,5]};
settings.randiter  = 100; % random search iters before optimization
settings.nosession = 1; % optimization iterations - more equals less chance of ending up in a local minimum
settings.overwrite = 1;

% directory mappings
if ~p.testing
    addpath(genpath(fullfile(rootdir, 'tools_analysis'))); % add tools folder to path (includes LBA scripts)
    p.save_path = fullfile(lbadatadir, 'results');
    if ~exist(p.save_path,'dir')
        mkdir(p.save_path);
    end
end

% get the data
if p.testing
   warning('you are running in test mode; saving into current folder');
   p.save_path = fullfile(pwd,'test_results');
   if ~exist(p.save_path,'dir')
       mkdir(p.save_path);
   end
   
   t.alldata = load(t.test_data_name);
   t.data = t.alldata.d.subjects(t.subject,:);
else
    t.fileinfo = dir(fullfile(lbadatadir,p.data_name));
    t.datapath = fullfile(lbadatadir,t.fileinfo.name);
    
    % get the data
    t.alldata = load(t.datapath);
    t.data = t.alldata.d.subjects;
end

%% enter the data
data = t.data; % here's the data
%% Prepare for parallel processing
if ~p.testing && ~t.local
    fprintf('prepping for parallel processing %s\n', mfilename);
    % OUT OF DATE 28 FEB 2020
    S = cbu_scheduler('custom',{'compute',30,5,36000});
    jobdir = ''; error('fix this');
    if ~exist(jobdir,'dir')
        mkdir(jobdir);
    else
        eval(sprintf('!rm -r --interactive=never %s*',jobdir));
    end
    
    S.JobStorageLocation = jobdir;
    
    
    % up to date - tells the scheduler the directories needed for the analysis
    % (so will need to duplicate your own genpaths)
    error('be generous')
    dependencies_path = {  cant genpath - must add each individually?      
        };
end

funanon = @fun_fitbehav_LBA_PD;

%% create the jobs
fprintf('creating jobs %s\n', mfilename);

J = [];
ind = 0;

for imod = 1:length(p.design_space) % for the number of combinations of free parameters you specified earlier (model variants)
        p.numfile = num2str(imod);
        t.subj_file_name = sprintf(p.save_name,p.numfile);
        settings.savename = fullfile(p.save_path,t.subj_file_name);  
        settings.modfeat  = p.design_space{imod}; % saves the model (which free parameters)
        settings.rseed = p.rng_seed; % fixed for reproducibility
        settings.data = data;

        % then, depending on overwrite setting, add this job
        if ~exist(settings.savename,'file')||settings.overwrite
            ind = ind + 1;
            J(ind).task = funanon; % run the fitting function as the task
            J(ind).n_return_values = 0;
            J(ind).input_args = {settings};
            J(ind).depends_on = 0; % if this requires previous scripts to be run (i.e. 'ind' =1 needs be be done before 'ind' =2)
        end

end; clear imod ind;

%% submit the jobs
fprintf('submitting jobs from %s\n', mfilename);

if p.testing % test on one subject
    warning('you are testing locally on subject %1.0f',t.subject);
    for ijob = 1:length(J)
        fun_fitbehav_LBA_PD(J(ijob).input_args{1},data(t.subject))
    end; clear ijob;
elseif t.local && ~p.testing % run locally
    warning('you are running locally, on all participants');
    for ijob = 1:length(J)
        fun_fitbehav_LBA_PD(J(ijob).input_args{1},data)
    end; clear ijob;
else % submit the jobs to the cluster
    cbu_qsub(J, S, dependencies_path) % then submit
end

% [status, id]=debrief_cluster(S.JobStorageLocation);



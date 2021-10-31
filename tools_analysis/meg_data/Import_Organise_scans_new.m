addpath /hpc-software/matlab/cbu/
addpath(genpath('/group/woolgar-lab/projects/Dorian/evaccum/evaccum/matlab'))

%numbers refer to the PID number
subjs = [01 02 03 04];
fld_tar = '/group/woolgar-lab/projects/Dorian/evaccum/evaccum-matlab/data/meg_pilot_1/megdata/';

%decide whether to update (0 - e.g. after scanning new participants) or whether
%to start from scratch (1);
overwrite = 0;

% if not overwrite select only the new subjects
if ~overwrite
    exdir=arrayfun(@(x) exist([fld_tar,num2str(x),'/'],'dir'),subjs);
    subjs =subjs(exdir==0);
end

%% Prepare for parallel processing
fprintf('prepping %s for parallel processing\n', mfilename);

S = cbu_scheduler('custom',{'compute',30,20,345600});
if ~exist(jobdir,'dir')
    mkdir(jobdir);
else
    eval(sprintf('!rm -r --interactive=never %s*',jobdir));
end

S.JobStorageLocation = jobdir;


% tells the scheduler the directories needed for the analysis (no
% genpath apparently, so need to specify subdirectories too
dependencies_path = {
    jobdir
    toolsdir
    fullfile(toolsdir,'Templates') % this gets included in local runs with genpath
    p.save_path
    };

%% create the jobs
fprintf('creating jobs %s\n', mfilename);

J = [];
ind = 0;

for imod = 1:length(p.design_space) % for the number of combinations of free parameters you specified earlier (model variants)
    p.numfile = num2str(imod);
    t.model_file_name = sprintf(p.save_name,p.numfile);
    settings.savename = fullfile(p.save_path,t.model_file_name);
    settings.modfeat  = p.design_space{imod}; % saves the model (which free parameters)
    settings.rseed = p.rng_seed; % fixed for reproducibility
    if p.testing
        settings.data = data(t.subject);
    else
        settings.data = data; % note that this is fine for a small amount of data, but you will need to make it pull data files individually per job with bigger data sets or it will be very slow
    end
    
    sb = 1:numel(subjs)

       baseF{sb} = [fld_tar,num2str(subjs(sb)),'/'];

       ID = getMEGID(sprintf('DM_evaccumpilot_%s',num2str(subjs(sb))));
       tidyup_evaccum(baseF{sb},ID,overwrite);
    
    % then, depending on overwrite setting, add this job
    if ~exist(settings.savename,'file')||settings.overwrite
        ind = ind + 1;
        J(ind).task = @fun_fitbehav_LBA_PD; % run the fitting function as the task
        J(ind).n_return_values = 0;
        J(ind).input_args = {settings};
        J(ind).depends_on = 0; % if this requires previous scripts to be run (i.e. 'ind' =1 needs be be done before 'ind' =2)
    end
    
end; clear imod ind;

%% submit the jobs
fprintf('submitting jobs from %s\n', mfilename);

if t.local || p.testing % run locally
    warning('you are running locally\n');
    if p.testing
        warning('you are testing locally on subject %1.0f\n',t.subject);
    end
    for ijob = 1:length(J)
        fun_fitbehav_LBA_PD(J(ijob).input_args{1}); % specification of {1} is to make sure the settings are entered in the right format
    end; clear ijob;
else % submit the jobs to the cluster
    cbu_qsub(J, S, dependencies_path) % then submit
end
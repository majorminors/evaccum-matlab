
%% Set parameters

design_space={[1,3],[1,4],[1,3,4],[1,3,4,5],[1,2],[1,2,3],[1,2,4],[1,2,3,4],[1,2,3,4,5],[1,5],[1,3,5],[1,4,5],[1,2,5],[1,2,3,5],[1,2,4,5]};
%design_space={[1,5],[1,3,5],[1,4,5],[1,2,5],[1,2,3,5],[1,2,4,5]};
settings.randiter  = 100;%random search iters before optimization
settings.nosession = 25;%optimization iterations
settings.overwrite = 0;

savename = '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Model/fitBehav/results/Mod_%s_PD_LBA_rtc_1.mat';
droot = '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Behav/';


addpath(genpath('/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Preprocessing/'));
addpath(genpath('/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Model/'));
addpath(genpath('/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Behav/'));


%% Get sample IDs
% --------------------
sname= getnames(droot,[],'*_MEGtrg_RT.mat');
settings.sname = sname;

%NOTE: add missing subjects (faster than fit the model again for all the
%subjects) %%% COMMENT OUT if need to re-fit all the sample
% settings.sname = {sname{38:40}};

%% Prepare for parallel processing

clc

  
S = cbu_scheduler();
S.NumWorkers = 37;
S.SubmitArguments = '-l mem=20GB -l walltime=96:00:00';
S.JobStorageLocation = '/home/at07/matlab/jobs/LBA/';

dependencies_path = {
    '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Preprocessing/'
    '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Preprocessing/lib/'
    
    '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Behav/'
    '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Behav/MEGBehav/'
    '/imaging/at07/Matlab/Projects/CBU2016/RDKUnc/Model/'
    '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Model/fitBehav/'
    '/imaging/at07/Matlab/Projects/CBU2015/RDKUnc/Model/LBA/Templates/'
    
    };

funanon = @fun_fitbehav_LBA_PD;

%% create the jobs
J = [];
ind = 0;

for imod = 1:length(design_space)
        numfile = num2str(imod);
        settings.savename = sprintf(savename,numfile);  
        settings.modfeat  = design_space{imod};
        settings.rseed = 17;% for reproducibility

        if ~exist(settings.savename,'file')||settings.overwrite
            ind = ind + 1;
            J(ind).task = funanon; %#ok<*SAGROW>
            J(ind).n_return_values = 0;
            J(ind).input_args = {settings};
            J(ind).depends_on = 0;
        end

end

%% submit the jobs to the cluster
 !rm -r '/home/at07/matlab/jobs/LBA/*'
 
 cbu_qsub(J, S, dependencies_path)

% [status, id]=debrief_cluster(S.JobStorageLocation);


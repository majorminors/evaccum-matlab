%%
clear all
clear classes

close all
clc
addpath /hpc-software/matlab/cbu/
addpath(genpath('/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/'))

t.local = 1;

root    = '/group/woolgar-lab/projects/Dorian/EvAccum/';
droot   = fullfile(root, 'data/meg_pilot_1/megdata/');
infld   = fullfile(droot, '%s/MEEG/MaxfilterOutput/');
tarfld  = fullfile(droot, '%s/MEEG/Preprocess/');
behav   = '/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/behavioural/%s/%s_EvAccum.mat';
megrt   = '/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/behavioural/%s/%s_MEGRTs.mat';
ICA     = fullfile(droot, '%s/MEEG/ICAOutput/ICA%s.mat');
outpt   = fullfile(droot, '%s/MEEG/Preprocess/MSL_%s.mat');
addpath(droot);

jobdir = fullfile(droot,'scheduled_jobs'); % where you'll save any scheduled jobs (i.e. running on the scheduler)


% droot = '/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/megdata/';
% infld = '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/MEGData/%s/MEEG/MaxfilterOutput/';
% tarfld= '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/MEGData/%s/MEEG/Preprocess/';
% behav = '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Behav/MEGBehav/%s_MEGRT.mat';
% megrt = '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Behav/%s_MEGtrg_RT.mat';
% ICA   = '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/MEGData/%s/MEEG/ICAOutput/ICA%s_PD.mat';
% outpt = '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/MEGData/%s/MEEG/Preprocess/MSL_%s.mat';
% addpath(droot);


[dname,IDnum] = getnames(droot,7);
bname         = getnames('/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/behavioural/',10);
bname = unique({bname{:}});

%% parallel processing
locpar = 0;%1 local 0 cluster

%%
clc


filename    = 'Run%d_%s_trans.fif';
outfirst    = 'Run%d_%s_trans.mat';
outfilename =  '%s_trans.mat';

settings.fld = tarfld;

settings.freshstart = 1;%delete everything before starting again -not needed for 1st preprocessing run
settings.overwrite = 0;
settings.ICAoverwrite = 0;

dependencies_path = {
    '/neuro/meg_pd_1.2/'
    '/hpc-software/matlab/cbu/'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/meg_data/'
    '/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/megdata/'
    '/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/behavioural/'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/external'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/functions'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/functions/adminfunc'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/functions/guifunc'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/functions/javachatfunc'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/functions/miscfunc'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/functions/popfunc'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/functions/resources'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/functions/sigprocfunc'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/functions/statistics'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/functions/studyfunc'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/functions/timefreqfunc'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/plugins'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/plugins/ADJUST1.1.1'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/plugins/dipfit2.3'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/plugins/dipfit2.3/standard_BEM'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/plugins/dipfit2.3/standard_BEM/elec'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/plugins/dipfit2.3/standard_BEM/skin'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/plugins/dipfit2.3/standard_BESA'
    '/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/eeglab13_5_4b/plugins/firfilt1.6.1'
    jobdir
    };

% dependencies_path = {
%     '/neuro/meg_pd_1.2/'
%     '/hpc-software/matlab/cbu/'
%     '/imaging/at07/Matlab/CommonScripts'
%     '/imaging/at07/Matlab/CommonScripts/lib/'
%     '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Preprocessing/'
%      '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Preprocessing/lib/'
%     '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/MEGData/'
%     '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Behav/'
%     '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Behav/MEGBehav/'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/external'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/functions'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/functions/adminfunc'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/functions/guifunc'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/functions/javachatfunc'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/functions/miscfunc'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/functions/popfunc'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/functions/resources'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/functions/sigprocfunc'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/functions/statistics'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/functions/studyfunc'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/functions/timefreqfunc'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/plugins'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/plugins/ADJUST1.1.1'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/plugins/dipfit2.3'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/plugins/dipfit2.3/standard_BEM'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/plugins/dipfit2.3/standard_BEM/elec'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/plugins/dipfit2.3/standard_BEM/skin'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/plugins/dipfit2.3/standard_BESA'
%     '/imaging/at07/Matlab/Toolboxes/eeglab13_5_4b/plugins/firfilt1.6.1'
%     };


funanon = @fun_rdk_preproc_PD;


%%
J = [];
ind = 0;

for subi = 3%1:length(dname)
    
       %ID = getMEGID(sprintf('AT_RDK2_%s',IDnum{subi}));
       ID = getMEGID(sprintf('DM_evaccumpilot_%s',IDnum{subi}));
       
       settings.outfname= sprintf([tarfld,outfilename],dname{subi},dname{subi});
       settings.dname = dname{subi};
       %settings.ctr   = ID.ctr;
       settings.bEEG  = ID.bad_eeg;
       settings.bMEG  = ID.bad_meg;
       settings.behav = sprintf(behav,bname{subi});
       settings.svbeh = sprintf(megrt,bname{subi}); 
       settings.ICA   = sprintf(ICA,dname{subi},dname{subi}); 
  
       settings.infname = {};settings.outfirst={};

        for runi = 1:12
            
            if exist(sprintf([infld filename], dname{subi}, runi, dname{subi}),'file')
                settings.infname{end+1} = sprintf([infld filename], dname{subi}, runi, dname{subi});
                settings.outfirst{end+1}= sprintf([tarfld outfirst], dname{subi}, runi, dname{subi});
            end
        end       

        if ~exist(sprintf(outpt,dname{subi},dname{subi}),'file') || settings.overwrite
        %if ~exist(sprintf(outpt, dname{subi},ID.subj),'file')||settings.overwrite %don't repeat if exists
            ind = ind + 1;
            J(ind).task = funanon; %#ok<*SAGROW>
            J(ind).n_return_values = 0;
            J(ind).input_args = {settings};
            J(ind).depends_on = 0;
       end
    
    
end

if locpar
    % ParType = 0;  % Fun on Login machines (not generally advised!)
    % ParType = 1;   % Run maxfilter call on Compute machines using spmd (faster)
    ParType = 2;   % Run on multiple Compute machines using parfar (best, but less feedback if crashes)
    
    % open matlabpool if required
    % matlabpool close force CBU_Cluster
    if ParType
        if  isempty(gcp('nocreate'))%matlabpool('size')==0;
            P = cbupool(6);
            P.ResourceTemplate='-l nodes=^N^,mem=22GB,walltime=10:00:00';
            parpool(P);
        end
    end
    
    parfor subi = 1:length(dname)
        
        fun_rdk_preproc_PD(J(subi).input_args{1});
        
    end
    
    if ParType
        parpool close force CBU_Cluster
    end
    
    
    
elseif locpar == 0
    S = cbu_scheduler('custom',{'compute',30,20,345600});
    %S = cbu_scheduler();
    %S.NumWorkers = 37;
    %S.SubmitArguments = '-l mem=20GB -l walltime=96:00:00'; %20GB 96:00:00 more than 10GB recommended
    
    %if ~exist('/home/at07/matlab/jobs/rdkICA/','dir'); mkdir('/home/at07/matlab/jobs/rdkICA/');end
    %S.JobStorageLocation = '/home/at07/matlab/jobs/rdkICA/';
    if ~exist(jobdir,'dir')
        mkdir(jobdir);
    else
        eval(sprintf('!rm -r --interactive=never %s*',jobdir));
    end
    
    S.JobStorageLocation = jobdir;
    
    %%
    %!rm -r /home/at07/matlab/jobs/rdkICA/*
if t.local % run locally
    warning('you are running locally\n');
    for ijob = 1:length(J)
        fun_rdk_preproc_PD(J(ijob).input_args{1}); % specification of {1} is to make sure the settings are entered in the right format
    end; clear ijob;
else % submit the jobs to the cluster
    cbu_qsub(J, S, dependencies_path) % then submit
end
    
end
% old pipeline (meg_pilot_1; meg_pilot_2)
%1 - Import_Organise_scans
%2 - Run_maxfilterRDK_PD -> somehow it fails when ran through the cluster
%    but it works fine on parpool
%3 - ConvertT1_PD + manual re-alignment
%4 - Run_MEGtrg_local
%5 - Run_RDK_preproc_PD

%% new pipeline
clear all


rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';
scriptdir = fullfile(rootdir,'tools_analysis'); cd(scriptdir)
datadir = fullfile(rootdir,'data','meg_pilot_4'); addpath(datadir);
runLocal = 1;
runBehav = 0;
subjectRange = 1;%[34 35 36 37];%0;%[-1 7]; % 0 does all; array like [-1 4] does 4 to end
jobdir = fullfile(rootdir,'job_logging','megtriggers01');
% functionToRun = @a1_importAndOrganiseScans; additionalParams={datadir,0}; % datadir top level and overwrite (0|1)
% functionToRun = @a2_megTriggers; additionalParams={datadir,scriptdir,0}; % datadir top level, analysis tools/bin/libs
% functionToRun = @a3_maxFilter; additionalParams={datadir,scriptdir,1}; % datadir top level, analysis tools/bin/libs and overwrite (0|1)
functionToRun = @a4_preProc01_atypical_artefacts; additionalParams={datadir,scriptdir,1};
% functionToRun = @a4_preProc02_filtering; additionalParams={datadir,scriptdir,1}; % run locally if manual
% functionToRun = @a4_preProc03_ica; additionalParams={datadir,scriptdir,0}; % run locally if manual


allSubjects = importParticipants();
if ~subjectRange; subjectRange = 1:numel(allSubjects); end
isneg = @(val) val < 0 ;
if any(isneg(subjectRange)); subjectRange = subjectRange(~isneg(subjectRange)):1:numel(allSubjects); end

% make an object to hold cbuscheduler jobs
clear J;
J(1:numel(allSubjects))=deal(struct('task',[],'n_return_values',[], ...
    'input_args',[],'depends_on',[]));

for subjectidx = subjectRange
    disp('starting with subject: ');
    
    thisSubject = allSubjects{subjectidx}; % let's make this easier to reference
    disp(thisSubject.id)
    % let's skip participants who are not useable
    if ~thisSubject.usable
        disp('marked not usable - skipping')
       continue
    end
        
    if runLocal
        warning('running locally');
        functionToRun(thisSubject,additionalParams{:}); cd(scriptdir);
    else
        disp('compiling job for submission');
        J(subjectidx).task=functionToRun;
        J(subjectidx).n_return_values=0;
        J(subjectidx).input_args={thisSubject,additionalParams{:}};
        J(subjectidx).depends_on=0;
        J(subjectidx).AutoAttachFiles=false;
    end
    
    
end

    

if ~runLocal
    
    % if we're skipping some participants, we'll need to clear those from jobs
    % first because we prepopulated J
    clearArray = [];
    for i = 1:numel(J)
        if isempty(J(i).task)
            clearArray = [clearArray,i];
        end
    end
    J(clearArray) = [];
    fprintf('you are looking at %.0f participants\n',numel(J));
    
    disp('creating scheduler object');
    % create a scheduler object
    clear S;
    S = cbu_scheduler('custom',{'compute',40,4,28800}); % cutsom params: compute job, 46 workers, 12 GB per worker, 14400 secs = 4 hours
    %     S.SubmitArguments=[S.SubmitArguments ' --exclusive=user']; % when we were testing a possible memory issue
    
    if ~exist(jobdir,'dir')
        disp('making job directory');
        mkdir(jobdir);
    else
        disp('job directory detected---crashing');
        error('you should really learn how to just delete these');
        %(sprintf('!rm -r --interactive=never %s*',jobdir)); not working
    end
    
    S.JobStorageLocation = jobdir;
    % stop it traversing the directory structure of all the functions in my
    % scripts
    jobopts.AutoAttachFiles=false;
    
    disp('submitting jobs');
    % submit the jobs
    cbu_qsub(J,S,{rootdir,scriptdir},jobopts);
end

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
subjectRange = 0;%[34 35 36 37];%0;%[-1 7]; % 0 does all; array like [-1 4] does 4 to end
jobdir = fullfile(rootdir,'job_logging','rsa_01');
% functionToRun = -@a1_importAndOrganiseScans; additionalParams={datadir,0};
% functionToRun = @a2_megTriggers; additionalParams={datadir,scriptdir,0}; 
% functionToRun = @a3_maxFilter; additionalParams={datadir,scriptdir,1};
% functionToRun = @a4_preProc_filtering; additionalParams={datadir,scriptdir,1,1};
% functionToRun = @a4_preProc_ica; additionalParams={datadir,scriptdir,0,{'f' '.fif' 'i'},0}; % run locally if manual
% functionToRun = @a4_preProc_atypical_artefacts; additionalParams={datadir,scriptdir,1,{'if' '.fif' 'C'},0,0}; % run locally if manual
% b1_data_inspection is a standalone file for local exploration
% functionToRun = @b2_aggregate_data; additionalParams={datadir,scriptdir,0};
functionToRun = @b2_run_rsa_analysis; additionalParams={datadir,scriptdir,1};

allSubjects = importParticipants();
if ~subjectRange; subjectRange = 1:numel(allSubjects); end
isneg = @(val) val < 0 ;
if any(isneg(subjectRange)); subjectRange = subjectRange(~isneg(subjectRange)):1:numel(allSubjects); end

% make an object to hold cbuscheduler jobs
clear J;
J(1:numel(allSubjects))=deal(struct('task',[],'n_return_values',[], ...
    'input_args',[],'depends_on',[]));

for subjectidx = subjectRange
    
    thisSubject = allSubjects{subjectidx}; % let's make this easier to reference
    
    fprintf('%s (idx: %.0f)\n', thisSubject.id, subjectidx);
    
    % let's skip participants who are not useable
    if ~thisSubject.usable
        disp('xxxxx marked not usable - skipping xxxxx')
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
    S = cbu_scheduler('custom',{'compute',12,12,28800}); % cutsom params: compute job, 46 workers, 12 GB per worker, 14400 secs = 4 hours
    %     S.SubmitArguments=[S.SubmitArguments ' --exclusive=user']; % when we were testing a possible memory issue
    
    if ~exist(jobdir,'dir')
        disp('making job directory');
        mkdir(jobdir);
    else
        warning('job directory detected');
        response = input('Do you want to overwrite? (y/n): ','s');
        if strcmpi(response,'y') || strcmpi(response,'yes')
            disp('continuing...');
            system(['rm -rf ' jobdir]);
            if exist(jobdir,'dir'); error('I couldnt delete it :('); end
            disp('re-making job directory');
            mkdir(jobdir);
        else
            disp('exiting...');
            error('aborted by user');
            % Handle the case when the user chooses not to continue
        end
    end
    
    S.JobStorageLocation = jobdir;
    % stop it traversing the directory structure of all the functions in my
    % scripts
    jobopts.AutoAttachFiles=false;
    
    disp('submitting jobs');
    % submit the jobs
    cbu_qsub(J,S,{rootdir,scriptdir},jobopts);
end

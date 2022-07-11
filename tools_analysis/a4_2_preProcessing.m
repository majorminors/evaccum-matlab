function a4_2_preProcessing(thisSubject,datadir,toolsdir)

% our preprocessing will do the following
% 1. convert our fif files
% 2. run a line and bandpass filter on the data
% 5. re-reference if required
% 3. do an ICA if required
% 4. remove problematic components automatically or manually or both
% 5. re-reference (again!??) if required
% 6. downsample
% 7. extract trial information and epoch using that

addpath /hpc-software/matlab/cbu/
addpath(genpath('/neuro/meg_pd_1.2/'));       % FIFACCESS toolbox if use some meg_misc functions below
addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/'));
addpath(genpath(fullfile(toolsdir,'eeglab13_5_4b')))
addpath(genpath(fullfile(toolsdir,'meg_data')))
addpath(fullfile(toolsdir,'..','..','..','Toolboxes','osl','osl-core'))


rootDir = fullfile(datadir,thisSubject.id);
inputFolder   = fullfile(rootDir,'MaxfilterOutput'); % where is the input? expects maxfiltered data, but you could do this on data that has not been maxfiltered
outputFolder  = fullfile(rootDir,'Preprocess'); % where should we put the output?
behavDir   = fullfile(datadir,'behavioural'); % where is our behavioural data?
behavfile = fullfile(behavDir,[thisSubject.id '_Evaccum.mat']); % what is it called?
megrtfile = fullfile(behavDir,[thisSubject.id '_MEGRTs.mat']); % and since I have put my MEG triggers somewhere else, we will specify that too
addpath(rootDir);

doFilter = 1;
doICA = 1;
doReRef = 0;
overwrite = 1;

randomSeed = 7; % a random seed for ICA decomposition with oslafrica (use this same seed for reproducability)

inputFiles = {}; % init this
outputFiles= {}; % init this

%% first lets put together the filenames of our input and output data
for runNum = 1:numel(thisSubject.meg_runs)
    inputFiles{runNum} = sprintf([inputFolder,filesep,thisSubject.meg_labs{runNum},'_trans.fif']);
    outputFiles{runNum} = sprintf([outputFolder,filesep,thisSubject.meg_labs{runNum},'_trans.mat']);
end

%% if we haven't already, lets convert fif files
% we might have already converted these, but let's convert again and place
% them in the Preprocess folder
for runi = 1:numel(inputFiles)
    if ~exist(outputFiles{runi},'file') || overwrite
        S=[];
        S.dataset       = inputFiles{runi};
        S.outfile       = outputFiles{runi};
        S.save          = 0;
        S.reviewtrials  = 0;
        S.channels      = 'all';
        S.continuous    = 1;
        S.checkboundary = 0;
        
        spm_eeg_convert(S);
    end
end


%% Line and highpass filter
if doFilter
    inputFiles = outputFiles;
    % electric line 50Hz in the uk needs to be removed
    for runi = 1:length(inputFiles)
        % let's save these into a new folder
        [pathstr,name,ext] = fileparts(inputFiles{runi});
        % by default spm_eeg_filter saves with prefix 'f', and we filter twice
        % so we'll end up with 'ff', so we look for those or make them
        outputFile = sprintf('%s/ff%s%s',pathstr,name,ext);
        if ~exist(outputFile,'file') || overwrite
            
            %highpass
            S   = [];
            S.D = inputFiles{runi};
            S.band = 'high';
            S.freq = 0.5; % removes the drift
            D = spm_eeg_filter(S);
            
            %line filter
            S   = [];
            S.D = D;
            S.band = 'stop';
            S.freq = [49 51];% This defines the notch filter frequency range, i.e. around 50Hz
            D = spm_eeg_filter(S);
            
            
            outputFiles{runi} = D.fullfile;
            D.save
        else
            outputFiles{runi} = outputFile;
        end
    end
end

%% Set bad channels and Re-reference
if doReRef
    inputFiles = outputFiles;
    % re-referencing to the average
    for runi = 1:length(inputFiles)
        
        [pathstr,name,ext] = fileparts(outputFiles{runi});
        inputFile=sprintf('%s/M%s%s',pathstr,name,ext);
        
        if (~exist(inputFile,'file') || overwrite)
            
            %add our bad channels
            
            D = spm_eeg_load(inputFiles{runi});
            
            %remove old badchannels
            D=D.badchannels(D.badchannels,0);
            %identify eeg badchannels
            badch = [];
            if ~isempty(thisSubject.bad_eeg)
                for ib = 1:numel(thisSubject.bad_eeg)
                    if thisSubject.bad_eeg(ib)<10
                        badch(ib) = D.indchannel(['EEG00',num2str(thisSubject.bad_eeg(ib))]);
                    else
                        badch(ib) = D.indchannel(['EEG0',num2str(thisSubject.bad_eeg(ib))]);
                    end
                end
            end
            
            if ~isempty(thisSubject.bad_meg)
                for ib = 1:numel(thisSubject.bad_meg)
                    if thisSubject.bad_meg(ib)<1000
                        badch(end+1) = D.indchannel(['MEG0',num2str(thisSubject.bad_meg(ib))]);
                    else
                        badch(end+1) = D.indchannel(['MEG',num2str(thisSubject.bad_meg(ib))]);
                    end
                end
            end
            %D = units(D, D.indchantype('EEG'), 'uV');
            %update with new badchannels
            D=D.badchannels(badch,1);
            D.save;
            clear D
            
            % and re-reference
            
            D = reref(inputFiles{runi});
            outputFiles{runi} = D{2}.Dfname;
            
        else
            outputFiles{runi}=inputFile;
        end
    end
end

%% do ICA
if doICA
    inputFiles = outputFiles;
    for runi = 1:length(inputFiles)
        
        [pathstr,name,ext] = fileparts(inputFiles{runi});
        outputFile=sprintf('%s/ica%s%s',pathstr,name,ext);
        
        % we need to prepare an EEG layout for africa
        D = spm_eeg_load(inputFiles{runi});
        cfg = [];
        cfg.output = 'eeg_layout.lay';
        cfg.elec = D.sensors('EEG');
        eeg_layout = ft_prepare_layout(cfg);
        
        if ~exist(outputFile,'file') || overwrite || doICA == 2
            
            rng(randomSeed); % random seed for reproducing ICA decomposition
            
            osl_startup
            
            % https://ohba-analysis.github.io/osl-docs/matlab/osl_example_africa.html#!
            
            % if we don't add "D = osl_africa..." it will save it as whatever
            % the D was that was loaded!
            
            if doICA == 1
                
                D = spm_eeg_load(inputFiles{runi});
                D = D.copy(outputFile); % so it doesn't save over the previous step
                
                % first automatically remove the bad components
                % with no modality specified will look for MEG
                osl_africa(D,...
                    'do_ica',true,...
                    'used_maxfilter',true,...
                    'artefact_channels',{'EOG','ECG'}...
                    );
                % now will look for EEG
                osl_africa(D,...
                    'do_ica',true,...
                    'auto_artefact_chans_corr_thresh',.35,...
                    'eeg_layout',eeg_layout.cfg.output,...
                    'modality',{'EEG'},...
                    'used_maxfilter',true,...
                    'artefact_channels',{'EOG','ECG'}...
                    );
                % automatically does artefact identification (default), and removes them
                % this creates D.ica
                
                
            elseif doICA == 2
                
                D = spm_eeg_load(inputFiles{runi});
                
                % we then want to redo the artefact rejection to make changes to the assignment if we like
                osl_africa(D,...
                    'do_ident','manual',...
                    'used_maxfilter',true,...
                    'artefact_channels',{'EOG','ECG'}...
                    );
                % now will look for EEG
                osl_africa(D,...
                    'do_ident','manual',...
                    'eeg_layout',eeg_layout.cfg.output,...
                    'modality',{'EEG'},...
                    'used_maxfilter',true,...
                    'artefact_channels',{'EOG','ECG'}...
                    );
                
            end
            
            osl_shutdown
            
            outputFiles{runi} = D.fullfile;
            
        else
            outputFiles{runi} = outputFile;
        end
        
    end
end

%% downsample and epoch
inputFiles = outputFiles;

% we need the behavioural data from the MEGtriggers script
if ~exist(megrtfile,'file')
    error('Prepare behavioural files first!')
end

% trial info files needed for epoching (put together in megtriggers script)
trls = getnames(outputFolder,[],'run*_trl.mat');
trls = {trls{:,1}};

load(megrtfile,'megBehaviouralData'); % load behavioural file from MRGtrgRT script
numtrials = megBehaviouralData(:,end-1:end); % pull the last two columns of megBehaviouralData: what run the trial belonged to, and a unique number for the trial respectively
numtrials = str2double(numtrials); % convert to double

%Xfiles=preprocPipeline(inputFiles,trls,numtrials);

% %% Re-reference
% Xfiles = outputf{1};clear outf;

% [pathstr,name,ext] = fileparts(Xfiles);
% outf=sprintf('%s/M%s%s',pathstr,name,ext);
% %if ~exist(outf,'file') || overwrite
%
% D = reref(Xfiles);
%
%
%
% Xfiles = D{2}.Dfname{1};

return
end





%% ancillary functions


function D = reref(fname)

matlabbatch{1}.spm.meeg.preproc.prepare.D = {fname};
matlabbatch{1}.spm.meeg.preproc.prepare.task{1}.avref.fname = 'avref_montage.mat';
matlabbatch{2}.spm.meeg.preproc.montage.D = {fname};
matlabbatch{2}.spm.meeg.preproc.montage.mode.write.montspec.montage.montagefile(1) = cfg_dep('Prepare: Average reference montage', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','avrefname'));
matlabbatch{2}.spm.meeg.preproc.montage.mode.write.montspec.montage.keepothers = 1;
matlabbatch{2}.spm.meeg.preproc.montage.mode.write.blocksize = 655360;
matlabbatch{2}.spm.meeg.preproc.montage.mode.write.prefix = 'M';
D = spm_jobman('run',matlabbatch);
spm_unlink('avref_montage.mat')

return
end


function Xfiles=preprocPipeline(infname,trls,numtrials)

spm_jobman('initcfg');


%% epoch + merge

for runFile = 1:numel(infname)
    
    clear matlabbatch;
    
    %Epoch+Rearrange labels
    matlabbatch{1}.spm.meeg.preproc.epoch.D(1) = cellstr(infname{runFile});
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.trlfile = cellstr(trls{runFile});
    
    matlabbatch{1}.spm.meeg.preproc.epoch.bc = 0;
    matlabbatch{1}.spm.meeg.preproc.epoch.eventpadding = 0;
    matlabbatch{1}.spm.meeg.preproc.epoch.prefix = 'e';
    
    output = spm_jobman('run',matlabbatch);
    
    %% Add  info needed for MVPA
    D = spm_eeg_load(output{1}.Dfname{1});
    %info to store: [trial num, block num, here you can add anything you wish]
    thesetrials = numtrials(find(numtrials(:,1) == runFile),2); % find the trial numbers related to this runFile which is the index for each file - files are the Maxfilter outputs ...trans.fif for each run - so what trials for this run
    info = num2cell([thesetrials';runFile.*ones(1,numel(thesetrials))],1); % put that into 'info'
    D = trialtag(D, ':', info) ;save(D); % tack 'info' onto the trial data
    
    clear matlabbatch;
    matlabbatch{1}.spm.meeg.preproc.prepare.D(1) = cellstr(output{1}.Dfname);
    
    %matlabbatch{2}.spm.meeg.preproc.prepare.D(1) = cfg_dep('Epoching: Epoched Datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
    matlabbatch{1}.spm.meeg.preproc.prepare.task{1}.sortconditions.label = {
        'HcHr'
        'EcHr'
        'HcEr'
        'EcEr'
        }';
    
    matlabbatch{2}.spm.meeg.preproc.downsample.D(1) = cfg_dep('Prepare: Prepared Datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
    matlabbatch{2}.spm.meeg.preproc.downsample.fsample_new = 500;
    matlabbatch{2}.spm.meeg.preproc.downsample.method = 'resample';
    matlabbatch{2}.spm.meeg.preproc.downsample.prefix = 'd';
    
    
    output = spm_jobman('run',matlabbatch);
    
    
    
    
    %get ready for another round
    clear matlabbatch;
    outfiles{runFile} = output{1}.Dfname; %#ok
    
    %check that number of trials corresponds
    D = spm_eeg_load(outfiles{runFile}{1});
    trl = load(trls{runFile},'trl');
    if numel(D.conditions)~= size(trl.trl,1); error([D.fname,' discrepant number of trials!']);end
    clear trl D;
    
    
    
end


%% merge
S=[];

%S.D needs to be  a char array
S.D = char([outfiles{:}]);

S.recode = 'same';
S.prefix = 'merged_';
D = spm_eeg_merge(S);



Xfiles = D.fullfile;

return
end


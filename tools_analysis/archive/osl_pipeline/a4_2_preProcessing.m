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
addpath('/neuro/meg_pd_1.2/');       % FIFACCESS toolbox if use some meg_misc functions below
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/');
addpath(fullfile(toolsdir,'fieldtrip-20190410')); ft_defaults;
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

% each step adds prefixes, so we will have to search for those files if they exist
doFilter = 0; filterDone = 1; % we filter first
doReRef = 0; rerefDone = 1; % then reref
doICA = 0; icaDone = 1; % then do the ICA
doEpoch = 0; epochDone = 1; % then epoch
overwrite = 0;

randomSeed = 7; % a random seed for ICA decomposition with oslafrica (use this same seed for reproducability)

if all([doFilter; filterDone])==true; error('filterDone will cause problems if you are also doing a filter'); end
if all([doReRef; rerefDone])==true; error('rerefDone will cause problems if you are also doing rereferencing'); end
if all([doICA; icaDone])==true; error('icaDone will cause problems if you are also doing ica'); end
if all([doEpoch; epochDone])==true; error('epochDone will cause problems if you are also epoching'); end

inputFiles = {}; % init this
outputFiles= {}; % init this

%% first lets put together the filenames of our input and output data
for runNum = 1:numel(thisSubject.meg_runs)
    inputFiles{runNum} = sprintf([inputFolder,filesep,thisSubject.meg_labs{runNum},'_trans.fif']);
    outputFiles{runNum} = sprintf([outputFolder,filesep,thisSubject.meg_labs{runNum},'_trans.mat']);
end

%% if we haven't already, lets convert fif files
for runi = 1:numel(inputFiles)
    if ~exist(outputFiles{runi},'file') || overwrite
        fprintf('converting %.0f of %.0f raw files\n',runi,numel(inputFiles))
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

% now let's catch up to where we're up to with the preprocessing
filterPrefix = 'ff'; % by default spm_eeg_filter saves with prefix 'f', and we filter twice so we'll end up with 'ff'
rerefPrefix = 'M'; % we rereference with spm montage, so it will add an 'M' to the files
ICAPrefix = 'ica'; % we add this to the file ourselves
epochPrefix = 'e'; % we also choose this

if filterDone; outputFiles = addPrefix(outputFiles,filterPrefix);end
if rerefDone; outputFiles = addPrefix(outputFiles,rerefPrefix);end
if icaDone; outputFiles = addPrefix(outputFiles,ICAPrefix);end
if epochDone; outputFiles = addPrefix(outputFiles,epochPrefix);end

%% Line and highpass filter
if doFilter
    % old output files become the input files
    inputFiles = outputFiles;
    % and new output files will have the filter prefix appended to the
    outputFiles = addPrefix(inputFiles,filterPrefix);
    % electric line 50Hz in the uk needs to be removed
    for runi = 1:length(inputFiles)
        fprintf('filtering %.0f of %.0f raw files\n',runi,numel(inputFiles))
        if ~exist(outputFile{runi},'file') || overwrite
            
            %highpass
            S   = [];
            S.D = inputFiles{runi};
            S.band = 'high';
            S.freq = 0.5; % removes the drift
            D = spm_eeg_filter(S); % by default spm_eeg_filter saves with prefix 'f'
            
            %line filter
            S   = [];
            S.D = D;
            S.band = 'stop';
            S.freq = [49 51];% This defines the notch filter frequency range, i.e. around 50Hz
            D = spm_eeg_filter(S); % by default spm_eeg_filter saves with prefix 'f'
            
            
            outputFiles{runi} = D.fullfile;
            D.save
        end
    end
end

%% Set bad channels and Re-reference
%
% we rereference with spm montage, so it will add an 'M' to the files
if doReRef
    % old output files become the input files
    inputFiles = outputFiles;
    % and new output files will have the reref prefix appended to the input files
    outputFiles = addPrefix(inputFiles,rerefPrefix);
    % re-referencing to the average
    for runi = 1:length(inputFiles)
        
        fprintf('rereferencing %.0f of %.0f raw files\n',runi,numel(inputFiles))
        if (~exist(outputFile{runi},'file') || overwrite)
            
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
            
        end
    end
end


%% do ICA

if doICA
    % old output files become the input files
    inputFiles = outputFiles;
    % and new output files will have the ica prefix appended to the input files
    outputFiles = addPrefix(inputFiles,ICAPrefix);
    for runi = 1:length(inputFiles)
        
        fprintf('ica for %.0f of %.0f raw files\n',runi,numel(inputFiles))
        
        % we need to prepare an EEG layout for africa
        D = spm_eeg_load(inputFiles{runi});
        cfg = [];
        cfg.output = 'eeg_layout.lay';
        cfg.elec = D.sensors('EEG');
        eeg_layout = ft_prepare_layout(cfg);
        
        if ~exist(outputFiles{runi},'file') || overwrite || doICA == 2
            
            rng(randomSeed); % random seed for reproducing ICA decomposition
            
            osl_startup
            
            % https://ohba-analysis.github.io/osl-docs/matlab/osl_example_africa.html#!
            
            % it is important that bad epochs are excluded prior to running AFRICA
            % look for bad channels
            D = osl_detect_artefacts(D,'badchannels',true,'badtimes',false,'modalities',modalities);
            % then look for bad segments
            D = osl_detect_artefacts(D,'badchannels',false,'badtimes',true,'modalities',modalities);

            
            if doICA == 1
                
                % osl_africa usually doesn't autosave after running, but it
                % does for ICA, because it's time consuming
                % so we load our input files
                D = spm_eeg_load(inputFiles{runi});
                % then we'll make a new copy with a different name so it doesn't save over the previous step
                D = D.copy(outputFiles{runi}); 
             
                % now each time we call osl_africa, call it as:
                % D = osl_africa...
                % because of the autosave, the D in memory won't be the
                % same D saved to disk!
                
                % first automatically remove the bad components
                % with no modality specified will look for MEG
                D = osl_africa(D,...
                    'do_ica',true,...
                    'used_maxfilter',true,...
                    'artefact_channels',{'EOG','ECG'}...
                    );
                % now will look for EEG
                D = osl_africa(D,...
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
                
                % WIP: this does a manual identification, but I haven't
                % checked it saves properly and whatnot
                
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
            
        end
        
    end
end

%% epoch and merge

if doEpoch
    % old output files become the input files
    inputFiles = outputFiles;
    % and new output files will have the epoch prefix appended to the input files
    outputFiles = addPrefix(inputFiles,epochPrefix); % this is actually unnecessary
    
    % we need the behavioural data from the MEGtriggers script
    if ~exist(megrtfile,'file')
        error('Prepare behavioural files first!')
    end
    
    % trial info files needed for epoching (put together in megtriggers script)
    trls = getnames(outputFolder,[],'run*_trl.mat');
    trls = arrayfun(@(x)(strcat([outputFolder filesep],x)), trls);
    
    load(megrtfile,'megBehaviouralData'); % load behavioural file from MRGtrgRT script
    numtrials = megBehaviouralData(:,end-1:end); % pull the last two columns of megBehaviouralData: what run the trial belonged to, and a unique number for the trial respectively
    numtrials = str2double(numtrials); % convert to double
    
    spm_jobman('initcfg');
    
    for runFile = 1:length(inputFiles)
        
        clear matlabbatch;
        
        %Epoch+Rearrange labels
        matlabbatch{1}.spm.meeg.preproc.epoch.D(1) = cellstr(inputFiles{runFile});
        matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.trlfile = cellstr(trls{runFile});
        
        matlabbatch{1}.spm.meeg.preproc.epoch.bc = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.eventpadding = 0;
        matlabbatch{1}.spm.meeg.preproc.epoch.prefix = epochPrefix;
        
        output = spm_jobman('run',matlabbatch);
        
        %% Add  info needed for MVPA
        D = spm_eeg_load(output{1}.Dfname{1});
        %info to store: [trial num, block num, here you can add anything you wish]
        thesetrials = numtrials(find(numtrials(:,1) == runFile),2); % find the trial numbers related to this runFile which is the index for each file - files are the Maxfilter outputs ...trans.fif for each run - so what trials for this run
        if strcmp(thisSubject.id,'S15')
            if runFile == 2
                % don't tag those first missing trials
                thesetrials = thesetrials(size(thesetrials,1) - size(D,3) + 1:end);
                skipTrlNums = 1;
            end
        end
        info = num2cell([thesetrials';runFile.*ones(1,numel(thesetrials))],1); % put that into 'info'
        D = trialtag(D, ':', info) ;save(D); % tack 'info' onto the trial data
        
        %     clear matlabbatch;
        %     matlabbatch{1}.spm.meeg.preproc.prepare.D(1) = cellstr(output{1}.Dfname);
        %
        %     %matlabbatch{2}.spm.meeg.preproc.prepare.D(1) = cfg_dep('Epoching: Epoched Datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
        %     matlabbatch{1}.spm.meeg.preproc.prepare.task{1}.sortconditions.label = {
        %         'HcHr'
        %         'EcHr'
        %         'HcEr'
        %         'EcEr'
        %         }';
        %
        %     matlabbatch{2}.spm.meeg.preproc.downsample.D(1) = cfg_dep('Prepare: Prepared Datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
        %     matlabbatch{2}.spm.meeg.preproc.downsample.fsample_new = 500;
        %     matlabbatch{2}.spm.meeg.preproc.downsample.method = 'resample';
        %     matlabbatch{2}.spm.meeg.preproc.downsample.prefix = 'd';
        %
        %
        %     output = spm_jobman('run',matlabbatch);
        %
        
        
        
        %get ready for another round
        clear matlabbatch;
        outputFiles{runFile} = output{1}.Dfname; %#ok
        
        %check that number of trials corresponds
        D = spm_eeg_load(outputFiles{runFile}{1});
        trl = load(trls{runFile},'trl');
        if ~exist('skipTrlNums','var'); skipTrlNums = 0; end
        if numel(D.conditions)~= size(trl.trl,1) && ~skipTrlNums; error([D.fname,' discrepant number of trials!']);end
        clear trl D;
    end
    
    
    
end


%% merge
% make a temp folder to work in, because I can't figure out how to
% determine the name or save location of spm_eeg_merge
tmpFld = [outputFolder filesep 'tmpMerge'];
mkdir(tmpFld)
cd(tmpFld)

S=[];

% collect the files to be merged
if epochDone
    S.D = char(outputFiles{:});
else
    % S.D needs to be  a char array, and if we epoch we collect these from
    % output{1}.Dfname - result of matlab batch processing
    S.D = char([outputFiles{:}]);
end

% add options and perform the merge
S.recode = 'same';
S.prefix = 'merged_';
D = spm_eeg_merge(S); % will save in the current folder and apparently with the first name in S.D?

% so since itll be saved same filename as the first file, lets grab the
% 'Run1' part and change it to something information
% started work trying to make this a bit name agnostic, but gave up for
% time
% prefixes = regexp(D.fname,[S.prefix '(\w*)Run'],'tokens'); % grab all the prefixes - everything before 'Run'
% toReplace = regexp(D.fname,[prefixes{:}{:} '(\w*)_trans'],'tokens'); % grab whatever is between the prefixes and '_trans'
% newFilename = replace(D.fname,toReplace{:}{:},'allRuns'); % replace that with something more informative
% copyfile(D.fname,[outputFolder filesep newFilename],'f'); % copy the file from here to where it should go

system('mv -f * ..'); % move everyhting here one directory up

cd(toolsdir) % go back home
rmdir(tmpFld,'s') % remove that temp folder we were working in

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



function outFiles = addPrefix(files,prefix)

for fileIdx = 1:length(files)
    [pathstr,name,ext] = fileparts(files{fileIdx});
    outFiles{fileIdx} = sprintf(['%s/' prefix '%s%s'],pathstr,name,ext);
end

return
end


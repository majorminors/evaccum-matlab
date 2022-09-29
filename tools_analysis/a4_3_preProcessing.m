function a4_3_preProcessing(thisSubject,datadir,toolsdir)

addpath /hpc-software/matlab/cbu/
addpath('/neuro/meg_pd_1.2/');       % FIFACCESS toolbox if use some meg_misc functions below
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/');
addpath(fullfile(toolsdir,'fieldtrip-20190410')); ft_defaults;
% addpath(genpath(fullfile(toolsdir,'eeglab13_5_4b')))
addpath(genpath(fullfile(toolsdir,'meg_data')))
addpath(fullfile(toolsdir,'..','..','..','Toolboxes','osl','osl-core'))

rootDir = fullfile(datadir,thisSubject.id);
inputFolder   = fullfile(rootDir,'MaxfilterOutput'); % where is the input? expects maxfiltered data, but you could do this on data that has not been maxfiltered
outputFolder  = fullfile(rootDir,'Preprocess'); % where should we put the output?
behavDir   = fullfile(datadir,'behavioural'); % where is our behavioural data?
behavfile = fullfile(behavDir,[thisSubject.id '_Evaccum.mat']); % what is it called?
megrtfile = fullfile(behavDir,[thisSubject.id '_MEGRTs.mat']); % and since I have put my MEG triggers somewhere else, we will specify that too
trlfile  = fullfile(outputFolder, 'run%s_trl.mat'); % this contains our trial related info for epoching - should exist (generated with megTriggers script)
addpath(rootDir);


manualChecks = 1; % enable manual checking and selecting throughout
overwrite = 1; 

% if any of these are 1, then it will just do it, if 0 will assume it has been done and skip processing but still search for the appropriately prefixed file 
doConvert = 1; % convert to SPM M/EEG object?
doFilter = 1; % filter line noise and do a high pass filter?
doArtefacts = 1; % detect bad channels and trials?
doICA = 1; % run our ICA? if 2, will skip ica, and just do manual checks (if manualChecks is also on)
doEpoch = 1; % epoch the data?

randomSeed = 7; % a random seed for ICA decomposition with oslafrica (use this same seed for reproducability)
filterPrefix = 'f'; % most programs default to an f, but note we filter twice so will have a file with ff prepended
artefactPrefix = 'b'
icaPrefix = 'ica'; % we add this to the file ourselves
epochPrefix = 'e'; % we also choose this

%% convert and merge

inputFiles = {}; % init this

% first we'll convert and move our files to the preprocessing folder
% maxfilter in our pipeline has first run sss then transformed to a
% the coordinate frames of the middle run for easy averaging, so we will do
% this not on the sss file, but the trans file
for runNum = 1:numel(thisSubject.meg_runs)
    % input files from input folder
    tmpIn = sprintf([inputFolder,filesep,thisSubject.meg_labs{runNum},'_trans.fif']);
    % output files will go into output folder
    tmpOut = sprintf([outputFolder,filesep,thisSubject.meg_labs{runNum},'_trans.mat']);
    % if we haven't already, lets convert fif files to SPM M/EEG objects
    if ~exist(tmpOut,'file') || (overwrite && doConvert)
        fprintf('converting %.0f of %.0f raw files\n',runi,numel(thisSubject.meg_runs))
        S=[];
        S.dataset       = tmpIn;
        S.outfile       = tmpOut;
        S.save          = 0;
        S.reviewtrials  = 0;
        S.channels      = 'all';
        S.continuous    = 1;
        S.checkboundary = 0;
        
        spm_eeg_convert(S);
    end % doConvert if statement

    % collate these
    inputFiles{runNum} = tmpOut;

    % also, we need to collate the epoching trl information, which is currently in different files according to run
    thisTrlFile = sprintf(trlfile,num2str(runi));
    thisEpochInfo = load(thisTrlFile); % should have vars 'trl' and 'conditionlabels' - this is how both osl and spm want it
    epochInfo.trl = [epochInfo.trl;thisEpochInfo.trl];
    epochInfo.conditionlabels = [epochInfo.conditionlabels;thisEpochInfo.conditionlabels];

end; clear thisEpochInfo thisTrlFile % run loop

% now let's merge the files into one
% make a temp folder to work in, because I can't figure out how to
% determine the name or save location of spm_eeg_merge
tmpFld = [outputFolder filesep 'tmpMerge'];
mkdir(tmpFld)
cd(tmpFld)
S=[];
% collect the files to be merged
S.D = char(inputFiles{:}); % might be char([outputFiles{:}])
% add options and perform the merge
S.recode = 'same';
S.prefix = 'merged_';
D = spm_eeg_merge(S); % will save in the current folder and apparently with the first name in S.D?
% so let's move it and rename it something more informative
system(['for file in *; do mv -f "$file" "../allRuns_trans.${file##*.}"; done']); % bash loop that moves everything one directory up, and changes the filename (but keeps the extension)
cd(toolsdir) % go back home
rmdir(tmpFld,'s') % remove that temp folder we were working in
% now let's grab our newly merged file
inputFile = sprintf([outputFolder,filesep,'allRuns_trans.mat']);
clear tmpIn tmpOut tmpFld

% so we now have an SPM M/EEG object saved
% if we load:
% D = spm_eeg_load(inputFile) 
% we can view some information about D just by running it alone
% D
% D is actually an "object" and therefore uses functions called 'methods' we can use to query information about it
% the internal structure is not accessible directly when working with it
% D.ntrials % gives you the number of trials
% D.conditions % shows a list of condition labels per trial, this should be 'undefined' in continuous data, but will have meaningful names in epoched data
% D.condlist % shows the list of unique conditions, again this should be 'undefined' in continuous data, but will have meaningful names in epoched data
% D.chanlabels % order and names of channels, you should be able to see channel names corresponding to the MEG sensors ('MEG1234') alongsize EOG, ECG and trigger channels ('STIXXX')
% D.chantype % type of channel, such as eeg ('EEG'), magnetometers ('MEGMAG') and planar gradiometers ('MEGPLANAR')
% This will show the size of the data matrix:
% D.size % size of the data matrix in channels, samples, and trials
% respectively (continuous data will have only one trial, epoched data will
% reflect your number of epochs).
% we can also index into D with D(sensor_number,timepoint,trial), which
% will return the field strength of the sensor at the timepoint/trial
% % For the full list of methods performing operations with the object, try:
% methods('meeg')
% and for specific methods:
% help meeg/method_name
% you can also query with function syntax:
% ntrials(D)

% we are going to swap now to use OSL for some of the next bits of the pipeline
% though note that we'll use spm_eeg functions occasionally and can swap bits and
% pieces out for other toolboxes (e.g. fieldtrop/eeglab)

osl_startup

if manualChecks
    % we can now use OSLVIEW to look manually at continuous data (not epoched data)
    % we can also interactively flag bad channels and time periods with it
    D = spm_eeg_load(inputFile);
    D = oslview(D);
end

% you'd see there, at a minimum low frequency noise and electric line noise (50Hz in the UK plus its harmonics)

%% filtering
% so let's filter that out
% note that it's probably best to epoch AFTER filtering, so we don't end up with edge effect problems
if doFilter
    % we filter twice, so we add the filterPrefix twice because spm_eeg_filter saves with a prefix
    outputFile = addPrefix(inputFile,[filterPrefix filterPrefix]);
    if ~exist(outputFile,'file') || overwrite
        % we'll do a band pass to start, though note that removing noise between
        % the specified frequency range means we can't look at it!
        % if we want to do spm_eeg_filter
        S   = [];
        S.D = inputFile;
        S.band = 'high';
        S.freq = 0.5; % removes the drift
        S.prefix = filterPrefix; % f is the default
        D = spm_eeg_filter(S); % by default spm_eeg_filter saves with prefix 'f'
        % we can also use osl
%             D=spm_eeg_load(inputFiles{runi});
%             D=osl_filter(D,[0.1 120],'prefix',filterPrefix); % this filter range is just the OSL default! change for your usecase
        
        % now line noise, plus any harmonics, but remember that filtering is blind to the origin of your signal, it might remove
        % both line noise and any other neural signals at the same frequency, for example real gamma activity
        % good to double-check data after filtering to avoid surprises later. oslview will let you do this
%             D=spm_eeg_load(inputFiles{runi});
%             D=osl_filter(D,-[48 52],'prefix',''); % removes line noise using a notch filter
%             D=osl_filter(D,-[98 102],'prefix',filterPrefix); % removes harmonic of line noise using a notch filter
        % similarly in spm_eeg
        S   = [];
        S.D = D;
        S.band = 'stop';
        S.freq = [49 51];% This defines the notch filter frequency range, i.e. around 50Hz
        S.prefix = filterPrefix; % f is the default
        D = spm_eeg_filter(S); % by default spm_eeg_filter saves with prefix 'f'
    end % exist/overwrite statment
end % doFilter if statment
% update inputFile
inputFile = outputFile; % now filtering prefixes have been prepended


%% downsampling
% we can also downsample here, if we want to speed things up or save space
% we'd just use spm_eeg to do that like this:
% S.D=inputFile;
% S.fsample_new = 150; % in Hz
% D = spm_eeg_downsample (S);
% note that there's an interaction between downsampling and filtering though, so we might want to downsample earlier---maxfilter can also downsample for you FY
    
if manualChecks
    % let's check the data again
    D = spm_eeg_load(inputFile); % load our first dataset
    oslview(D);
end

%% artefact detection pt 1
% we probably still have artefacts in our data, so we'll want to remove those
% we definitely want to do this before we run an ICA, and exclude bad
% timepoints or the ICA might do weird stuff: https://ohba-analysis.github.io/osl-docs/matlab/osl_example_africa.html
if doArtefacts
    % add the artefactPrefix
    outputFile = addPrefix(inputFile,artefactPrefix);
    if ~exist(outputFile,'file') || overwrite
        S=[];
        S.D=inputFile;
        S.outfile=outputFile;
        D=spm_eeg_copy(S); % copy D object and add prefix
        
        modalities = {'MEGMAG','MEGPLANAR','EEG'};
        
        % look for bad channels
        D = osl_detect_artefacts(D,'badchannels',true,'badtimes',false,'modalities',modalities);
        
        % then look for bad time segments
        D = osl_detect_artefacts(D,'badchannels',false,'badtimes',true,'modalities',modalities);
        
        D.save;
        
    end % exist/overwrite if
end % doArtefacts
% update inputFile
inputFile = outputFile; % now artefact prefix has been prepended
    
if manualChecks
    % this doesn't remove bad timepoints (segments) or channels - just
    % marks them for future manipulation
    % we can view these and we can also use oslview to manually mark and
    % remove artefacts
    D = spm_eeg_load(inputFiles); % load our dataset
    D = oslview(D);
    % if we want to do manual removal, remember to save!
    D.save();
end


%% ICA!
% so now we'll use an ICA program (osl_africa here) to decompose the sensor data into a set of independent components (time courses), and remove artefactual ones (like eyeblinks and heartbeats etc).

if doICA
    % add the icaPrefix
    outputFile = addPrefix(inputFile,icaPrefix);
    if ~exist(outputFile,'file') || overwrite

    rng(randomSeed); % random seed for reproducing ICA decomposition

    % we need to prepare an EEG layout for africa
    D = spm_eeg_load(inputFile);
    cfg = [];
    cfg.output = 'eeg_layout.lay';
    cfg.elec = D.sensors('EEG');
    eeg_layout = ft_prepare_layout(cfg);

    % osl_africa usually doesn't autosave after running, but it
    % does for ICA, because it's time consuming
    % so we load our input files
    D = spm_eeg_load(inputFile);
    
    % then we'll make a new copy with a different name so it doesn't save over the previous step
    D = D.copy(outputFile); 

    % now each time we call osl_africa, call it as:
    % D = osl_africa...
    % because of the autosave, the D in memory won't be the
    % same D saved to disk!

    if doICA < 2 % if less than 2 we'll (re)run the ica, else we'll skip (with the intention of allowing us to just run the manual checks on the existing ica
    
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
        D.save(); % do we need this after auto ica? not sure

    end % check to see if running ICA

    if manualChecks

            % WIP: this does a manual identification, but I haven't
            % checked it works/saves properly and whatnot
            
            D = spm_eeg_load(outputFiles); % load this, in case we're just re-running for doing manual checks

            % we then want to redo the artefact rejection to make changes to the assignment if we like
            osl_africa(D,...
                'do_ica',false,...
                'do_ident','manual',...
                'used_maxfilter',true,...
                'artefact_channels',{'EOG','ECG'}...
                );
            % now will look for EEG
            osl_africa(D,...
                'do_ica',false,...
                'do_ident','manual',...
                'eeg_layout',eeg_layout.cfg.output,...
                'modality',{'EEG'},...
                'used_maxfilter',true,...
                'artefact_channels',{'EOG','ECG'}...
                );

            D.save(); % we do need to save the MEEG object to commit marked bad components to disk after the manual ica

        end % end manual checks

    end % exist/overwrite if
end % doICA
% update inputFile
inputFile = outputFile; % now ica prefix has been prepended

if manualChecks
    % osl africa saves an 'online montage' (a linear combination of the original sensor data) that we can explore
    D = spm_eeg_load(inputFile); % load our first dataset
    has_montage(D) % should show available montages
    D_pre_ica=D.montage('switch',0); % grab the original data
    D_post_ica=D.montage('switch',1); % grab the data transformed by the ica
    % we can plot this
    EOGsensor = 308;
    sensorToCheck = 306;
    figure;
    subplot(2,1,1);
    plot(D_pre_ica.time(1:10000),D_pre_ica(EOGsensor,1:10000)); % takes first 10000 sample points
    title('ECG channel')
    xlim([10 20]);
    xlabel('Time (s)');

    subplot(2,1,2);
    plot(D_pre_ica.time(1:10000),D_pre_ica(sensorToCheck,1:10000)); % takes first 10000 sample points
    title('ECG contaminated channel')
    xlim([10 20]);
    hold on;
    plot(D_post_ica.time(1:10000),D_post_ica(sensorToCheck,1:10000),'r');
    xlim([10 20]);
    xlabel('Time (s)');
    legend({'pre ICA' 'post ICA'});
    
    clear EOGsensor sensorToCheck
end

%% epoch the data
% so now we want to epoch the data
% we can do this with SPM (spm.meeg.preproc.epoch) or OSL (osl_epoch)
% both require the same inputs: a trl file with the epoch information and a list of the condition labels
% this pipeline has an megtriggers script that very manually walks through the meg data and codes the trigger events, since this has been very troublesome in the past. we'll use the output of that and feed it to osl, rather than using osl's native trigger recognition etc
% this lets us have more visibility on the triggers and deal with the inevitable partipant-specific weirdness, while also taking advantage of the osl recognition of the bad channel/segment markings which to this point have not been removed

if doEpoch
    % add the icaPrefix
    outputFile = addPrefix(inputFile,epochPrefix);
    if ~exist(outputFile,'file') || overwrite

        % first we'll load up the continuous data in case we want to check it later
        D_continuous=spm_eeg_load(inputFile);
        D_continuous=D_continuous.montage('switch',0); % make sure we select the data, and not the ica montage --- epoching will carry over to the ica montage
        
        if manualChecks
            % osl (I think?) lets us visualise what'll happen
            % trials go from 'o' to 'x'
            % bad trials are in black
            report.trial_timings(inputFile, epochInfo);
        end

        S=[];
        S = epochInfo; % we created this when we were converting and merging
        S.D = inputFile;
        D = osl_epoch(S); % note that osl_epoch prefixes the file it creates with an e!

        if manualChecks
            % should mirror our earlier check
            report.trial_timings(D);
            disp('Bad trial indices:');
            disp(D.badtrials);
        end

    end % exist/overwrite if
end % doEpoch
% update inputFile
inputFile = outputFile; % now epoch prefix has been prepended

%% artefact detection pt 2
% now we'll re-run our artefact detection on our epoched data
if doArtefacts
    % add the artefactPrefix
    outputFile = addPrefix(inputFile,artefactPrefix);
    if ~exist(outputFile,'file') || overwrite
        S=[];
        S.D=inputFile;
        S.outfile=outputFile;
        D=spm_eeg_copy(S); % copy D object and add prefix
        
        modalities = {'MEGMAG','MEGPLANAR','EEG'};
        
        % look for bad channels
        D = osl_detect_artefacts(D,'badchannels',true,'badtimes',false,'modalities',modalities);
        
        % then look for bad time segments
        D = osl_detect_artefacts(D,'badchannels',false,'badtimes',true,'modalities',modalities);
        
        D.save;

        if manualChecks
            % osl also has a manual tool for this
            % just click 'quit' on channels you want to skip (like EOG)
            D=osl_rejectvisual(D,[-0.2 0.4]);
            D.save();
        end
        
    end % exist/overwrite if
end % doArtefacts
% update inputFile
inputFile = outputFile; % now artefact prefix has been prepended

%% final checks
% so now we'll have a dataset that has two montages (raw and ica), epoched into trials, with good and bad trials marked
% we can get the indices of good trials: D.indtrial('condition of interest','good')
% and we can use that to get good trials in the ica montage of the data: D=D.montage('switch',1);
% we'll do all of that checking in a different script though

osl_shutdown


return
end

function outFiles = addPrefix(files,prefix)

for fileIdx = 1:length(files)
    [pathstr,name,ext] = fileparts(files{fileIdx});
    outFiles{fileIdx} = sprintf(['%s/' prefix '%s%s'],pathstr,name,ext);
end

return
end

function a4_3_preProcessing(thisSubject,datadir,toolsdir,runLocal)
% this is the best version of my osl pipeline
% it will produce an spm meeg object with the raw data and 2 montages
% 1 is the meg data transformed by ica
% 2 is the eeg data transformed by ica
% on the way, it filters, does a bit of automatic artefact removal besides
% but I didn't trust the automatic ica from osl, and found their manual
% very difficult to understand and trying to plot the results was also
% difficult and didn't look right so I swapped to fieldtrip which has more
% support both online and at the CBU
% I do recommend reading through the osl tutorials for some good ideas about
% what to do and why, but probably implement those things with fieldtrip or some stapling
% together of your favourite tools

%% this was an attempt to do all preprocessing in one file, toggling on and off bits
% that was hard, so just use it to run preprocessing and don't try to
% toggle

addpath /hpc-software/matlab/cbu/
addpath('/neuro/meg_pd_1.2/');       % FIFACCESS toolbox if use some meg_misc functions below
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/');
addpath(fullfile(toolsdir,'fieldtrip-20190410')); ft_defaults; % not needed if running locally, but needed if running on the cluster for some reason
% addpath(genpath(fullfile(toolsdir,'eeglab13_5_4b')))
addpath(genpath(fullfile(toolsdir,'meg_data')))
oslDir = fullfile('/imaging/woolgar/projects/Dorian/Toolboxes/osl'); % osl wants this if run in shared user mode
addpath(fullfile(oslDir,'osl-core')); % osl_startup can error without this
if runLocal; oslUserMode = 'user'; elseif ~runLocal; oslUserMode = 'shared';end

rootDir = fullfile(datadir,thisSubject.id);
inputFolder   = fullfile(rootDir,'MaxfilterOutput'); % where is the input? expects maxfiltered data, but you could do this on data that has not been maxfiltered
outputFolder  = fullfile(rootDir,'Preprocess'); % where should we put the output?
behavDir   = fullfile(datadir,'behavioural'); % where is our behavioural data?
behavfile = fullfile(behavDir,[thisSubject.id '_Evaccum.mat']); % what is it called?
megrtfile = fullfile(behavDir,[thisSubject.id '_MEGRTs.mat']); % and since I have put my MEG triggers somewhere else, we will specify that too
trlfile  = fullfile(outputFolder, 'run%s_trl.mat'); % this contains our trial related info for epoching - should exist (generated with megTriggers script)
addpath(rootDir);


manualChecks = 0; % enable manual checking and selecting throughout
overwrite = 0;

% if any of these are 1, then it will just do it, if 0 will assume it has been done and skip processing but still search for the appropriately prefixed file
doConvert = 1; % convert to SPM M/EEG object?
doFilter = 1; % filter line noise and do a high pass filter?
doArtefacts = 1; % detect bad channels and trials?
doICA = 1; % run our ICA? if 2 (or higher) will also do a manual ica rejection (if manualChecks is on)
doEpoch = 1; % epoch the data?

randomSeed = 7; % a random seed for ICA decomposition with oslafrica (use this same seed for reproducability)
filterPrefix = 'f'; % most programs default to an f, but note we filter twice so will have a file with ff prepended
artefactPrefix = 'b';
icaPrefix = 'ica'; % we add this to the file ourselves
epochPrefix = 'e'; % we also choose this

%% print subject

disp('this is subject:')
disp(thisSubject.id);

%% convert and merge

inputFiles = {}; % init this

% first we'll convert our files from fif to spm meeg object and move our files to the preprocessing folder
% maxfilter in our pipeline has first run sss then transformed to a
% the coordinate frames of the middle run for easy averaging, so we will do
% this not on the sss file, but the trans file
% NOTE: I've made an error here and used the trans file, but this is
% actually maxfilter adjusted to the default coord frames NOT the
% coord frames of the middle run as I suggest above---read about maxfilter
% to understand the difference
for runNum = 1:numel(thisSubject.meg_runs)
    fprintf('converting %.0f of %.0f raw files\n',runNum,numel(thisSubject.meg_runs))
    % input files from input folder
    tmpIn = sprintf([inputFolder,filesep,thisSubject.meg_labs{runNum},'_trans.fif']);
    % output files will go into output folder
    tmpOut = sprintf([outputFolder,filesep,thisSubject.meg_labs{runNum},'_trans.mat']);
    % if we haven't already, lets convert fif files to SPM M/EEG objects
    if ~exist(tmpOut,'file') || (overwrite && doConvert)
        S=[];
        S.dataset       = tmpIn;
        S.outfile       = tmpOut;
        S.save          = 0;
        S.reviewtrials  = 0;
        S.channels      = 'all';
        S.continuous    = 1;
        S.checkboundary = 0;
        
        spm_eeg_convert(S);
    else
        disp('file exists and not overwriting')
    end % doConvert if statement
    
    % collate these
    inputFiles{runNum} = tmpOut;
    
    % also, we need to collate the epoching trl information, which is currently in different files according to run
    thisTrlFile = sprintf(trlfile,num2str(runNum));
    thisEpochInfo = load(thisTrlFile); % should have vars 'trl' and 'conditionlabels' - this is how both osl and spm want it
    % to use our little trlfile merge function, we need an index---which
    % trl rows belong to what file
    thisEpochInfo.runcode = repmat(runNum,numel(thisEpochInfo.conditionlabels),1);
    % now combine them
    if ~exist('epochInfo','var')
        epochInfo = thisEpochInfo;
    else
        epochInfo.trl = [epochInfo.trl;thisEpochInfo.trl];
        epochInfo.conditionlabels = [epochInfo.conditionlabels;thisEpochInfo.conditionlabels];
        epochInfo.runcode = [epochInfo.runcode;thisEpochInfo.runcode];
    end
end; clear thisEpochInfo thisTrlFile % run loop

% now let's merge those files into one
inputFile = sprintf([outputFolder,filesep,'allRuns_trans.mat']);
newTrlFile = sprintf([outputFolder,filesep,'allRuns_trl.mat']);
disp('merging files')
if ~exist(inputFile,'file') || ~exist(newTrlFile,'file') || (overwrite && doConvert)
    % this function will automatically save D as the second vargin
    [D,epochInfo] = continuousMeegMerger(inputFiles,inputFile,epochInfo);
    % let's save trl info/epoch info in a format that spm would like it
    trl = epochInfo.trl;
    conditionlabels = epochInfo.conditionlabels;
    save(newTrlFile,'trl','conditionlabels'); clear trl conditionlabels
else
    disp('file exists and not overwriting')
end
clear tmpIn tmpOut

% so we now have an SPM M/EEG object saved
if manualChecks
    % put a breakpoint here and you can have a play
    % this will be boring if you haven't epoched yet (if data is
    % continuous), but principles apply throughout the pipeline
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
end

% we are going to swap now to use OSL for some of the next bits of the pipeline
% though note that we'll use spm_eeg functions occasionally and can swap bits and
% pieces out for other toolboxes (e.g. fieldtrop/eeglab)

osl_startup(oslDir,oslUserMode)

if manualChecks
    % we can now use OSLVIEW to look manually at continuous data (not epoched data)
    % we can also interactively flag bad channels and time periods with it
    D_raw = spm_eeg_load(inputFile);
    if strcmp(D_raw.type,'continuous')
        %         oslview(D_raw); % can also: D = oslview(D); if you want to edit the data now
    else
        disp('not continuous data, cant use olsview')
    end
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
        S.freq = 0.5; % removes the slow drift
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
    else
        disp('file exists and not overwriting')
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
% note that there's an interaction between downsampling and filtering though
% so we might want to downsample earlier---maxfilter can also downsample
% for you FYI
% but generally, we probably want to filter first to avoid any aliasing issues (https://en.wikipedia.org/wiki/Nyquist_frequency#Aliasing).

if manualChecks
    % let's check the data again
    D_filtered = spm_eeg_load(inputFile);
    if strcmp(D_filtered.type,'continuous')
        %oslview(D_raw); % if you'd like to compare these
%         oslview(D_filtered); % can also: D = oslview(D); if you want to edit the data now
    else
        disp('not continuous data, cant use olsview')
    end
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
    else
        disp('file exists and not overwriting')
    end % exist/overwrite if
end % doArtefacts
% update inputFile
inputFile = outputFile; % now artefact prefix has been prepended

if manualChecks
    % this doesn't remove bad timepoints (segments) or channels - just
    % marks them for future manipulation
    % we can view these and we can also use oslview to manually mark and
    % remove artefacts
    D_deartefacted = spm_eeg_load(inputFile); % load our dataset
    if strcmp(D_deartefacted.type,'continuous')
        %oslview(D_raw); % if you'd like to compare these
        %oslview(D_filtered); % if you'd like to compare these
%         D_deartefacted = oslview(D_deartefacted);
    else
        disp('not continuous data, cant use olsview')
    end
    % if we want to do manual removal, remember to save!
    D_deartefacted.save();
end


%% ICA!
% so now we'll use an ICA program (osl_africa here) to decompose the sensor data into a set of independent components (time courses), and remove artefactual ones (like eyeblinks and heartbeats etc).

if doICA
    % add the icaPrefix
    outputFile = addPrefix(inputFile,icaPrefix);
    disp('doing ica')
    if ~exist(outputFile,'file') || overwrite
        
        rng(randomSeed); % random seed for reproducing ICA decomposition
        
        % we need to prepare an EEG layout for africa
        % I think spm_eeg_convert does this, but I can't figure out how to
        % access it, so this is easier
        D = spm_eeg_load(inputFile);
        cfg = [];
        cfg.output = [outputFolder filesep 'eeg_layout.lay'];
        cfg.elec = D.sensors('EEG');
        eeg_layout = ft_prepare_layout(cfg);
        
        % osl_africa usually doesn't autosave after running, but it
        % does for ICA, because it's time consuming
        % so we load our input files
        D = spm_eeg_load(inputFile);
        
        % then we'll make a new copy with a different name so it doesn't save over the previous step
        disp('copying file with new filename')
        D = D.copy(outputFile);
        
        % now each time we call osl_africa, call it as:
        % D = osl_africa...
        % because of the autosave, the D in memory won't be the
        % same D saved to disk!
        
        
        % with no modality specified will look for MEG
        D = osl_africa(D,...
            'do_ica',true,...
            'used_maxfilter',true,...
            'artefact_channels',{'EOG','ECG'},...
            'montagename','ols africa ica denoised meg'...
            );
        % now will look for EEG
        D = osl_africa(D,...
            'do_ica',true,...
            'auto_artefact_chans_corr_thresh',.35,...
            'eeg_layout',eeg_layout.cfg.output,...
            'modality',{'EEG'},...
            'used_maxfilter',true,...
            'artefact_channels',{'EOG','ECG'},...
            'montagename','ols africa ica denoised eeg'...
            );
        % automatically does artefact identification (default), and removes them
        
        % this creates D.ica
        D.save(); % do we need this after auto ica? not sure
        
    else
        disp('file exists and not overwriting')
    end % exist/overwrite if
end % doICA
% update inputFile
inputFile = outputFile; % now ica prefix has been prepended

if manualChecks
    
    if doICA > 1
        disp('now doing manual ica')
        % WIP: this does a manual identification, but I haven't
        % checked it works/saves properly and whatnot
        
        D = spm_eeg_load(inputFile); % load this, in case we're just re-running for doing manual checks
        cfg = [];
        cfg.output = [outputFolder filesep 'eeg_layout.lay'];
        cfg.elec = D.sensors('EEG');
        eeg_layout = ft_prepare_layout(cfg);
        
        % we then want to redo the artefact rejection to make changes to the assignment if we like
        D = osl_africa(D,...
            'do_ica',false,...
            'do_ident','manual',...
            'used_maxfilter',true,...
            'artefact_channels',{'EOG','ECG'},...
            'montagename','ols africa ica denoised meg'...
            );
        % now will look for EEG
        D = osl_africa(D,...
            'do_ica',false,...
            'do_ident','manual',...
            'eeg_layout',eeg_layout.cfg.output,...
            'modality',{'EEG'},...
            'used_maxfilter',true,...
            'artefact_channels',{'EOG','ECG'},...
            'montagename','ols africa ica denoised eeg'...
            );
        
        D.save(); % we do need to save the MEEG object to commit marked bad components to disk after the manual ica
    end
    
    % osl africa saves an 'online montage' (a linear combination of the original sensor data) that we can explore
    D = spm_eeg_load(inputFile); % load our first dataset
    has_montage(D) % should show available montages
    
    D_pre_ica = D.montage('switch',0); % grab the original data
    D_post_ica = D.montage('switch',2); % grab the data transformed by the ica
    
    %oslview(D_raw); % if you'd like to compare these
    %oslview(D_filtered); % if you'd like to compare these
    %oslview(D_deartefacted); % if you'd like to compare these
    
    % but now we can start to plot interesting stuff
    artefactSensor = 2; % change this to an ECG or EOG sensor - something that is artefactual
    sensorToCheck = 4; % change this to a sensor near your artefact sensor - we will compare it to the artefact before and after ICA
    % note that the way we've done this means there are two ica montages,
    % one for eeg and one for meg, so if we have selected the wrong
    % montage+sensor combination, the plots will only show the artefact
    % sensor
    figure;
    subplot(2,1,1);
    plot(D_pre_ica.time(1:10000),D_pre_ica(artefactSensor,1:10000)); % takes first 10000 sample points
    title('artefact channel')
    xlim([0 10]);
    xlabel('Time (s)');
    
    subplot(2,1,2);
    plot(D_pre_ica.time(1:10000),D_pre_ica(sensorToCheck,1:10000)); % takes first 10000 sample points
    title('artefact contaminated channel')
    xlim([0 10]);
    hold on;
    plot(D_post_ica.time(1:10000),D_post_ica(sensorToCheck,1:10000),'r');
    xlim([0 10]);
    xlabel('Time (s)');
    legend({'pre ICA' 'post ICA'});
    
    clear artefactSensor sensorToCheck
    close all
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
    disp('epoching')
    if ~exist(outputFile,'file') || overwrite
        
        % first we'll load up the continuous data in case we want to check it later
        D_continuous = spm_eeg_load(inputFile);
        D_continuous = D_continuous.montage('switch',0); % make sure we select the data, and not the ica montage --- epoching will carry over to the ica montage
        
        S=[];
        S = epochInfo; % we created this when we were converting and merging
        S.D = inputFile;
        D = osl_epoch(S); % note that osl_epoch prefixes the file it creates with an e!
        
        if manualChecks
            % osl (I think?) lets us visualise what'll happen
            % trials go from 'o' to 'x'
            % bad trials are in black
            % this takes bloody ages mate
%             report.trial_timings(D);
%             report.bad_trials(D);
%             report.bad_channels(D,{'MEGMAG','MEGPLANAR','EEG'});
            disp('index of bad trials:');
            disp(D.badtrials);
        end
        
    else
        disp('file exists and not overwriting')
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
        
    else
        disp('file exists and not overwriting')
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
disp('done')

return
end

function outFiles = addPrefix(files,prefix)

if iscell(files)
    for fileIdx = 1:length(files)
        [pathstr,name,ext] = fileparts(files{fileIdx});
        outFiles{fileIdx} = sprintf(['%s/' prefix '%s%s'],pathstr,name,ext);
    end
else
    [pathstr,name,ext] = fileparts(files);
    outFiles = sprintf(['%s/' prefix '%s%s'],pathstr,name,ext);
end



return
end

function [Dnew, trlsMerged] = continuousMeegMerger(spmfilelist,newfile,trls)
% based on osl_meeg_merge() by Adam Baker 2015
% Merge continuous MEEG objects in filelist to a new single object newfile
% also optionally combine your trl structures for those objects
% trls must have 3 fields:
% trls.trl and trls.conditionlabels as per expected input to
% spm_eeg_epoch trlfile
% plus trls.runcode which indexes the runs as per spmfile list
% so 3 spmfiles, runcode would index which of trls.trl/conditionlabels
% belong to which file


D = cell(size(spmfilelist));
for i = 1:length(spmfilelist)
    D{i} = spm_eeg_load(spmfilelist{i});
end

num_samples = cellfun(@nsamples,D);

end_samples = cumsum(num_samples);
start_samples = [1 end_samples(1:end-1)+1];

if exist('trls','var')
    trlsMerged.trl = zeros(size(trls.trl));
    trlsMerged.conditionlabels = cell(size(trls.conditionlabels));
    for i = 1:length(spmfilelist)
        trlIdx = find(trls.runcode == i);
        
        theseTrlConds = trls.conditionlabels(trlIdx);
        trlsMerged.conditionlabels(trlIdx) = theseTrlConds;
        
        theseTrlTimes = trls.trl(trlIdx,:);
        theseTrlTimes(:,1:2) = theseTrlTimes(:,1:2)+(start_samples(i)-1);
        trlsMerged.trl(trlIdx,:) = theseTrlTimes;
    end
else
    trlsMerged = NaN;
end

% Copy data:
Dnew = clone(D{1},newfile,[D{1}.nchannels,sum(num_samples),1]);
for i = 1:length(spmfilelist)
    Dnew(:,start_samples(i):end_samples(i),:) = D{i}(:,:,:);
end

% Ensure bad channels are consistent:
Dnew = badchannels(Dnew,1:Dnew.nchannels,0);
Dnew = badchannels(Dnew,unique(cell2mat(cellfun(@badchannels,D,'uniformoutput',0))),1);

% Correct events:
ev = cellfun(@events,D,'uniformoutput',0);
for i = 1:length(spmfilelist)
    for j = 1:numel(ev{i})
        ev{i}(j).time = ev{i}(j).time + (start_samples(i) - 1)/Dnew.fsample;
    end
end
ev = cat(1,ev{:});
Dnew = events(Dnew,1,ev);

Dnew.save;

return
end

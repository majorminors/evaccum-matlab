function a4_preProc_atypical_artefacts(thisSubject,datadir,toolsdir,overwrite,fileMods,ica,manual)
% here we will:
% 1. view the data if we want (just uncomment and run locally)
% 2. do some basic artefact detection either manually or automatically
% this can either produce a de-artefacted file for subsequent ica (for
% which we should remove crazy artefacts so it doesn't use up components
% and also because ica can reintroduce weird stuff), or it can be run on
% ica-ed data to do aretefact checking on our trials after the fact to
% reject any remaining weird trials before analysis
% I coded this, but then switched to Olaf's ica scripts for mne python, in
% which he rejects artefacts pre-ica, so we don't need to do this pre-ica
% anymore


if ~exist('overwrite','var'); overwrite = 0; end
if ~exist('ica','var'); ica = 0; end
% do manual if unspecified (or don't do it if unspecified but doing ica)
if ~exist('manual','var'); if ica == 1; manual = 0; else; manual = 1; end; end
% fileMods = cell array with letters in file in, extension in file in,
% letters to add in file out
if ~exist('fileMods','var'); fileMods = {'' '.mat' 'C'}; end

% set up paths
addpath /hpc-software/matlab/cbu/
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/');
ftDir = fullfile(toolsdir,'..','..','..','Toolboxes','fieldtrip'); % I pass this to my saving function too
addpath(ftDir); ft_defaults;
addpath(genpath(fullfile(toolsdir,'lib')))
rootDir = fullfile(datadir,thisSubject.id);
preProcDir = fullfile(rootDir,'Preprocess');
addpath(rootDir);

% file definition
% sprintf to enter missing data ('%s')
if isempty(fileMods{1}); inputFileName = fullfile(preProcDir,'run%s_transmid%s'); else
    inputFileName = fullfile(preProcDir,['run%s_' fileMods{1} '_transmid%s']); end
outputFileName = fullfile(preProcDir,['run%s_' fileMods{3} fileMods{1} '_transmid%s']);
trlfile = fullfile(preProcDir, 'run%s_trl.mat'); % this contains our trial related info for epoching - should exist (generated with megTriggers script)
behavfile = fullfile(datadir,'behavioural',[thisSubject.id '_MEGRTs.mat']); % grab our behavioural data + output from a3_megtriggers

%% print subject for logging

disp('this is subject:')
disp(thisSubject.id);

% first we'll load and epoch it all, so we can do other stuff while we wait
goodRuns = [];
for runNum = 1:numel(thisSubject.meg_runs)
    
    thisInFile = sprintf(inputFileName,num2str(runNum),fileMods{2});
    if ica; ext = '.fif'; else; ext = '.mat'; end
    thisOutFile{runNum} = sprintf(outputFileName,num2str(runNum),ext); clear ext
    
    if exist(thisInFile,'file')
        if ~exist(thisOutFile{runNum},'file') || overwrite
            
            fprintf('loading and epoching %.0f of %.0f run files\n',runNum,numel(thisSubject.meg_runs))
            
            %% set up a fieldtrip config
            cfg = [];
            cfg.continuous = 'yes'; % this data is not epoched, it's continuous
            if strcmp(fileMods{2},'.fif')
                cfg.dataset = thisInFile;
                rawData{runNum} = ft_preprocessing(cfg); % load it
            elseif strcmp(fileMods{2},'.mat')
                rawData{runNum} = ft_preprocessing(cfg,loadFtData(thisInFile)); % load it
            end
            
            %% we can look at the data
            % e.g. to see if we've missed anything in filtering
            % or do manual artefact detection at the outset if you know something
            % weird happened
            %         cfg.channel = {'MEG'}; cfg.ylim = [-5e-11 5e-11];% grab just these, so we aren't plagued by all the STI channels and triggers etc
            %         cfg.channel = {'EEG'}; cfg.ylim = [-0.0005 0.0005];
            %         cfg.blocksize = 10; % how many seconds of data to show at once
            % choose one of
            %         cfg.viewmode = 'butterfly';
            %         cfg.viewmode = 'vertical';
            %         cfg = ft_databrowser(cfg, rawData{runNum});
            
            % we can also look at the frequencies
            %         cfg = [];
            %         cfg.output = 'pow';
            %         cfg.channel = 'EEG';
            %         cfg.method = 'mtmfft';
            %         cfg.taper = 'hanning';
            %         cfg.pad = 'nextpow2';
            %         cfg.foilim = [0.3 rawData{runNum}.fsample/2];
            %         freqChk = ft_freqanalysis(cfg,rawData{runNum});
            %         figure; plot(freqChk.freq,mag2db(freqChk.powspctrm));
            
            %% NOTE: if you looked at your data, be sure to re-run the previous ft config section before moving on
            
            
            %% now we prep the files for some preliminary artefact removal
            % I do this for all runs before artefact checking, so we just have one chunk of loading time
            % (which we can use to go get coffee or something)
            
            if ica
                % we want to remove atypical artefacts, prior to any ICA
                % this is (a) because we don't want to waste components on obvious nonsense and
                % (b) because ICA can re-introduce crazy artefacts to our data (see e.g.
                % https://ohba-analysis.github.io/osl-docs/matlab/osl_example_africa.html#12)
                % for ica, we actually prefer continuous data, so we'll cut the data into arbitrary
                % chunks of 1 second, so we can compare time segments, then stitch it back together again,
                cfg = [];
                cfg.length = 1;
                cfg.overlap = 0;
                epochedData{runNum} = ft_redefinetrial(cfg,rawData{runNum});
                % a reference for how to put back together again, though I do this when we save:
                %         cfg = [];
                %         cfg.continuous = 'yes';
                %         stitchedData{runNum} = ft_redefinetrial(cfg,epochedData{runNum});
                
            else
                % for decoding, we probably don't need to worry about ICA,
                % and for post ica processing, the same, so here we
                % can just cut this immediately into our trials of interest and
                % work with that
                
                % I have already produced a trial definition in spm/osl trl file
                % format, so I will not do it as per fieldtrip tutorials, but instead
                % import that file in the fieldtrip format (which is very similar)
                cfg = [];
                % once, fieldtrip had a cfg.trialfun option that let you
                % specify a custom trial function to define trials.
                % it doesn't have this anymore, but we can use it in the
                % same way. it just needs to return something we can put in
                % cfg.trl that does [trialStart trialEnd trialOffset; etc]
                % for each trial
                % I just arbitrarily add these to cfg struct so the function has something to read to find the files it needs
                cfg.trlfile = sprintf(trlfile,num2str(runNum));
                cfg.behavfile = behavfile;
                cfg.trl = trlFromFile(cfg);
                % and that should work---if you're using an older version
                % of fieldtrip, then you'd do basically the same thing by
                % doing:
                % cfg.trialfun = 'trlFromFile';
                epochedData{runNum} = ft_redefinetrial(cfg,rawData{runNum});
                
            end
            
            goodRuns = [goodRuns, runNum]; % collect this, so we can index into runs that we just processed later
        else
            disp(['the file:\n' thisOutFile{runNum} '\nexists, skipping'])
        end
    else
        warning(['the file:\n' thisInFile '\ndoes not exist, skipping'])
    end % end infile check
end

clear rawData

%% now loop through whatever runs we have prepped and do some rejection
count = 1;
while count <= numel(goodRuns)
    goodIdx = goodRuns(count);
    fprintf('artefact checking %.0f of %.0f run files\n',count,numel(goodRuns))
    
    %% prep layout
    % make an eeg layout
    cfg = [];
    cfg.elec = epochedData{goodIdx}.elec;
    eeglayout = ft_prepare_layout(cfg);
    % and an meg layout
    cfg = [];
    cfg.grad = epochedData{goodIdx}.grad;
    meglayout = ft_prepare_layout(cfg);
    % combine them
    cfg = [];
    combinedLayout = ft_appendlayout(cfg, eeglayout, meglayout);
    clear meglayout eeglayout
    
    if manual
        %% manual/visual rejection
        
        cfg = [];
        cfg.method = 'summary';
        cfg.keepchannel = 'yes'; % no is default
        % if we leave the default for keepchannel, then channels we remove from
        % the display will be removed from the data. this is annoying,
        % particularly because the e.g. eog channels will sometimes be quite
        % different from neighbours and you'll remove them by accident.
        % also if we remove channels then this would
        % change the channel numbers (e.g. removing EEG004 would make
        % EEG005 number 4 in later analyses) so keep that in mind
        % if this is pre-ica processing, or you're looking to remove some
        % eog/ecg artefacts I recommend eithers replacing with nans, zeroes, or better still
        % just use the output of maxfilter to eventually remove bad
        % channels, rather than trying to work out if you should remove them by hand
        if ica
            % ICA is better on continuous data, so lets fill trials with nans,
            % rather than removing them, though this will only work on
            % scripts that can handle nans (like fieldtrip's
            % implementation, but not mne without some editing of the
            % 'picks')
            cfg.keeptrial = 'nan';
        else
            % let's just ditch bad trials
            cfg.keeptrial = 'no'; % no is default
        end
        cfg.metric = 'zvalue';
        cfg.layout = combinedLayout;
        
        
        quitCheckLoop = 0;
        while ~quitCheckLoop
            
            
            % let's subselect our channels of interest
            cfg.channel = {'MEG' 'EEG'};
            % but note that if we were removing channels, it would assume we
            % wanted to remove any channels not sub-selected here (e.g. EOG
            % channels), so don't do this if you're removing bad channels by
            % hand or replacing them by nans etc or have finished with those channels)
            cleanedData = ft_rejectvisual(cfg, epochedData{goodIdx});
            
            % lets check we didn't make a mistake
            inputResponse = input('happy with that ([y]/n)? ','s'); % input response as a string
            if isempty(inputResponse); inputResponse = 'y'; end % a default response
            
            if strcmp(inputResponse,'y')
                disp('happy, saving and moving on')
                % lets loop through a prompt asking if we actually want to save that or
                % re-do it
                % you can send it to the cluster if the data isn't too
                % big---saves a bit of time                 on saving
                % disp('sending to cluster to save')
                % sendToCluster(@saveThis, {thisOutFile{goodIdx} cleanedData}, {fullfile(rootDir,'Preprocess') ftDir})
                % disp('sent')                
                disp('saving')
                saveThis(thisOutFile{goodIdx}, cleanedData)
                count = count+1;
                quitCheckLoop = 1;
            elseif strcmp(inputResponse,'n')
                disp('not happy, lets go again on that run')
                quitCheckLoop = 1;
                else
                disp('invalid response, try again')
            end
        end; clear quitLoop inputResponse
        count = count+1;
    else
        %% automatic rejection via z-value
        count = count+1;
        % so, fieldtrip has a nice tutorial for this:
        % https://www.fieldtriptoolbox.org/tutorial/automatic_artifact_rejection/
        % we really might want to do the squid jump removal, because we know the CBU has a
        % funny squid (at the time of writing) so an example would be:
        %         cfg = [];
        %         % if you want to view it to see if you want to change my defaults:
        %         cfg.artfctdef.zvalue.interactive = 'yes';
        %         cfg.artfctdef.zvalue.channel = 'MEG';
        %         cfg.artfctdef.zvalue.cutoff = 20;
        %         cfg.artfctdef.zvalue.trlpadding = 0;
        %         cfg.artfctdef.zvalue.artpadding = 0;
        %         cfg.artfctdef.zvalue.fltpadding = 0;
        %         cfg.artfctdef.zvalue.cumulative = 'yes';
        %         cfg.artfctdef.zvalue.medianfilter = 'yes';
        %         cfg.artfctdef.zvalue.medianfiltord = 9;
        %         cfg.artfctdef.zvalue.absdiff = 'yes';
        %         [cfg, artifact_jump] = ft_artifact_zvalue(cfg, epochedData{goodIdx});
        
        % but I'm much less enthusiastic about removing muscle artefacts
        % automatically, because this really seems like it should be done by hand (or,
        % honestly, not at all:
        % https://www.biorxiv.org/content/10.1101/2022.12.03.518987v1.full.pdf)
        
        % so we aren't going to do *any of that*, and instead just remove any
        % segments that seem crazy in relation to the others---anything
        % that e.g. the ICA might see and use as a component, or worse,
        % introduce into the data on application (e.g.
        % https://ohba-analysis.github.io/osl-docs/matlab/osl_example_africa.html#12)
        % we'll essentially just use the same kind of process as we would
        % for visual rejection, made possible by fieldtrips recent
        % implementation of badchannel and badsegment detection
        
        cfg = [];
        % here's an example of badchannels (NO CHANNELS SPECIFIED!)
        % but we won't worry about that for now
        %         cfg.layout = combinedLayout;
        %         cfg.method = 'distance';
        %         cfg.neighbours = ft_prepare_neighbours(cfg, epochedData{goodIdx});
        %         cfg = ft_badchannel(cfg, epochedData{goodIdx});
        %         badMegChans = cfg.badchannel;
        
        if ica
            cfg.artfctdef.reject = 'nan'; % for ica we want continuous data, so just make them nans
        else
            cfg.artfctdef.reject = 'complete'; % otherwise remove 'complete' trials
        end
        
        cfg.metric = 'zvalue';
        cfg.threshold = 4;
        cfg.channel = 'MEG';
        cfg = ft_badsegment(cfg, epochedData{goodIdx});
        cleanedData = ft_rejectartifact(cfg, epochedData{goodIdx});
        cfg.channel = 'EEG';
        cfg = ft_badsegment(cfg, cleanedData);
        disp('dont worry about warnings about nans---this wont spread, its contained to the trial youre naning')
        cleanedData = ft_rejectartifact(cfg, cleanedData);
        disp('saving')
        saveThis(thisOutFile{goodIdx}, cleanedData)
        count = count+1;
        
    end % manual toggle
    
end


end

function saveThis(savePath, data)
% requires the save directory and ft on the path

disp('saving data:')
disp(data)
disp('to:')
disp(savePath)
if endsWith(savePath,'fif')
    ft_defaults;
    cfg = [];
    cfg.continuous = 'yes';
    stitchedData = ft_redefinetrial(cfg,data);
    fieldtrip2fiff(savePath,stitchedData);
else
    save(savePath,'data','-v7.3')
end
disp('done')

end
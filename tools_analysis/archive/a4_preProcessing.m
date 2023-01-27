function a4_4_preProcessing(thisSubject,datadir,toolsdir)
% here we apply various filters to the data

% set up paths
addpath /hpc-software/matlab/cbu/
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/');
addpath(fullfile(toolsdir,'fieldtrip-20190410')); ft_defaults;
addpath(genpath(fullfile(toolsdir,'lib')))
rootDir = fullfile(datadir,thisSubject.id);
maxFilterFolder   = fullfile(rootDir,'MaxfilterOutput'); % where is the input? expects maxfiltered data, but you could do this on data that has not been maxfiltered
preProcFolder  = fullfile(rootDir,'Preprocess'); % where should we put the output?
behavDir   = fullfile(datadir,'behavioural'); % where is our behavioural data?
behavfile = fullfile(behavDir,[thisSubject.id '_Evaccum.mat']); % what is it called?
megrtfile = fullfile(behavDir,[thisSubject.id '_MEGRTs.mat']); % and since I have put my MEG triggers somewhere else, we will specify that too
trlfile  = fullfile(preProcFolder, 'run%s_trl.mat'); % this contains our trial related info for epoching - should exist (generated with megTriggers script)
addpath(rootDir);

%% print subject for logging

disp('this is subject:')
disp(thisSubject.id);

for runNum = 1:numel(thisSubject.meg_runs)
    fprintf('preprocessing %.0f of %.0f files\n',runNum,numel(thisSubject.meg_runs))
    
    %% set up a fieldtrip config
    cfg = [];
    cfg.continuous = 'yes'; % this data is not epoched, it's continuous
    % get the path to the file
    % here I'm using the transmid output of my maxfiltering---this is the
    % run transformed to the coordinates of the middle run, so all my runs
    % are comparable
    cfg.dataset = sprintf([maxFilterFolder,filesep,thisSubject.meg_labs{runNum},'_transmid.fif']);
    thisHdr = ft_read_header(cfg.dataset); % grab some of the details (e.g. channels etc) for later
    
    %% we can look at the data
    % e.g. to see if we need to filter harmonics of line noise (see filtering section below)
    % or do manual artefact detection at the outset if you know something
    % weird happened
%     cfg.channel = {'MEG'}; cfg.ylim = [-5e-11 5e-11];% grab just these, so we aren't plagued by all the STI channels and triggers etc
% %     cfg.channel = {'EEG'}; cfg.ylim = [-0.0005 0.0005];
%     cfg.blocksize = 10; % how many seconds of data to show at once
%     % choose one of
%     cfg.viewmode = 'butterfly';
% %     cfg.viewmode = 'vertical';
%     cfg = ft_databrowser(cfg);

%% NOTE: if you looked at your data, be sure to re-run the previous ft config section before moving on

    %% now we want to define the trials
    % I have already produced a trial definition in spm/osl trl file
    % format, so I will not do it as per fieldtrip tutorials, but instead
    % import that file in the fieldtrip format (which is very similar)
    cfg.trialfun = 'trlFromFile'; % this is my own trial function that simply reads in from the spm/osl trl file I made in a3_megTriggers
    cfg.filename = sprintf(trlfile,num2str(runNum)); % I just arbitrarily add this to cfg struct so the function has something to read to find the file it needs
    cfg = ft_definetrial(cfg);
    % so I'll end up with epoch start, epoch end, offset, and condition codes

    %% now filter the data
    
    % we might want to high-pass filter to remove slow drift potentials/fields that
    % the faster signals of interest ride on. But this risks losing slow
    % fluctuations of brain potential, and introducing artefacts (ringing).
    % temporal maxfilter high-pass filters already, so see if you need to do more of that.
    % I have done my maxfilter at a 10 sec data buffer, which corresponds to  a cut-off freq of 0.1Hz (1/10s)
    % I'm pretty sure my old script filtered again at 0.5 but I think I'll
    % just leave this for now
    cfg.hpfilter = 'no';
    
    % we also might want to low pass here, the benefit being that it might
    % remove high-frequency variance that is often noise in relation to the
    % typically slower dynamics of interest. But this risks masking
    % high-frequency (e.g. gamma band) activity, or even just useful
    % information about the shape of some responses. I'm not going to
    % because I'm pretty sure I care about gamma.
    cfg.lpfilter = 'no';
    
    % what we almost certainly want to do is notch filter any line noise (+
    % harmonics of the line noise if necessary). Worth examining the data
    % at this point to see if you have line noise/harmonics. Obviously this
    % means you lose any information at this frequency.
    % fieldtrip has a specific kind of filter for this: DFT (a sharp
    % 'discrete fourier transform' but note this needs the data to have
    % lots of data either side (5-10 secs). should be ok if each run starts
    % with the 10-20 second recording that is standard at cbu), but will
    % throw an error if you don't have enough/will warn you that it's using
    % data mirroring to achieve this
    % alternatively you can use their notch filter, they call a band-stop
    % filter, so e.g. cfg.bsfilter & cfg.bsfreq = [49 51];
    cfg.dftfilter = 'yes';
    cfg.dftfreq = '50'; % we'll just do the line noise, but we could add harmonics here (.e.g [50 100])---the default will do [50 100 150]
    % it's worth noting that ica could probably pick up line noise and
    % harmonics too
    
    % note all filters default to a two-pass filter (back and forward)
    % which results in a zero-phase shift of erp components, but you can
    % change this in fieldtrip---both order and direction
    
end





end

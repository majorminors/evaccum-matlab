function a4_preProc02_filtering(thisSubject,datadir,toolsdir,manual,overwrite)
% here we will:
% 2. define our trials
% 3. filter the line noise and maybe other filters
% 2. save the file as a fif, so we can run Olafs ica on it

if ~exist('overwrite','var'); overwrite = 0; end

% set up paths
addpath /hpc-software/matlab/cbu/
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/');
addpath(fullfile(toolsdir,'..','..','..','Toolboxes','fieldtrip')); ft_defaults;
addpath(genpath(fullfile(toolsdir,'lib')))
rootDir = fullfile(datadir,thisSubject.id);
addpath(rootDir);

% file definition

inputFileName = fullfile(rootDir,'Preprocess','run%s_C_transmid.mat'); % I'll use thisSubject.meg_labs to cycle through with sprintf
outputFileName = fullfile(rootDir,'Preprocess','run%s_fC_transmid.fif'); % run number to full this out
trlfile  = fullfile(rootDir,'Preprocess', 'run%s_trl.mat'); % this contains our trial related info for epoching - should exist (generated with megTriggers script)
behavfile = fullfile(datadir,'behavioural',[thisSubject.id '_MEGRTs.mat']); % grab our behavioural data + output from a3_megtriggers

%% print subject for logging

disp('this is subject:')
disp(thisSubject.id);

for runNum = 1:numel(thisSubject.meg_runs)
    
    clear thisFile
    thisFile = sprintf(outputFileName,runNum);
    
    if ~exist(thisFile,'file') || overwrite
        
        fprintf('filtering %.0f of %.0f run files\n',runNum,numel(thisSubject.meg_runs))
        %% now filter the data
        
        % what we almost certainly want to do is notch filter any line noise (+
        % harmonics of the line noise if necessary). Worth examining the data
        % at this point to see if you have line noise/harmonics. Obviously this
        % means you lose any information at this frequency.
        % fieldtrip has a specific kind of filter for this: DFT (a sharp
        % 'discrete fourier transform' but note this is sharpest with
        % lots of padding either side (5-10 secs) and also can have weird effects if the line noise
        % isn't constant (https://www.fieldtriptoolbox.org/faq/why_is_there_a_residual_50hz_line-noise_component_after_applying_a_dft_filter/)
        % as such, you can instead use their notch filter, they call a band-stop
        % filter, so e.g. cfg.bsfilter & cfg.bsfreq = [49 51]; but this
        % comes with it's own considerations (e.g. you lose all that
        % information)
        cfg.dftfilter = 'yes';
        cfg.dftfreq = '50'; % we'll just do the line noise, but we could add harmonics here (.e.g [50 100])---the default will do [50 100 150]
        cfg.padding = '5';
        % it's worth noting that ica could probably pick up line noise and
        % harmonics too
        
        % note all filters default to a two-pass filter (back and forward)
        % which results in a zero-phase shift of erp components, but you can
        % change this in fieldtrip---both order and direction
        
        %% now we run the filtering
        filteredData = ft_preprocessing(cfg);
        
        % you might want other filters, but think about whether you should do
        % these pre or post ICA
        % you probably want a different padding for this than the dft filter, so you'd run it
        % seperately, clearing the cfg.dftfilter/freq and overridding the
        % padding, then running ft_preprocessing again.
        
        % we might want to high-pass filter to remove slow drift potentials/fields that
        % the faster signals of interest ride on. But this risks losing slow
        % fluctuations of brain potential, and introducing artefacts (ringing).
        % temporal maxfilter high-pass filters already, so you might not need
        % to do more of this.
        % I have done my maxfilter at a 10 sec data buffer, which corresponds to  a cut-off freq of 0.1Hz (1/10s)
        % I'm pretty sure my old script filtered again at 0.5 but I think I'll
        % just leave this for now
        %         cfg.hpfilter = 'no';
        
        % we also might want to low pass here, the benefit being that it might
        % remove high-frequency variance that is often noise in relation to the
        % typically slower dynamics of interest. But this risks masking
        % high-frequency (e.g. gamma band) activity, or even just useful
        % information about the shape of some responses. I'm not going to
        % because I'm pretty sure I care about gamma.
        %         cfg.lpfilter = 'no';
        
        %% save this
        
        disp('saving');
        save(thisFile,'filteredData','-v7.3');
        clear filteredData
        
    else
        disp([outputFileName ' exists, skipping'])
    end
    
end

end
function a4_preProc_filtering(thisSubject,datadir,toolsdir,overwrite,ica)
% here we will:
% 1. define our trials
% 2. filter the line noise and maybe apply other filters
% 3. save as either a fif (for mne ica) or fieldtrip .mat, based on the
% value of ica

if ~exist('overwrite','var'); overwrite = 0; end

% set up paths
addpath /hpc-software/matlab/cbu/
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/');
addpath(fullfile(toolsdir,'..','..','..','Toolboxes','fieldtrip')); ft_defaults;
addpath(genpath(fullfile(toolsdir,'lib')))
rootDir = fullfile(datadir,thisSubject.id);
preProcDir = fullfile(rootDir,'Preprocess');
addpath(rootDir);

% file definition

inputFileName = fullfile(rootDir,'MaxfilterOutput','Run%s_transmid.fif'); % I'll use thisSubject.meg_labs to cycle through with sprintf
outputFileName = fullfile(preProcDir,'run%s_f_transmid%s'); % requires extension and run number
trlfile  = fullfile(preProcDir, 'run%s_trl.mat'); % this contains our trial related info for epoching - should exist (generated with megTriggers script)
behavfile = fullfile(datadir,'behavioural',[thisSubject.id '_MEGRTs.mat']); % grab our behavioural data + output from a3_megtriggers

%% print subject for logging

disp('this is subject:')
disp(thisSubject.id);

for runNum = 1:numel(thisSubject.meg_runs)
    
    thisInFile = sprintf(inputFileName,num2str(runNum));
    if ica; ext = '.fif'; else; ext = '.mat'; end
    thisOutFile = sprintf(outputFileName,num2str(runNum),ext); clear ext
    
    if exist(thisInFile,'file')
        if ~exist(thisOutFile,'file') || overwrite
            
            fprintf('filtering %.0f of %.0f run files\n',runNum,numel(thisSubject.meg_runs))
            
            %% set up a fieldtrip config to load + check the data
            cfg = [];
            cfg.continuous = 'yes'; % this data is not epoched, it's continuous
            % let's get rid of all the junk channels
            % you wouldn't do this if you were planning to use e.g. the STI
            % channels to define triggers later (but I've already done that in a2_megTriggers)
            cfg.channel = {'EOG' 'ECG' 'EEG' 'MEG'};
            if endsWith(thisInFile,'.fif')
                % either load from fif file by specifying path
                cfg.dataset = thisInFile;
                rawData = ft_preprocessing(cfg); % load it
            elseif endsWith(thisInFile,'.mat')
                % or load from existing ft data object by loading the file and
                % passing the object into ft_preprocessing
                rawData = ft_preprocessing(cfg,loadFtData(thisInFile)); % load it
            end
            %% view the frequencies
            
            % you can use this to view: frequency on x, power spectrum on y
            % I have a big spike at 50 (line noise) for example and the odd
            % harmonics in particular (150, 250 etc)
            plotFreqs(rawData,'EEG')
            print([thisOutFile '_EEG_pre.png'], '-dpng');
            plotFreqs(rawData,'MEG')
            print([thisOutFile '_MEG_pre.png'], '-dpng');
            
            %% NOTE: if you looked at your data, be sure to re-run the previous ft config section before moving on
            
            %% band stop filter
            % so what we almost certainly want to do is filter any line noise (+
            % harmonics of the line noise if necessary). Typically, we do
            % this with a notch (a.k.a. bandstop) filter.
            %         cfg = [];
            %         cfg.bsfilter = 'yes';
            %         cfg.bsfreq = [49 51];
            %         filteredData = ft_preprocessing(cfg,rawData);
            % Obviously this means you lose any information at these freqs.
            % Also, this will introduce a symmetric smearing of the signal in
            % the time domain (aka ringing/transients)
            % so it's not really actually that clear to me why we'd do this
            % if you plot pre and post notch filter with my plotFreq()
            % function, you'll see and be dismayed by the smearing
            
            %% dft filter
            % fieldtrip has a specific kind of filter for this: DFT (a sharp
            % 'discrete fourier transform' but note this is sharpest with
            % lots of padding either side (5-10 secs) and also can have weird effects if the line noise
            % isn't constant (https://www.fieldtriptoolbox.org/faq/why_is_there_a_residual_50hz_line-noise_component_after_applying_a_dft_filter/)
            % so where typically you'd filter BEFORE epoching, (or at least
            % stitch the epochs together again before filtering), with this we
            % don't want long data segments for the line noise to fluctuate over
            % the way I've done it below doesn't entirely eliminate the line
            % noise (probably because it's fluctation over time as linked above)
            % but certainly substantially reduces it (and seems better
            % than the smearing from the notch filter!)
            
            % so first we cut the data up into arbitrary chunks
            cfg = [];
            cfg.length = 1; % 1 = 1 second epochs
            cfg.overlap = 0;
            segmentedRawData = ft_redefinetrial(cfg,rawData);
            % then apply the filter
            cfg = [];
            cfg.dftfilter = 'yes';
            cfg.dftfreq = [50 100 150];
            % for dftfreq [50 100 150] is the default, but we might be able to
            % improve performance by adding neighbouring frequencies to get rid
            % of residual line noise (e.g. if the line noise fluctuates)
            % the sampling of neighbouring frequencies will depend on the padding parameter
            % but I couldn't get any variation of this to really improve things much
            % if I was going to try, remember to pad with 5 or 10 seconds of
            % data
            %         cfg.padding = '10';
            % alternatively, we can try spectrum interpolation, which is a
            % newish feature (2019) and performs much better on my data,
            % although at the cost of introducing some ringing
            cfg.dftreplace = 'neighbour';
            cfg.dftbandwidth = [2 2 2];
            cfg.dftneighbourwidth = [2 2 2];
            segmentedFilteredData = ft_preprocessing(cfg,segmentedRawData);
            % in the end, I was left mightily unimpressed by all this, and now
            % very much understand why some people just don't bother
            % if you're doing decoding, you probably don't need to worry about
            % this too much, but if you're doing univariate analyis I'd be pretty
            % worried that this was making things worse
            
            % now we'll put that segmented data back together again in case we want to do other kinds of filtering:
            clear segmentedRawData
            cfg = [];
            cfg.continuous = 'yes';
            filteredData = ft_redefinetrial(cfg,segmentedFilteredData);
            clear segmentedFilteredData
            
            %% other filters
            
            % it's worth noting that ica could probably pick up line noise and
            % harmonics too
            
            % note all filters default to a two-pass filter (back and forward)
            % which results in a zero-phase shift of erp components, but you can
            % change this in fieldtrip---both order and direction
            
            % you might want other filters, but think about whether you should do
            % these pre or post ICA
            % you probably want a different padding for this than the dft filter, so you'd run it
            % seperately, clearing the cfg.dftfilter/freq and overridding the
            % padding, then running ft_preprocessing again.
            % you also probably want to stitch together the epochs for other
            % kinds of filters. I do this below for the pre-ICA file
            
            % we might want to high-pass filter to remove slow drift potentials/fields that
            % the faster signals of interest ride on. But this risks losing slow
            % fluctuations of brain potential, and introducing artefacts (ringing).
            % temporal maxfilter high-pass filters already, so you might not need
            % to do more of this.
            % I have done my maxfilter at a 10 sec data buffer, which corresponds to  a cut-off freq of 0.1Hz (1/10s)
            % I'm pretty sure my old script filtered again at 0.5, and lots of people
            % seem to prefer that, but I'd leave this unless your effect isn't
            % visible, and in that case, return and see if this improves things
            %         cfg.hpfilter = 'no';
            
            % we also might want to low pass here, the benefit being that it might
            % remove high-frequency variance that is often noise in relation to the
            % typically slower dynamics of interest. But this risks masking
            % high-frequency (e.g. gamma band) activity, or even just useful
            % information about the shape of some responses. I'm not going to
            % because I'm pretty sure I care about gamma.
            %         cfg.lpfilter = 'no';
            
            %% plot the results
            plotFreqs(filteredData,'EEG')
            print([thisOutFile '_EEG_post.png'], '-dpng');
            plotFreqs(filteredData,'MEG')
            print([thisOutFile '_MEG_post.png'], '-dpng');
            
            
            %% save it
            
            disp('saving');
            saveThis(thisOutFile, filteredData)
            
            disp('done')
            
        else
            disp(['the file:\n' thisOutFile '\nexists, skipping'])
        end % end outfile check
    else
        warning(['the file:\n' thisInFile '\ndoes not exist, skipping'])
    end % end infile check
end % end runs



end % end function

function figHandle = plotFreqs(data,channels)
% this basically expect you'll do 'EEG' or 'MEG' as channels, but it might
% work with a cell array of strings like {'EEG MEG'}

cfg = [];
cfg.output = 'pow';
cfg.channel = channels;
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.pad = 'nextpow2';
cfg.foilim = [0.3 data.fsample/2];
freqChk = ft_freqanalysis(cfg,data);
figure; plot(freqChk.freq,mag2db(freqChk.powspctrm));
figHandle = gcf;

end

function saveThis(savePath, data)
% requires the save directory and ft on the path

disp('saving data:')
disp(data)
disp('to:')
disp(savePath)
if endsWith(savePath,'fif')
    ft_defaults;
    fieldtrip2fiff(savePath,data);
else
    save(savePath,'data','-v7.3')
end
disp('done')

end
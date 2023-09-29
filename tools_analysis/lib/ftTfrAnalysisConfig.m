function cfg = ftTfrAnalysisConfig(lockedTo,trials,freqs)

% first let's create a basic config
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.pad          = 'nextpow2'; % ft recommends nextpow2 over their default 'maxperlen' as more efficient for FFT computation
cfg.method       = 'mtmconvol';

switch freqs
    case 'high'
        % multitaper good for high freqs (30Hz plus). see: https://www.fieldtriptoolbox.org/tutorial/timefrequencyanalysis/#multitapers
        % it is the default for mtmconvol: mtm = multitaper method, so no
        % need to set cfg.taper
        % see also:
        % https://mailman.science.ru.nl/pipermail/fieldtrip/2007-August/001327.html
        % for suggested settings
        % t_ftimwin of 0.25 s
        % -> Rayleigh frequency of 4 Hz
        % Spectral concentration (tapsmofrq) of plus/minus 12 Hz (up to 20 Hz)
        cfg.foi          = 30:4:150; % do 1 to 30 Hz in steps of 4 Hz
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.25; % length of time window = 0.25 sec
        cfg.tapsmofrq    = 12; % width of frequency smoothing in Hz
        % 
    case 'low'
        % hanning is good for low freqs, up to gamma (30Hz)
        % we will use the fieldtrip tutorial settings: https://www.fieldtriptoolbox.org/tutorial/timefrequencyanalysis/#hanning-taper-fixed-window-length
        cfg.taper        = 'hanning';
        cfg.foi          = 2:2:30; % do 1 to 30 Hz in steps of 2 Hz
        % t_ftimewin: 1/length of time window in sec must equal cfg.foi Hz steps. an integer number of cycles must fit in the time window.
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5; % length of time window = 0.5 sec
end

% if we want to filter by trials
if exist('trials','var')
    cfg.trials = trials;
end

% a timewindow that slides from the start of the period of interest to the end
if strcmp(lockedTo,'response')
    cfg.toi          = -0.6:0.05:0.2;
elseif strcmp(lockedTo,'onset')
    cfg.toi          = -0.5:0.05:1.5;
end

return
end
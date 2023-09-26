function cfg = ftTfrAnalysisConfig(lockedTo,trials)

% first let's create a basic config
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
% cfg.pad          = 'nextpow2'; % ft recommends nextpow2 over their default 'maxperlen' as more efficient for FFT computation
cfg.foi          = 1:2:100; % do 2 to 100 Hz in steps of 2 Hz
% t_ftimewin: 1/length of time window in sec must equal cfg.foi Hz steps. an integer number of cycles must fit in the time window.
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.1; % length of time window = 0.5 sec

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
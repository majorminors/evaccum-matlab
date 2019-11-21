%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);

% set up variables
rootdir = '\\cbsu\data\Group\Woolgar-Lab\projects\EvAccum'; % root directory - used to inform directory mappings
p.resp_keys = {'RB','RY'}; % only accepts two response options
p.continue_key = {'RR'};
p.MEG_enabled = 1;
p.quitkey = 'q';

% directory mapping
addpath(genpath(fullfile(rootdir, 'tools'))); % add tools folder to path (includes moving_dots function which is required for dot motion, as well as an external copy of subfunctions for backwards compatibility with MATLAB)
save_file = fullfile(rootdir, mfilename);
save(save_file);

% psychtoolbox setup
AssertOpenGL; % check Psychtoolbox (on OpenGL) and Screen() is working
KbName('UnifyKeyNames'); % makes key mappings compatible (mac/win)

% set up MEG
MEG = MEGSynchClass; % call MEG function
MEG.SendTrigger(0); % make sure all triggers are off

try
    
    fprintf('running %s\n', mfilename);
    
    t.queuekeys = KbName(p.quitkey); % define the keys the queue cares about
    t.queuekeylist = zeros(1,256); % create a list of all possible keys (all 'turned off' i.e. zeroes)
    t.queuekeylist(t.queuekeys) = 1; % 'turn on' the keys we care about in the list (make them ones)
    KbQueueCreate([], t.queuekeylist); % initialises queue to collect response information from the list we made (not listening for response yet)
    KbQueueStart(); % starts delivering keypress info to the queue
    
    fprintf('entering response waiter function\n');
    response_waiter(p,MEG)
    KbQueueFlush(); % flush the response queue from the response waiter
    
    fprintf('entering megtest function\n');
    KbQueueFlush(); % flush the response queue so any accidental presses recorded in the cue period won't affect responses in the dots period
    [t.pressed, t.firstPress] = megtest_function(p,MEG); % pull time of first flip for dots, as well as information from KBQueueCheck from moving_dots
    
    MEG.LastButtonPress
    t.firstPress{1}
    MEG.TimeOfLastButtonPress
    t.firstPress{2}
    
    save(save_file);
    
catch err
    save(save_file);
    %MEG.delete;
    rethrow(err);

end
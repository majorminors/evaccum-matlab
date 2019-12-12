%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);

% set up variables
rootdir = '\\cbsu\data\Group\Woolgar-Lab\projects\EvAccum'; % root directory - used to inform directory mappings
p.resp_keys = {'a','s'}; % only accepts two response options
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
MEG = MEGSynchClass(1); % call MEG function
MEG.SendTrigger(0); % make sure all triggers are off
MEG.Keys = {'a','s'};

try
    
    fprintf('running %s\n', mfilename);
    
%     t.queuekeys = KbName(p.quitkey); % define the keys the queue cares about
%     t.queuekeylist = zeros(1,256); % create a list of all possible keys (all 'turned off' i.e. zeroes)
%     t.queuekeylist(t.queuekeys) = 1; % 'turn on' the keys we care about in the list (make them ones)
%     KbQueueCreate([], t.queuekeylist); % initialises queue to collect response information from the list we made (not listening for response yet)
%     KbQueueStart(); % starts delivering keypress info to the queue
    
    fprintf('entering response waiter function\n');
    if p.MEG_enabled == 1
        MEG.WaitForButtonPress(0); % reset MEG button press to empty
    end
    response_waiter(p,MEG)
%     KbQueueFlush(); % flush the response queue from the response waiter
    
    fprintf('entering megtest function\n');
    if p.MEG_enabled == 1
        MEG.WaitForButtonPress(0); % reset MEG button press to empty
    end
%     KbQueueFlush(); % flush the response queue so any accidental presses recorded in the cue period won't affect responses in the dots period
    
    %% this is where the function would be
if p.MEG_enabled == 1
    MEG.SendTrigger(1); % send a trigger for trial onset
    MEG.ResetClock; % reset the timer
    button_pressed = 0; % a counter to make sure we catch the first time a button was pressed
    MEG.WaitForButtonPress(1); % listen for button press
    pause(0.005); % quick pause before we reset triggers
    MEG.SendTrigger(0); % reset triggers
end

fprintf('entering function response loop\n');
               WaitSecs(1)
for n = 1:10
n

    fprintf('button coding block begins\n');
    if ~isempty(MEG.LastButtonPress) && ~button_pressed % check for a keypress in the MEG key wait function every frame, if a key hasn't been pressed yet
        button_pressed = 0; % record that a key has been pressed this trial
        MEG.SendTrigger(2); % send a trigger
        firstPress{1} = MEG.LastButtonPress; % record the key pressed
        firstPress{2} = MEG.TimeOfLastButtonPress; % record the time of the key pressed
        fprintf('recorded response as\n'); firstPress{1}
        pressed = NaN; % just put something in here
        pause(0.005); % quick pause before resetting
        MEG.SendTrigger(0); % reset the triggers
        WaitSecs(1)
%         fpidx = 0;
%     elseif ~strcmp(MEG.LastButtonPress,firstPress{1}) && button_pressed > 0 % if it's not the first time a key has been pressed
%         button_pressed = 2; % record that another key has been pressed this trial
%         fpidx = fpidx+1;
%         firstPress{3+fpidx} = MEG.LastButtonPress; % record the key pressed
%         firstPress{4+fpidx} = MEG.TimeOfLastButtonPress; % record the time of the key pressed
%         fprintf('recorded a subsequent response as\n'); firstPress{3+fpidx}
    end
    
end
    
    %% function would end here
    
    MEG.LastButtonPress
    firstPress{1}
    MEG.TimeOfLastButtonPress
    firstPress{2}
    
    save(save_file);
    
catch err
    save(save_file);
    %MEG.delete;
    rethrow(err);

end
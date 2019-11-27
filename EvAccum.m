%% Matching motion coherence to direction cue in MEG
% Dorian Minors
% Created: JUN19
% Last Edit: 18NOV19

% trial settings all saved in 'p'

% data related information saved in 'd'
%   d.stim_mat_all contains trial condition matrices for each block
%       (d.stim_mat_all(:,:,[block number])
%   rows indicate trials
%   columns are explained when we define the matrix in 'define stimulus
%       parameters', below
%   in this script two blocks represent two different response
%       mappings (to counterbalance two button-presses)

% other trial specific variables are in 't' in case something goes wrong
%   and we want to see them

% note: response keys for the trial based off whether participant id is odd or even

% note: will swap response keys according to p.keyswap

% note: will provide breaks before blocks specified in p.breakblocks

% note: d.correct will record correct (1), incorrect (0), and nil (-1) responses.
%       a nil response will occur if participant responds in cue period and fails
%       to respond again in the dots presentation period.

% note: for backwards compatibility (to 2016), subfunctions are disabled
%       in scripts (although functions can have subfunctions and scripts
%       instead use external functions.

% note: for compatibility with testing lab computers, KbQueueWait() only accepts
%       one argument - but from at least 2018, KbQueueWait can accept two,
%       which is commented out next to any call to KbQueueWait here.

%% subclasses

%% MEGSynchClass
% a class of functions for the MEG interface with the National Instruments PCI 6503 card (MRC CBU)

%% subfunctions - can use these at bottom instead of external functions from at least R2018a:

%% dots_onset_time, pressed, firstPress] = moving_dots(p,dots,MEG,exp_trial)
% creates a cloud of moving dots, then creates a fixation by putting a small
% black square in a bigger white square, then flips the screen
%% response_waiter(p,MEG)
% will wait for button responses before continuing
% updated 18NOV2018
% since previous versions of MATLAB require functions to be external to
% script, then if you're using this, check it's up to date with the
% external one
%% f = keyswap_training(p,dots,d,MEG)
% will run some training on the EvAccum paradigm in a sandbox when called
%% pix = angle2pix(p,ang)
% calculates pixel size from visual angles, assuming isotropic (square) pixels
%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);
p = struct(); % est structure for parameter values
d = struct(); % est structure for trial data
t = struct(); % another structure for untidy trial specific floating variables that we might want to interrogate later if we mess up

% set up variables
rootdir = 'C:\Users\doria\Google Drive\04 Research\05 Evidence Accumulation\01 EvAccum Code';%'\\cbsu\data\Group\Woolgar-Lab\projects\EvAccum'; % root directory - used to inform directory mappings
p.screen_num = 0; % screen to display experiment on (0 unless multiple screens)
p.fullscreen_enabled = 0; % 1 is full screen, 0 is whatever you've set p.window_size to
p.testing_enabled = 1; % change to 0 if not testing (1 skips PTB synctests and sets number of trials and blocks to test values) - see '% test variables' below
p.training_enabled = 0; % if 0 (or any other than 1) will do nothing, if 1, initiates training protocol (reduce dots presentation time from 'p.training_dots_duration' to 'p.dots_duration' by one 'p.training_reduction' every 'p.training_interval') - see '% training variables' below
p.fix_trial_time = 1; % if 0 then trial will end on keypress, if 1 will go for duration of p.dots_duration
p.iti_on = 1; % if 1 will do an intertrial interval with fixation, if 0 (or anything other than 1) will not do iti
p.feedback_type = 1; % if 0 (or anything other than 1 or 2) no feedback, if 1 then trialwise feedback, if 2 then blockwise feedback
p.num_blocks = 20;
p.breakblocks = 0; % before which blocks should we initiate a break (0 for no breaks, otherwise to manipulate based on a fraction of blocks, use 'p.num_blocks' or if testing 'p.num_test_blocks')
p.keyswap = 2; % swaps keys at some point in experiment - 1 to not swap, 2 to swap once, 3 to swap twice etc (it's a division operation)
p.MEG_enabled = 0; % using MEG
p.MEG_emulator_enabled = 1; % using the emulator - be aware we can't quit using the quitkey with emulator

% check set up
if ~ismember(p.fullscreen_enabled,[0,1]); error('invalid value for p.fullscreen_enabled'); end % check if valid or error
if ~ismember(p.testing_enabled,[0,1]); error('invalid value for p.testing_enabled'); end % check p.testing_enabled is a valid number, or error
if ~ismember(p.training_enabled,[0,1]); error('invalid value for p.training_enabled'); end % check if valid or error
if ~ismember(p.fix_trial_time,[0,1]); error('invalid value for p.fix_trial_time'); end % check if valid or error
if ~ismember(p.iti_on,[0,1]); error('invalid value for p.iti_on'); end % check if valid or error
if ~ismember(p.feedback_type,[0,1,2]); error('invalid value for p.feedback_type'); end % check if valid or error
if ~ismember(p.MEG_enabled,[0,1]); error('invalid value for p.MEG_enabled'); end % check if valid or error
if ~ismember(p.MEG_emulator_enabled,[0,1]); error('invalid value for p.MEG_emulator_enabled'); end % check if valid or error
%if p.MEG_enabled == 1 && p.testing_enabled == 1; error('are you sure you want to be testing with MEG enabled? if so, comment out this line'); end
if p.MEG_enabled == 1 && p.training_enabled == 1; error('you cannot train with MEG enabled currently'); end
if p.MEG_emulator_enabled == 1 && p.MEG_enabled == 0
    warning('you cannot emulate MEG without enabling MEG - turning off emulation\n');
    WaitSecs(1);
    p.MEG_emulator_enabled = 0;
    fprintf('p.MEG_emulator_enabled is now off\n')
end

% directory mapping
addpath(genpath(fullfile(rootdir, 'tools'))); % add tools folder to path (includes moving_dots function which is required for dot motion, as well as an external copy of subfunctions for backwards compatibility with MATLAB)
stimdir = fullfile(rootdir, 'stimuli');
datadir = fullfile(rootdir, 'data'); % will make a data directory if none exists
if ~exist(datadir,'dir')
    mkdir(datadir);
end

% training variables
p.training_dots_duration = 5; % duration (secs) dot stimulus presented for for training trials (will reduce trial by trial to p.dots_duration)
p.training_reduction = 1; % by how much do we reduce the duration during training (seconds)
p.training_interval = 2; % how many trials should we train on before reducing the dots presentation time

% test variables
p.num_test_trials = 3;
p.num_test_blocks = 4;
if p.testing_enabled == 1
    p.PTBsynctests = 1; % PTB will skip synctests if 1
    p.PTBverbosity = 1; % PTB will only display critical warnings with 1
elseif p.testing_enabled == 0
    p.PTBsynctests = 0;
    p.PTBverbosity = 3; % default verbosity for PTB
end
Screen('Preference', 'SkipSyncTests', p.PTBsynctests);
Screen('Preference', 'Verbosity', p.PTBverbosity);

% psychtoolbox setup
AssertOpenGL; % check Psychtoolbox (on OpenGL) and Screen() is working
KbName('UnifyKeyNames'); % makes key mappings compatible (mac/win)
rng('shuffle'); % seed rng using date and time

% set up participant info and save
t.prompt = {'enter participant number:',... % prompt a dialog to enter subject info
    'enter easy coherence threshold(fm 0-1, higher is easier)',...
    'enter hard coherence threshold (fm 0-1, lower is harder)',...
    'enter hard matching threshold (between 0 and 90 degrees from cued direction)'};%',...
% 'enter easy matching threshold (between 0 and 90 degrees from cued direction)',...
% 'enter hard matching threshold (between 0 and 90 degrees from cued direction)'};
t.prompt_defaultans = {num2str(99), num2str(0.75), num2str(0.25), num2str(80)}; % default answers corresponding to prompts
t.prompt_rsp = inputdlg(t.prompt, 'enter participant info', 1, t.prompt_defaultans); % save dialog responses
d.participant_id = str2double(t.prompt_rsp{1}); % add subject number to 'd'
d.easy_coherence = str2double(t.prompt_rsp{2}); % add participant coherence thresholds to 'd'
d.hard_coherence = str2double(t.prompt_rsp{3}); % add participant coherence thresholds to 'd'
d.hard_rule = str2double(t.prompt_rsp{4}); % add participant matching thresholds to 'd'
d.easy_rule = 90-d.hard_rule; % add participant matching thresholds to 'd' - in this case, the easy rule is the inverse of the hard rule
% d.easy_rule = str2double(t.prompt_rsp{4}); % add participant matching thresholds to 'd'
% d.hard_rule = str2double(t.prompt_rsp{5}); % add participant matching thresholds to 'd'

% check participant info has been entered correctly for the script
if isnan(d.participant_id)
    error('no participant number entered');
elseif isnan(d.easy_coherence) || d.easy_coherence > 1 || d.easy_coherence < 0
    error('invalid participant coherence threshold (easy)');
elseif isnan(d.hard_coherence) || d.hard_coherence > 1 || d.easy_coherence < 0
    error('invalid participant coherence threshold (hard)')
elseif isnan(d.easy_rule) || d.easy_rule > 90 || d.easy_rule < 0
    error('invalid participant matching threshold (easy)');
elseif isnan(d.hard_rule) || d.hard_rule > 90 || d.hard_rule < 0
    error('invalid participant matching threshold (hard)')
end

% create a save file
save_file_name = [num2str(d.participant_id,'S%02d'),'_',mfilename];
save_file = fullfile(datadir, save_file_name);
if exist([save_file '.mat'],'file') % check if the file already exists and throw a warning if it does
    warning('the following save file already exists - overwrite? (y/n)\n %s.mat', save_file);
    while 1 % loop forever until y or n
        ListenChar(2);
        [secs,keyCode] = KbWait; % wait for response
        key_name = KbName(keyCode); % find out name of key that was pressed
        if strcmp(key_name, 'y')
            fprintf('instructed to overwrite:\n %s.mat\n overwriting and continuing with %s\n', save_file, mfilename)
            ListenChar(0);
            clear secs keyCode key_name
            break % break the loop and continue
        elseif strcmp(key_name, 'n')
            ListenChar(0);
            clear secs keyCode key_name
            error('instructed not to overwrite:\n %s.mat\n aborting %s\n', save_file, mfilename); % error out
        end
    end % end response loop
end % end check save file exist
save(save_file); % save all data to a .mat file

%% define experiment parameters

fprintf('defining exp params for %s\n', mfilename);

% define keys
if p.MEG_enabled == 0
    p.resp_keys = {'a','s'}; % only accepts two response options
elseif p.MEG_enabled == 1 % what keys in the MEG
    if p.MEG_emulator_enabled == 0
        p.resp_keys = {'RB','RY'}; % only accepts two response options
    elseif p.MEG_emulator_enabled == 1
        p.resp_keys = {'LY','RB'}; % LY and RB correspond to a and s on the keyboard
    end
    %p.continue_key = {'RR'};
end
p.quitkey = {'q'}; % this is watched by KbQueue regardless of p.MEG_enabled
t.keyswapper = 0; % will use this variable to mark a keyswap event (code currently at commencement of block loop)
% establish response keys for the trial based off whether participant id is odd or even
if mod(d.participant_id,2) % if pid not divisible by 2 (i.e. leaves a modulus after division)
    p.resp_keys = {p.resp_keys{1},p.resp_keys{2}}; % essentially do nothing - keep resp keys the same
elseif ~mod(d.participant_id,2) % if pid is divisible by 2 (i.e. does not leave a modulus after division)
    p.resp_keys = {p.resp_keys{2},p.resp_keys{1}}; % then swap the response keys
end


% define display info
p.bg_colour = [0 0 0]; % needs to be the same as the cue stimuli background colour (unless transparent)
p.text_colour = [255 255 255]; % colour of instructional text
p.matching_cue_1 = 'BLUE'; % variable used to indicate response keys - this the upward arrow of the doublesided arrow cue in stimdir
p.matching_cue_2 = 'ORANGE'; % variable used to indicate response keys - this the downward arrow of the doublesided arrow cue in stimdir
p.cue_colour_blue = [121 181 240]; % colour of cue, for text formatting
p.cue_colour_orange = [240 181 121]; % colour of cue, for text formatting
p.text_size = 40; % size of text
p.window_size = [0 0 1200 800]; % size of window when ~p.fullscreen_enabled
p.screen_width = 35;   % Screen width in cm
p.screen_height = 50;    % Screen height in cm
p.screen_distance = 50; % Screen distance from participant in cm
p.visual_angle_cue = 15; % visual angle of the cue expressed as a decimal - determines size
p.visual_angle_dots = 0.15; % visual angle of the dots expressed as a decimal - determines size

% timing info
p.min_cue_time = 0.5; % minimum period to display cue (participants can't continue during this time)
p.iti_time = 0.3; % inter trial inteval time
p.MEG_onset_trigger_time = 0.1; % time to let the MEG trigger reach full strength - recommend 100ms but must be smaller than p.iti_time
p.dots_duration = 1.5; % seconds for the dot cloud to be displayed
p.feedback_time = 0.5; % period to display feedback after response
p.keyswap_inform_time = 1; % minumum period to display keyswap notification
p.break_inform_time = 1; % minumum period to display break notification (stop participants from accidentally continuing)
if p.MEG_onset_trigger_time>p.iti_time; error('p.MEG_onset_trigger_time is larger than p.iti_time'); end

% trial settings (*p.stim_mat* = parameter required to calculate stimulus condition matrix)
t.takeabreak = 0; % will use this variable to mark a break event (code currently at commencement of block loop)
p.num_trials_per_block = 64; % *p.stim_mat* - as many as unique conditions
p.num_cues = 4; % *p.stim_mat*
p.num_motion_coherence = 8; % *p.stim_mat* - number of coherent directions
p.cue_directions = 45:90:315; % *p.stim_mat* - refers to the direction of the upward arrow of the doublesided arrow cue in stimdir
p.dot_motion_directions = union([p.cue_directions+d.easy_rule],[p.cue_directions+d.hard_rule]); % *p.stim_mat* - adds easy rule and hard rule to each cue, then puts them in a vector sorted low to high

% lets check all those parameters
t.view_p = struct2table(p, 'AsArray', true);
disp(t.view_p);
warning('happy with all this? (y/n)\n %s.mat', save_file);
while 1 % loop forever until y or n
    ListenChar(2);
    [secs,keyCode] = KbWait; % wait for response
    key_name = KbName(keyCode); % find out name of key that was pressed
    if strcmp(key_name, 'y')
        fprintf('happy with parameters\n continuing with %s\n', mfilename)
        ListenChar(0);
        clear secs keyCode key_name
        break % break the loop and continue
    elseif strcmp(key_name, 'n')
        ListenChar(0);
        clear secs keyCode key_name
        error('not happy with parameters\n aborting %s\n', mfilename); % error out
    end
end % end response loop
%% define stimuli parameters

fprintf('defining stimuli params for %s\n', mfilename);


% read in stimulus file for the cue
p.cue = imread(fullfile(stimdir, 'arrows_cue_colours_with_line.png'));

% dot cloud parameters for moving_dots function
%   note that for the following required params for moving_dots:
%       'dots.direction' is specified within the trial loop using outputs from
%           'p.stim_mat'
%       'dots.coherence' is also specified withing the trial loop using outputs
%           from 'p.stim_mat' and 'd.easy_coherence'/'d.hard_coherence'
dots.aperture_size = [10,10]; % [width,height] in degrees of the rectangular aperture dots are displayed in
dots.centre = [0,0]; % [x,y] centre of the dot cloud
dots.num_dots = 100; % number of dots
dots.colour = [255,255,255]; % colour of the dots in [r,g,b]
dots.visual_angle = p.visual_angle_dots; % visual angle of the dots expressed as a decimal - determines size
dots.speed = 5; % speed of dots in degrees per second
dots.lifetime = 5; % number of frames dots live for


% create matrix specifying stimulus conditions per trial:
%  1)  cue direction (1-4)
%  2)  cue direction in degrees
%  3)  dot motion direction condition (1-8)
%  4)  dot motion direction in degrees
%  5)  coherence difficulty (1 = easy, 2 = hard)
%  6)  matching distance from cue direction (absolute value) - used to calc match and match difficulty
%  7)  match blue (1) or match orange (2)
%  8)  matching difficulty (1 = easy, 2 = difficult)
%  9)  gives you a unique number for each trial condition
% 10)  gives a number based on 9 to identify each response for each trial condition
% note: if you try to test with two matching angles that are the same, you
%       will get an error, so make them different by at least 1 degree. this is
%       because we use 'union()' to calc matching distance from cue direction.

% create matrix
p.stim_mat = zeros(p.num_trials_per_block,10);

% create columns
p.stim_mat(:,1) = sort(repmat(1:p.num_cues,[1,p.num_trials_per_block/p.num_cues]));
p.stim_mat(:,2) = p.cue_directions(p.stim_mat(:,1));
p.stim_mat(:,3) = repmat(sort(repmat(1:p.num_motion_coherence,[1,p.num_trials_per_block/p.num_cues/p.num_motion_coherence])),[1,p.num_cues]);
p.stim_mat(:,4) = p.dot_motion_directions(p.stim_mat(:,3));
p.stim_mat(:,5) = repmat(sort(repmat(1:2,[1,p.num_trials_per_block/p.num_cues/p.num_motion_coherence/2])),[1,p.num_cues*p.num_motion_coherence]);
dist = abs([p.stim_mat(:,4)-360,p.stim_mat(:,4),p.stim_mat(:,4)+360]-repmat(p.stim_mat(:,2),[1,3]));
p.stim_mat(:,6) = min(dist,[],2);
p.stim_mat(:,7) = (p.stim_mat(:,6)>90)+1;
p.stim_mat(:,8) = ~((p.stim_mat(:,6)==min(p.stim_mat(:,6)))|(p.stim_mat(:,6)==max(p.stim_mat(:,6))))+1;
p.stim_mat(:,9) = 1:length(p.stim_mat(:,9));
p.stim_mat(:,10) = p.stim_mat(:,9)+p.num_trials_per_block;
% clear floating variables
clear dist;

% shuffle the trial condition order for each block into 'd.stim_mat_all' then re-sort by cue
for block=1:p.num_blocks
    d.stim_mat_all(:,:,block) = p.stim_mat(Shuffle(1:p.num_trials_per_block),:);
    d.stim_mat_all(:,:,block) = sortrows(d.stim_mat_all(:,:,block),1);
end
% clear floating variables
clear block;

%% set up MEG

% MEG trigger info
p.MEGtriggers.training = 255; % unique trigger to tell us when to ignore triggers sent during training
p.MEGtriggers.onsets = 9; % what column of p.stim_mat are you keeping your trigger information for onset in?
p.MEGtriggers.responses = 10; % what column p.stim_mat are you keeping your trigger information for responses in?

% invoke the MEG functions if p.MEG_enabled
if p.MEG_enabled == 1
    if p.MEG_emulator_enabled == 0
        MEG = MEGSynchClass;
    elseif p.MEG_emulator_enabled == 1
        MEG = MEGSynchClass(1); % call MEG function
        MEG.Keys = {'a','s'};
    end
    MEG.SendTrigger(0); % make sure all triggers are off
elseif p.MEG_enabled == 0
    MEG = 0; % just put something here so the switches and the functions we pass this to don't freak out
end

%% exp start

fprintf('running experiment %s\n', mfilename);

try
    
    % open screen
    if p.fullscreen_enabled % zero out p.window_size if p.fullscreen_enabled = 1
        p.window_size=[];
    end
    [p.win,p.rect] = Screen('OpenWindow',p.screen_num,p.bg_colour,p.window_size);
    Screen('BlendFunction',p.win,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA); % allows transparency in .png images
    Screen('TextSize', p.win, p.text_size); % set the text size
    % then need some info based on the screen for later
    p.frame_rate = 1/Screen('GetFlipInterval', p.win); % is Hz - for moving_dots
    p.resolution = p.rect([3,4]); % pull resolution info from p.rect - used to scale cue image and is passed to moving_dots to do the same
    HideCursor;
    WaitSecs(0.5); % warm up
    
    %% block start
    
    % changes number of blocks to testing amount if testing
    if p.testing_enabled == 1
        p.act_block_num = p.num_test_blocks;
        fprintf('testing (p.testing_enabled set to 1) - will run %u blocks\n', p.num_test_blocks);
    else
        p.act_block_num = p.num_blocks;
    end
    
    % repeat trials for each block
    for block = 1:p.act_block_num
        fprintf('entering block %u of %u\n',block, p.act_block_num); %report block number to command window
        
        % pick up trial condition order for this block
        p.stim_mat = d.stim_mat_all(:,:,block);
        
        % initiate break if block is a break block
        %t.breakblock = ismember(block,p.breakblocks);
        if ismember(block,p.breakblocks)
            t.takeabreak = 1; % break initiated (tells a screen to pop up)
        end
        
        % swap response keys halfway through blocks (rounded to nearest integer)
        if block == round(p.act_block_num/p.keyswap)+1
            p.resp_keys = {p.resp_keys{2},p.resp_keys{1}}; % then swap the response keys
            t.keyswapper = 1; % keyswap initiated (tells a screen to pop up informing participants and starts some training)
        end
        
        %% trials start
        
        % changes number of trials to testing amount for trial loop if testing
        if p.testing_enabled == 1
            p.act_trial_num = p.num_test_trials;
            fprintf('testing (p.testing_enabled set to 1) - will only run %u trials\n', p.num_test_trials);
        else
            p.act_trial_num = p.num_trials_per_block;
        end
        
        % open trial loop
        i = 0;
        while i < p.act_trial_num
            i = i + 1;
            fprintf('trial %u of %u\n',i,p.act_trial_num); % report trial number to command window
            
            %set up a queue to collect response info (or in the case of p.MEG_enabled, to watch for the quitkey)
            if ~p.MEG_emulator_enabled
                if p.MEG_enabled == 0
                    t.queuekeys = [KbName(p.resp_keys{1}), KbName(p.resp_keys{2}), KbName(p.quitkey)]; % define the keys the queue cares about
                elseif p.MEG_enabled == 1
                    t.queuekeys = KbName(p.quitkey); % define the keys the queue cares about
                end
                t.queuekeylist = zeros(1,256); % create a list of all possible keys (all 'turned off' i.e. zeroes)
                t.queuekeylist(t.queuekeys) = 1; % 'turn on' the keys we care about in the list (make them ones)
                KbQueueCreate([], t.queuekeylist); % initialises queue to collect response information from the list we made (not listening for response yet)
                KbQueueStart(); % starts delivering keypress info to the queue
            end
            
            % take a break
            if t.takeabreak == 1
                DrawFormattedText(p.win, sprintf('\n take a little break\n\n we are on block %u of %u\n',block, p.act_block_num), 'center', 'center', p.text_colour);
                Screen('Flip', p.win);
                WaitSecs(p.break_inform_time);
                DrawFormattedText(p.win,sprintf('\n take a little break\n\n we are on block %u of %u\n\n experimenter will continue',block, p.act_block_num), 'center', 'center', p.text_colour);
                Screen('Flip', p.win);
                response_waiter(p,MEG) % call response_waiter function
                if ~p.MEG_emulator_enabled; KbQueueFlush(); end % flush the response queue from the response waiter
                t.takeabreak = 2; % break event complete
            end
            
            % do keyswap and training
            if t.keyswapper == 1
                DrawFormattedText(p.win,'\n the response keys have swapped!\n\n\n we will do some practice\n\n\n ', 'center', 'center', p.text_colour);
                Screen('Flip', p.win);
                WaitSecs(p.keyswap_inform_time);
                DrawFormattedText(p.win,'\n the response keys have swapped!\n\n\n we will do some practice\n\n\n press either button to continue\n', 'center', 'center', p.text_colour);
                Screen('Flip', p.win);
                response_waiter(p,MEG) % call response_waiter function
                if ~p.MEG_emulator_enabled; KbQueueFlush(); end % flush the response queue from the response waiter
                t.keyswapper = 2; % keyswap complete
                fprintf('first we will run some practice on the new keys before we get into trial %u\n', i); % report that we're about to do some training
                if p.MEG_enabled == 1
                    MEG.SendTrigger(p.MEGtriggers.training); % send a trigger to inform we're training
                    pause(p.MEG_onset_trigger_time); % quick pause before we reset triggers
                    MEG.SendTrigger(0); % reset triggers
                end
                %% training function
                d.trainingdata = keyswap_training(p,dots,d,MEG); % run some training on the new keys
                %% training function ends
                if p.MEG_enabled == 1
                    MEG.SendTrigger(p.MEGtriggers.training); % send a trigger to inform we're training
                    pause(p.MEG_onset_trigger_time); % quick pause before we reset triggers
                    MEG.SendTrigger(0); % reset triggers
                end
            end
            
            %% present cue and response mapping
            if i == 1 || p.stim_mat(i,1) ~= p.stim_mat(i-1,1) || ~mod(i-1,8) % && p.stim_mat(i+1,1) == p.stim_mat(i,1) % if first trial, or cue changes (as currently blocked), or every 8 trials unless we're about to change cue then display cue
                
                % make the texture and scale it
                p.cue_tex = Screen('MakeTexture', p.win, p.cue);
                [t.tex_size1, t.tex_size2, t.tex_size3] = size(p.cue_tex); % get size of texture
                t.aspectratio = t.tex_size2/t.tex_size1; % get the aspect ratio of the image for scaling purposes
                t.imageheight = angle2pix(p,p.visual_angle_cue); % scale the height of the image using the desired visual angle
                t.imagewidth = t.imageheight .* t.aspectratio; % get the scaled width, constrained by the aspect ratio of the image
                
                % parameterise the rect to display cue in
                t.imgrect = [0 0 t.imagewidth t.imageheight]; % make a scaled rect for the cue
                t.rect = CenterRectOnPointd(t.imgrect,p.resolution(1,1)/2,p.resolution(1,2)/2); % offset it for the centre of the window
                
                % then display cue and response mapping
                Screen('DrawTexture', p.win, p.cue_tex, [], t.rect, p.stim_mat(i,2)); % draws the cue in the orientation specified in column 2 of p.stim_mat for the current trial
                % DrawFormattedText(p.win, sprintf('\n press %s for %s\n press %s for %s\n', p.resp_keys{1}, p.matching_cue_1, p.resp_keys{2}, p.matching_cue_2), 'center', t.rect(RectBottom)*1.01, p.text_colour);
                % draw the response mapping at the corners of the t.rect you just made depending on the orientation of the cue
                if p.stim_mat(i,1) == 1
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{1}), t.rect(RectRight)*1.01, t.rect(RectTop)*1.01, p.cue_colour_blue);
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{2}), t.rect(RectLeft)*1.01, t.rect(RectBottom)*1.01, p.cue_colour_orange);
                elseif p.stim_mat(i,1) == 2
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{1}), t.rect(RectRight)*1.01, t.rect(RectBottom)*1.01, p.cue_colour_blue);
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{2}), t.rect(RectLeft)*1.01, t.rect(RectTop)*1.01, p.cue_colour_orange);
                elseif p.stim_mat(i,1) == 3
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{1}), t.rect(RectLeft)*1.01, t.rect(RectBottom)*1.01, p.cue_colour_blue);
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{2}), t.rect(RectRight)*1.01, t.rect(RectTop)*1.01, p.cue_colour_orange);
                elseif p.stim_mat(i,1) == 4
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{1}), t.rect(RectLeft)*1.01, t.rect(RectTop)*1.01, p.cue_colour_blue);
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{2}), t.rect(RectRight)*1.01, t.rect(RectBottom)*1.01, p.cue_colour_orange);
                end
                d.cue_onset(block,i) = Screen('Flip', p.win); % pull the time of the screen flip from the flip function while flipping
                WaitSecs(p.min_cue_time);
                Screen('DrawTexture', p.win, p.cue_tex, [], t.rect, p.stim_mat(i,2)); % redraws the cue
                % redraw also the response mapping
                DrawFormattedText(p.win, sprintf('\n press either button to continue\n'), 'center', t.rect(RectBottom)*1.1, p.text_colour);
                if p.stim_mat(i,1) == 1
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{1}), t.rect(RectRight)*1.01, t.rect(RectTop)*1.01, p.cue_colour_blue);
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{2}), t.rect(RectLeft)*1.01, t.rect(RectBottom)*1.01, p.cue_colour_orange);
                elseif p.stim_mat(i,1) == 2
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{1}), t.rect(RectRight)*1.01, t.rect(RectBottom)*1.01, p.cue_colour_blue);
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{2}), t.rect(RectLeft)*1.01, t.rect(RectTop)*1.01, p.cue_colour_orange);
                elseif p.stim_mat(i,1) == 3
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{1}), t.rect(RectLeft)*1.01, t.rect(RectBottom)*1.01, p.cue_colour_blue);
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{2}), t.rect(RectRight)*1.01, t.rect(RectTop)*1.01, p.cue_colour_orange);
                elseif p.stim_mat(i,1) == 4
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{1}), t.rect(RectLeft)*1.01, t.rect(RectTop)*1.01, p.cue_colour_blue);
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{2}), t.rect(RectRight)*1.01, t.rect(RectBottom)*1.01, p.cue_colour_orange);
                end
                Screen('Flip', p.win);
                %% response waiter function
                response_waiter(p,MEG) % call response_waiter function
                %% response waiter function ends
                if ~p.MEG_emulator_enabled; KbQueueFlush(); end % flush the response queue from the response waiter
                if p.MEG_enabled == 1
                    WaitSecs(0.5);
                    MEG.WaitForButtonPress(0); % reset MEG button press to empty
                end
            end
            
            %% present dots
            
            % first need some info about the dots based on the trial order from 'p.stim_mat'
            dots.direction = p.stim_mat(i,4); % pulls direction in degrees from 'p.stim_mat' to give to moving_dots
            if p.stim_mat(i,5) == 1 % if this trial is supposed to be easy coherence (a '1' in 'p.stim_mat')
                dots.coherence = d.easy_coherence; % then give moving_dots value from 'd.easy_coherence'
            elseif p.stim_mat(i,5) == 2 % if the trial is supposed to be hard coherence (a '2' in 'p.stim_mat')
                dots.coherence = d.hard_coherence; % then instead give moving_dots value from 'd.hard_coherence'
            end
            
            % initiates training protocol (reduce dots presentation time from 'p.training_dots_duration' to 'p.dots_duration' by one 'p.training_reduction' every 'p.training_interval')
            if p.training_enabled == 1
                if block == 1 && i == 1
                    fprintf('training (p.training_enabled set to 1) - will reduce dots presentation from p.training_dots_duration (%u secs) to p.dots_duration (%u secs) trial by trial\n', p.training_dots_duration, p.dots_duration);
                    t.orig_dots_duration = p.dots_duration; % save experimental dots duration for later
                    p.dots_duration = p.training_dots_duration; % set dots duration to training value to start
                end
                if ~mod(i,p.training_interval) && p.training_dots_duration > t.orig_dots_duration % if we've hit a trial that is divisible by the training interval and the current value of training dots duration is greater than our final dots duration
                    p.dots_duration = p.training_dots_duration-p.training_reduction; % reduce dots duration
                    p.training_dots_duration = p.training_dots_duration-p.training_reduction; % save that value for the next trial
                    fprintf('reducing p.training_dots_duration (%u secs) by p.training_reduction (%u secs)\n', p.training_dots_duration, p.training_reduction);
                elseif ~mod(i,p.training_interval) && p.training_dots_duration <= t.orig_dots_duration % if we've hit a trial that is divisible by the training interval and the current value of training dots duration is less than or equal to our final dots duration
                    p.dots_duration = t.orig_dots_duration; % revert to the experimental dots duration
                    fprintf('dots duration now stable at p.dots_duration (%u secs)\n', p.dots_duration);
                end
            end
            
            % intertrial period - display fixation
            if p.iti_on == 1
                t.centre = p.resolution/2;
                t.sz_l = angle2pix(p,0.5/2); % this value (0.5/2) comes from p.fixation.size specified in movingdots.m
                t.iti_rect = [-t.sz_l+t.centre(1),-t.sz_l+t.centre(2),t.sz_l+t.centre(1),t.sz_l+t.centre(2)];
                t.sz_s = angle2pix(p,0.5/4); % this value (0.5/4) comes from p.fixation.size specified in movingdots.m
                t.iti_rect_sml = [-t.sz_s+t.centre(1),-t.sz_s+t.centre(2),t.sz_s+t.centre(1),t.sz_s+t.centre(2)];
                Screen('FillOval', p.win, [255,255,255],t.iti_rect);
                Screen('FillOval', p.win, [0,0,0],t.iti_rect_sml);
                Screen('Flip', p.win);
                if p.MEG_enabled == 1
                    MEG.SendTrigger(0); % reset triggers
                    WaitSecs(p.iti_time-p.MEG_onset_trigger_time);
                    MEG.SendTrigger(p.stim_mat(i,p.MEGtriggers.onsets)); % send a trigger for trial onset
                    WaitSecs(p.MEG_onset_trigger_time);
                    MEG.SendTrigger(0); % reset triggers
                elseif p.MEG_enabled == 0
                    WaitSecs(p.iti_time);
                end
                % fixation will remain in place until next flip called, else can call here %Screen('Flip', p.win);
            end
            
            % now run moving_dots
            if ~p.MEG_emulator_enabled; KbQueueFlush(); end % flush the response queue so any accidental presses recorded in the cue period won't affect responses in the dots period
            if p.MEG_enabled == 1
                MEG.WaitForButtonPress(0); % reset MEG button press to empty
            end
            %% moving dots function begins
            [d.dots_onset(block,i), t.pressed, t.firstPress] = moving_dots(p,dots,MEG,i); % pull time of first flip for dots, as well as information from KBQueueCheck from moving_dots
            
            %% deal with response and provide feedback (or abort if 'p.quitkey' pressed)
            
            % code response info
            if p.MEG_enabled == 0
                if nnz(t.firstPress) > 1 % (if number of non-zero elements is greater than 1 - i.e. if participant has responded more than once)
                    d.multiple_keypress_name(:,i,block) = KbName(t.firstPress); % record all the keypresses
                    d.multiple_keypress_time(:,i,block) = t.firstPress(t.firstPress>0)' - d.dots_onset(block,i); % record all the rts
                    t.firstPress(find(t.firstPress==max(t.firstPress))) = 0; % zero out the max value so only the first press is recorded in the matrix
                end
                d.resp_key_name{block,i} = KbName(t.firstPress); % get the name of the key used to respond - needs to be squiggly brackets or it wont work for no response
                d.resp_key_time(block,i) = sum(t.firstPress); % get the timing info of the key used to respond
                d.rt(block,i) = d.resp_key_time(block,i) - d.dots_onset(block,i); % rt is the timing of key info - time of dots onset (if you get minus values something's wrong with how we deal with nil/early responses)
            elseif p.MEG_enabled == 1
                %                 if exist('t.firstPress.multipress','var')
                %                     d.multiple_keypresses(i,:,block) = t.firstPress;
                %                 end
                d.resp_key_name(block,i) = t.firstPress{1}; % get response key from array
                d.rt(block,i) = t.firstPress{2}; % get response time from array - don't need to minus dots onset here because we're using the MEG timing functions
            end
            
            % create variable for correct response
            if p.stim_mat(i,7) == 1 % if trial matches blue
                d.correct_resp{block,i} = p.resp_keys{1};
                d.incorrect_resp{block,i} = p.resp_keys{2};
            elseif p.stim_mat(i,7) == 2 % if trial matches orange
                d.correct_resp{block,i} = p.resp_keys{2};
                d.incorrect_resp{block,i} = p.resp_keys{1};
            end
            
            % score response
            if strcmp(d.resp_key_name(block,i), d.correct_resp(block,i))
                d.correct(block,i) = 1; %correct trial
                t.feedback = 'correct';
            elseif strcmp(d.resp_key_name(block,i), d.incorrect_resp(block,i))
                d.correct(block,i) = 0; %incorrect trial
                t.feedback = 'incorrect';
            elseif strcmp(d.resp_key_name(block,i),p.quitkey)
                fclose('all');
                error('%s quit by user (p.quitkey pressed)\n', mfilename);
            else
                d.correct(block,i) = -1; % nil response
                d.rt(block,i) = 0;
                t.feedback = 'no valid input';
            end % end check correct
            
            % display some feedback if trialwise feedback on
            if p.feedback_type == 1
                DrawFormattedText(p.win, t.feedback, 'center', 'center', p.text_colour); %display feedback
                Screen('Flip', p.win);
                WaitSecs(p.feedback_time);
                Screen('Flip', p.win);
            end
            
            %% post trial clean up
            if ~p.MEG_emulator_enabled; KbQueueRelease(); end
            
            % save the trial data
            save(save_file); % save all data in '.mat' format
            
        end % trial loop
        
        % display some feedback if blockwise feedback on
        if p.feedback_type == 2
            t.block_pc = sum(d.correct(block,:)==1)/numel(d.correct(block,:));
            t.blockfeedback = round(100*t.block_pc);
            sprintf('\n %u percent correct\n',t.blockfeedback)
            DrawFormattedText(p.win, sprintf('\n you got %u percent correct\n',t.blockfeedback), 'center', 'center', p.text_colour); %display feedback
            Screen('Flip', p.win);
            WaitSecs(p.feedback_time);
            Screen('Flip', p.win);
        end
        
    end % end block loop
    
    %% wrap up
    
    
    % tell them it's over
    DrawFormattedText(p.win,'done!', 'center', 'center', p.text_colour); % tell them it's over!
    Screen('Flip', p.win);
    WaitSecs(1);
    Screen('Flip', p.win);
    
    % close screen
    ShowCursor;
    if ~p.MEG_emulator_enabled; KbQueueRelease(); end %KbReleaseWait();
    if p.MEG_enabled == 1; MEG.delete; end % stop MEG from limiting button presses
    clear block i ans; % clear specific indexes and stuff
    Screen('Close',p.win);
    
    fprintf('done running %s\n', mfilename);
    
catch err
    save(save_file);
    ShowCursor;
    if ~p.MEG_emulator_enabled; KbQueueRelease(); end %KbReleaseWait();
    %     if p.MEG_enabled == 1; MEG.delete; end % stop MEG from limiting button presses
    sca; %Screen('Close',p.win);
    rethrow(err);
end

%% subfunctions
% can put subfunctions here if running MATLAB at least R2018a
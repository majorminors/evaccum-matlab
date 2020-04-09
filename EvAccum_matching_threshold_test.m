%% Find matching threshold values for evidence accumulation task
% Dorian Minors
% Created: JUL19
% Last Edit: 13SEP19

% trial settings all saved in 'p'

% bugnotes: currently relying on 'sum' in training, rather than 'movsum' since that's a newer feature. This
%           makes things messier - have to set the minimum number of trials in the training protocol to be higher
%           than p.training_trials_per_level when you start assessing the
%           percent correct

% expects coherence thresholds to be saved in datadir, with specific naming
%   convention and variable names. see "set up participant info and save"
%   section of script

% data related information saved in 'd'
%   d.stim_mat_all contains trial condition matrices for each block
%       (d.stim_mat_all(:,:,[block number])
%   rows indicate trials
%   columns are explained when we define the matrix in 'define stimulus
%       parameters', below
%   in this script, the two blocks are two tests:
%       block 1 tests using easy coherence threshold.
%       block 2 tests using hard coherence threshold.

% other trial specific variables are in 't' in case something goes wrong
%   and we want to see them

% note: response keys for the trial based off whether participant id is odd or even

% note: block 1 tests using easy coherence threshold.
%       block 2 tests using hard coherence threshold.
%       stored as d.stim_mat_all(:,:,[block number])

% note: for backwards compatibility (to 2016), subfunctions are disabled
%       in scripts (although functions can have subfunctions and scripts
%       instead use external functions.

% note: for compatibility with testing lab computers, KbQueueWait() only accepts
%       one argument - but from at least 2018, KbQueueWait can accept two,
%       which is commented out next to any call to KbQueueWait here.

%% subfunctions (at bottom):
%% pix = angle2pix(p,ang)
% calculates pixel size from visual angles, assuming isotropic (square) pixels

% requires:
% p.screen_distance           distance from screen in cm
% p.screen_width              width of screen in cm
% p.resolution                number of pixels of p in horizontal direction - this needs to be calculated after the window is opened for the dots

%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);
p = struct(); % est structure for parameter values
d = struct(); % est structure for trial data
t = struct(); % another structure for untidy trial specific floating variables that we might want to interrogate later if we mess up

% initial settings
rootdir = pwd; % root directory - used to inform directory mappings
p.training_enabled = 0; % if 1, initiates training protocol (reduce dots presentation time from 'p.training_dots_duration' to 'p.dots_duration' by one 'p.training_reduction' every 'p.training_interval') - see '% training variables' below. will also suppress pre-block information about whether it's a hard or easy test
p.screen_width = 35;   % Screen width in cm
p.screen_height = 50;    % Screen height in cm
p.screen_distance = 50; % Screen distance from participant in cm

% general settings
p.manually_set_coherence = 0; % if 1, will include prompts to set coherence manually
p.screen_num = 0; % screen to display experiment on (0 unless multiple screens)
p.fullscreen_enabled = 1; % 1 is full screen, 0 is whatever you've set p.window_size to
p.testing_enabled = 0; % change to 0 if not testing (1 skips PTB synctests and sets number of trials and blocks to test values) - see '% test variables' below
p.fix_trial_time = 0; % if 0 then trial will end on keypress, if 1 will go for duration of p.dots_duration
p.num_blocks = 2; % each block currently feeds the two coherence values (block 1 is easy, block 2 is hard)

% for use in MEG
p.MEG_enabled = 0;
MEG = 0;

% check set up
if ~ismember(p.fullscreen_enabled,[0,1]); error('invalid value for p.fullscreen_enabled'); end % check if valid or error
if ~ismember(p.testing_enabled,[0,1]); error('invalid value for p.testing_enabled'); end % check p.testing_enabled is a valid number, or error
if ~ismember(p.training_enabled,[0,1]); error('invalid value for p.training_enabled'); end % check if valid or error
if ~ismember(p.fix_trial_time,[0,1]); error('invalid value for p.fix_trial_time'); end % check if valid or error
if p.num_blocks ~= 2; error('currently only two blocks are interpretable - check p.num_blocks'); end % check if valid or error

% directory mapping
addpath(genpath(fullfile(rootdir, 'tools_exp'))); % add tools folder to path (includes moving_dots function which is required for dot motion, as well as an external copy of subfunctions for backwards compatibility with MATLAB)
stimdir = fullfile(rootdir, 'tools_exp', 'stimuli');
datadir = fullfile(rootdir, 'data'); % will make a data directory if none exists
if ~exist(datadir,'dir')
    mkdir(datadir);
end

% training variables
p.training_dots_duration = 5; % duration (secs) dot stimulus presented for for training trials (will reduce trial by trial to p.dots_duration)
p.training_reduction = 1; % by how much do we reduce the duration during training (seconds)
p.training_interval = 2; % how many trials should we train on before reducing the dots presentation time
p.training_trials_per_level = 8; % how many trials to run each training level at (currently 3 levels - levels 1 and 2 will wait for p.training_percent_correct/p.training_trials_per_level while other levels will just do three and move on)
p.training_percent_correct = 0.75; % percent of correct trials requred to move on
p.training_rule_lvl_easy = 10; % what matching rule level for easy (training level 1)
p.training_rule_lvl_mid = 45; % what matching rule level for mid (training level 2)
p.training_rule_lvl_hard = 75; % what matching rule level for hard (training level 3)

% test variables
p.num_test_trials = 4;
p.num_test_blocks = 2;
if p.testing_enabled == 1
    p.PTBsynctests = 1; % PTB will skip synctests if 1
    p.PTBverbosity = 1; % PTB will only display critical warnings with 1
elseif p.testing_enabled == 0
    p.PTBsynctests = 0;
    p.PTBverbosity = 3; % default verbosity for PTB
end
Screen('Preference', 'SkipSyncTests', p.PTBsynctests);
Screen('Preference', 'Verbosity', p.PTBverbosity);
if p.testing_enabled > 1; error('invalid value for p.testing_enabled'); end % check p.testing_enabled is a valid number, or error

% psychtoolbox setup
AssertOpenGL; % check Psychtoolbox (on OpenGL) and Screen() is working
KbName('UnifyKeyNames'); % makes key mappings compatible (mac/win)
rng('shuffle'); % seed rng using date and time

% set up participant info and save
if p.manually_set_coherence
    t.prompt = {'enter participant number:','enter easy coherence threshold (fm 0-1, higher is easier)','enter hard coherence threshold (fm 0-1, lower is harder)'}; % prompt a dialog to enter subject info
    t.prompt_defaultans = {num2str(99), num2str(0.75), num2str(0.25)}; % default answers corresponding to prompts
    t.prompt_rsp = inputdlg(t.prompt, 'enter participant info', 1, t.prompt_defaultans); % save dialog responses
    d.participant_id = str2double(t.prompt_rsp{1}); % add subject number to 'd'
    d.easy_coherence = str2double(t.prompt_rsp{2}); % add participant coherence thresholds to 'd'
    d.hard_coherence = str2double(t.prompt_rsp{3}); % add participant coherence thresholds to 'd'
elseif ~p.manually_set_coherence
    t.prompt = {'enter participant number:'}; % prompt a dialog to enter subject info
    t.prompt_defaultans = {num2str(99)}; % default answers corresponding to prompts
    t.prompt_rsp = inputdlg(t.prompt, 'enter participant info', 1, t.prompt_defaultans); % save dialog responses
    d.participant_id = str2double(t.prompt_rsp{1}); % add subject number to 'd'
    tmp = load(fullfile(datadir,[num2str(d.participant_id,'S%02d'),'_EvAccum_coherence_threshold_test.mat']),'d');
    d.easy_coherence = tmp.d.easy_threshold;
    d.hard_coherence = tmp.d.hard_threshold;
    clear tmp;
end

% check participant info has been entered correctly for the script
if isnan(d.participant_id)
    error('no participant number entered');
elseif isnan(d.easy_coherence) || d.easy_coherence > 1
    error('invalid participant coherence threshold (easy)');
elseif isnan(d.hard_coherence) || d.hard_coherence > 1
    error('invalid participant coherence threshold (hard)')
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
p.quitkey = {'q'};
p.resp_keys = {'a','s'}; % only accepts two response options
% establish response keys for the trial based off whether participant id is odd or even
if mod(d.participant_id,2) % if pid not divisible by 2 (i.e. leaves a modulus after division)
    p.resp_keys = {p.resp_keys{1},p.resp_keys{2}}; % essentially do nothing - keep resp keys the same
elseif ~mod(d.participant_id,2) % if pid is divisible by 2 (i.e. does not leave a modulus after division)
    p.resp_keys = {p.resp_keys{2},p.resp_keys{1}}; % then swap the response keys
end

% define display info
p.bg_colour = [0 0 0]; % needs to be the same as the cue stimuli background colour (unless transparent)
p.text_colour= [255 255 255]; % colour of instructional text
p.cue_colour_one = [255 255 255]; %[121 181 240]; % colour of cue, for text formatting
p.cue_colour_two = [255 255 255]; %[240 181 121]; % colour of cue, for text formatting
p.text_size = 40; % size of text
p.window_size = [0 0 1200 800]; % size of window when ~p.fullscreen_enabled
p.visual_angle_cue = 10; % visual angle of the cue expressed as a decimal - determines size
p.visual_angle_dots = 0.15; % visual angle of the dots expressed as a decimal - determines size

% timing info
p.min_cue_time = 1; % minimum period to display cue (participants can't continue during this time)
p.dots_duration = 2; % seconds for the dot cloud to be displayed
p.min_resp_mapping_time = 1; % minimum period to display response mapping (participants can't continue during this time)
p.feedback_time = 0.5; % period to display feedback after response

% trial settings (*p.stim_mat* = parameter required to calculate stimulus condition matrix)
p.num_trials_per_block = 160; % *p.stim_mat* - must be divisible by p.num_cues && >= 15*p.num_points
p.num_cues = 4; % *p.stim_mat*
p.cue_directions = 45:90:315; % *p.stim_mat* - refers to the direction of the upward arrow of the doublesided arrow cue in stimdir
p.num_points = 10; % *p.stim_mat* - number of points to test participants on for each test
p.rule_points = 0:10:90; %union([0:5:20],[70:5:90]); % *p.stim_mat* - length(p.rule_points) must == p.num_points
if length(p.rule_points) ~= p.num_points; error('number of coherence points is not equal to the number of testing points you specified'); end % check length(p.rule_points) == p.num_points, or error

%% define stimuli parameters

fprintf('defining stimuli params for %s\n', mfilename);

% read in stimulus file for the cue
p.cue = imread(fullfile(stimdir, 'arrows_cue_with_line.png'));

% create matrix specifying stimulus conditions per trial:
%    1)  cue direction (1-4) - evenly allocates trials to cues
%    2)  blue cue direction in degrees for each trial - evenly adds cue
%        directions to trials in a similar manner to (1)
%    3)  coherence condition (1 = match blue, 2 = match orange) - evenly allocates
%        trials
%    4)  point condition (1-10) - each repeated
%        p.num_trials_per_block/p.num_points times
%    5)  matching points allocated to point conditions
%    6)  whether trial should add or subtract degrees from cues test (1 =
%        add, 2 = subtract) - currently 2 of each per cue (since four reps
%        of point conditions per cue)
%    7)  coherence direction from cued direction in degrees - calculated from cue direction and matching
%        point
%    8)  coherence direction in degrees incorporating blue match or orange match

% create matrix
p.stim_mat = zeros(p.num_trials_per_block,8);
% create columns
p.stim_mat(:,1) = sort(repmat(1:p.num_cues,[1,p.num_trials_per_block/p.num_cues]));
p.stim_mat(:,2) = p.cue_directions(p.stim_mat(:,1));
p.stim_mat(:,3) = repmat(sort(repmat(1:2,[1,p.num_trials_per_block/p.num_cues/2])),[1,p.num_cues]);
p.stim_mat(:,4) = repmat(1:p.num_points,[1,p.num_trials_per_block/p.num_points]);
p.stim_mat(:,5) = repmat(p.rule_points,[1,p.num_trials_per_block/p.num_points]);
p.stim_mat(:,6) = repmat(sort(repmat(1:2,[1,p.num_trials_per_block/p.num_cues/4])),[1,p.num_cues*2]);
add = (p.stim_mat(:,6)==1); % index all posns in (6) with 1s into 'add'
subtract = (p.stim_mat(:,6)==2); % index all posns in (6) with 2s into 'subtract'
p.stim_mat(add,7) = p.stim_mat(add,2)+p.stim_mat(add,5); % insert addition of matching point to cue direction into (7)
p.stim_mat(subtract,7) = p.stim_mat(subtract,2)-p.stim_mat(subtract,5); % insert subtraction of matching point to cue direction into (7)
temp = p.stim_mat(:,7)>360; % get all places where (7) > 360 degrees
p.stim_mat(temp,7) = p.stim_mat(temp,7)-360; % make any in (7) over 360 wrap around from 0 again
temp = p.stim_mat(:,7)<0; % get all places where (7) < 0 degrees
p.stim_mat(temp,7) = p.stim_mat(temp,7)+360; % make any in (7) under 0 wrap around from 360 again
temp = p.stim_mat(:,7)+180; % get 180 degrees on each cue
temp(temp>360) = temp(temp>360)-360; % make any over 360 wrap around from 0 again
p.stim_mat(:,8) = times(p.stim_mat(:,3)-1,temp); % make 1s and 2s from (3) into 0s and 1s, then times by temp (so becomes 0s and cue+180degs)
temp = (p.stim_mat(:,8)==0); % index all posns in (8) with 0s into 'temp'
p.stim_mat(temp,8) = p.stim_mat(temp,7); % insert coherence directions from (7) into (8) where there are 0s (using 'temp' index)
% clear floating variables
clear temp add subtract;

% shuffle the trial condition order for each block into 'd.stim_mat_all' then re-sort by cue
for block=1:p.num_blocks
    d.stim_mat_all(:,:,block) = p.stim_mat(Shuffle(1:p.num_trials_per_block),:);
    d.stim_mat_all(:,:,block) = sortrows(d.stim_mat_all(:,:,block),1);
end
% clear floating variables
clear block;

% dot cloud parameters for moving_dots function
%   note that for the following required params for moving_dots:
%       'dots.direction' is specified within the trial loop
%       'dots.coherence' is also specified withing the trial loop using outputs
%           from 'p.stim_mat' and 'd.easy_coherence'/'d.hard_coherence'
dots.aperture_size = [10,10]; % [width,height] in degrees of the rectangular aperture dots are displayed in
dots.centre = [0,0]; % [x,y] centre of the dot cloud
dots.num_dots = 100; % number of dots
dots.colour = [255,255,255]; % colour of the dots in [r,g,b]
dots.visual_angle = p.visual_angle_dots; % visual angle of the dots expressed as a decimal - determines size
dots.speed = 5; % speed of dots in degrees per second
dots.lifetime = 5; % number of frames dots live for

%% test start

fprintf('running matching threshold assessment (%s)\n', mfilename);

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
        
        % tell participant whether hard or easy test if training not enabled
        if block == 1 && p.training_enabled == 0
            DrawFormattedText(p.win,'first with easy dots', 'center', 'center', p.text_colour); %display feedback
            Screen('Flip', p.win);
            WaitSecs(1);
            Screen('Flip', p.win);
        elseif block == 2 && p.training_enabled == 0
            DrawFormattedText(p.win,'now with harder dots', 'center', 'center', p.text_colour); %display feedback
            Screen('Flip', p.win);
            WaitSecs(1);
            Screen('Flip', p.win);
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
        i = 0; % initialise trial index
        while i < p.act_trial_num
            i = i + 1;
            fprintf('trial %u of %u\n',i,p.act_trial_num); %report trial number to command window
            
            %set up a queue to collect response info
            t.queuekeys = [KbName(p.resp_keys{1}), KbName(p.resp_keys{2}), KbName(p.quitkey)]; % define the keys the queue cares about
            t.queuekeylist = zeros(1,256); % create a list of all possible keys (all 'turned off' i.e. zeroes)
            t.queuekeylist(t.queuekeys) = 1; % 'turn on' the keys we care about in the list (make them ones)
            KbQueueCreate([], t.queuekeylist); % initialises queue to collect response information from the list we made (not listening for response yet)
            KbQueueStart(); % starts delivering keypress info to the queue
            
            %% present cue and response mapping
            if i == 1 || p.stim_mat(i,1) ~= p.stim_mat(i-1,1) || ~mod(i-1,8) %|| i > p.act_trial_num-1 && ~mod(i,8) && p.stim_mat(i+1,1) == p.stim_mat(i,1) % if first trial, or cue changes (as currently blocked), or every 8 trials unless we're about to change cue then display cue
                
                % save the trial data
                save(save_file); % save all data in '.mat' format
                
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
                if p.stim_mat(i,1) == 1
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{1}), t.rect(RectRight)*1.01, t.rect(RectTop)*1.01, p.cue_colour_one);
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{2}), t.rect(RectLeft)*1.01, t.rect(RectBottom)*1.01, p.cue_colour_two);
                elseif p.stim_mat(i,1) == 2
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{1}), t.rect(RectRight)*1.01, t.rect(RectBottom)*1.01, p.cue_colour_one);
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{2}), t.rect(RectLeft)*1.01, t.rect(RectTop)*1.01, p.cue_colour_two);
                elseif p.stim_mat(i,1) == 3
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{1}), t.rect(RectLeft)*1.01, t.rect(RectBottom)*1.01, p.cue_colour_one);
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{2}), t.rect(RectRight)*1.01, t.rect(RectTop)*1.01, p.cue_colour_two);
                elseif p.stim_mat(i,1) == 4
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{1}), t.rect(RectLeft)*1.01, t.rect(RectTop)*1.01, p.cue_colour_one);
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{2}), t.rect(RectRight)*1.01, t.rect(RectBottom)*1.01, p.cue_colour_two);
                end
                d.cue_onset(block,i) = Screen('Flip', p.win); % pull the time of the screen flip from the flip function while flipping
                WaitSecs(p.min_cue_time);
                Screen('DrawTexture', p.win, p.cue_tex, [], t.rect, p.stim_mat(i,2)); % redraws the cue
                DrawFormattedText(p.win, sprintf('\n press either button to continue\n'), 'center', t.rect(RectBottom)*1.1, p.text_colour);
                if p.stim_mat(i,1) == 1
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{1}), t.rect(RectRight)*1.01, t.rect(RectTop)*1.01, p.cue_colour_one);
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{2}), t.rect(RectLeft)*1.01, t.rect(RectBottom)*1.01, p.cue_colour_two);
                elseif p.stim_mat(i,1) == 2
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{1}), t.rect(RectRight)*1.01, t.rect(RectBottom)*1.01, p.cue_colour_one);
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{2}), t.rect(RectLeft)*1.01, t.rect(RectTop)*1.01, p.cue_colour_two);
                elseif p.stim_mat(i,1) == 3
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{1}), t.rect(RectLeft)*1.01, t.rect(RectBottom)*1.01, p.cue_colour_one);
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{2}), t.rect(RectRight)*1.01, t.rect(RectTop)*1.01, p.cue_colour_two);
                elseif p.stim_mat(i,1) == 4
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{1}), t.rect(RectLeft)*1.01, t.rect(RectTop)*1.01, p.cue_colour_one);
                    DrawFormattedText(p.win, sprintf('\n %s \n', p.resp_keys{2}), t.rect(RectRight)*1.01, t.rect(RectBottom)*1.01, p.cue_colour_two);
                end
                Screen('Flip', p.win);
                t.waiting = []; % wait for user input
                while isempty(t.waiting)
                    t.waiting = KbQueueWait();%([],3);
                end
                
                KbQueueFlush(); % flush the response queue from the continue screen presses
            end
            
            %% present dots and collect response info
            
            % first need some info about the dots based on the trial order from 'p.stim_mat'
            dots.direction = p.stim_mat(i,8); % pulls direction in degrees from 'p.stim_mat' to give to moving_dots
            if block == 1
                dots.coherence = d.easy_coherence; % then give moving_dots value from 'd.easy_coherence'
            elseif block == 2
                dots.coherence = d.hard_coherence; % then give moving_dots value from 'd.hard_coherence'
            end
            
            % training protocol
            % 1) will reduce dots presentation time from
            %    'p.training_dots_duration' to 'p.dots_duration' by one
            %    'p.training_reduction' every 'p.training_interval' to
            %    familiarise participants with timing
            % 2) will then start training difficulty levels - first will
            %    wait for 'p.training_percent_correct' out of the last p.training_trials_per_level before
            %    moving on - thereafter each level will just run the number
            %    of training trials regardless of correct/incorrect
            if p.training_enabled == 1 % if initiated
                if block == 1 && i == 1
                    fprintf('training (p.training_enabled set to 1)\n will reduce dots presentation from p.training_dots_duration (%u secs) to p.dots_duration (%u secs) every %u trials\n', p.training_dots_duration, p.dots_duration, p.training_interval);
                    t.orig_dots_duration = p.dots_duration; % save experimental dots duration for later
                    p.dots_duration = p.training_dots_duration; % set dots duration to training value to start
                    t.training_placeholder = p.training_rule_lvl_easy;
                    t.training_level = 0; % counter variables
                    %d.correct = 0; % have to set this to a value for the first trial or the code will break (I think - I didn't
                end
                if ~mod(i,p.training_interval) && p.training_dots_duration > t.orig_dots_duration % if we've hit a trial that is divisible by the training interval and the current value of training dots duration is greater than our final dots duration
                    p.dots_duration = p.training_dots_duration-p.training_reduction; % reduce dots duration
                    p.training_dots_duration = p.training_dots_duration-p.training_reduction; % save that value for the next trial
                    t.training_placeholder = p.training_rule_lvl_easy;
                    fprintf('reducing p.training_dots_duration (%u secs) by p.training_reduction (%u secs)\n', p.training_dots_duration, p.training_reduction);
                elseif ~mod(i,p.training_interval) && p.training_dots_duration <= t.orig_dots_duration % if we've hit a trial that is divisible by the training interval and the current value of training dots duration is less than or equal to our final dots duration
                    p.dots_duration = t.orig_dots_duration; % revert to the experimental dots duration
                    fprintf('dots duration now stable at p.dots_duration (%u secs) - now exposing to angle changes\n', p.dots_duration);
                    if t.training_level == 0 && i >= 9
                        fprintf('training level %u at easy angle - require %u (percent) correct to continue\n', t.training_level, p.training_percent_correct);
                        t.training_placeholder = p.training_rule_lvl_easy;
                        t.training_correct_counter = sum(d.correct(block,i-p.training_trials_per_level:i-1)); % movsum(d.correct,[p.training_trials_per_level 0]);
                        if t.training_correct_counter >= p.training_trials_per_level*p.training_percent_correct
                            t.training_level = t.training_level+1;
                            t.trialcounter = i;
                        end
                    end
                    if t.training_level == 1
                        fprintf('training level %u at mid angle - require %u (percent) correct to continue\n', t.training_level, p.training_percent_correct); % fprintf('training level %u at mid angle - exposure for %u trials\n', t.training_level, p.training_trials_per_level);
                        t.training_placeholder = p.training_rule_lvl_mid;
                        t.training_correct_counter = sum(d.correct(block,i-p.training_trials_per_level:i-1)); % movsum(d.correct,[p.training_trials_per_level 0]);
                        if i == t.trialcounter+p.training_trials_per_level && t.training_correct_counter >= p.training_trials_per_level*p.training_percent_correct
                            t.training_level = t.training_level+1;
                            t.trialcounter = i;
                        end
                    end
                    if t.training_level == 2
                        fprintf('training level %u at hard angle - exposure for %u trials\n', t.training_level, p.training_trials_per_level);
                        t.training_placeholder = p.training_rule_lvl_hard;
                        if i == t.trialcounter+p.training_trials_per_level
                            t.training_level = t.training_level+1;
                        end
                    end
                    if t.training_level == 3
                        fclose('all');
                        error('training complete for %s\n', mfilename);
                    end
                end
                if i == 1 || mod(i,2)
                    t.trainer = p.stim_mat(i,2)+t.training_placeholder;
                elseif ~mod(i,2)
                    t.trainer = p.stim_mat(i,2)-t.training_placeholder;
                end
                if p.stim_mat(i,3) == 2
                    t.trainer = t.trainer+180;
                    if t.trainer > 360
                        t.trainer = t.trainer-360;
                    end
                end
                dots.direction = t.trainer;
            end
            
            % now run moving_dots
            KbQueueFlush(); % flush the response queue so any accidental presses recorded in the cue period won't affect responses in the dots period
            [d.dots_onset(block,i), t.pressed, t.firstPress] = moving_dots(p,dots); % pull time of first flip for dots, as well as information from KBQueueCheck from moving_dots
            
            %% deal with response and provide feedback (or abort if 'p.quitkey' pressed)
            
            % code response info
            d.resp_key_name{block,i} = KbName(t.firstPress); % get the name of the key used to respond - needs to be squiggly brackets or it wont work for no response
            d.resp_key_time(block,i) = sum(t.firstPress); % get the timing info of the key used to respond
            d.rt(block,i) = d.resp_key_time(block,i) - d.dots_onset(block,i); % rt is the timing of key info - time of dots onset (if you get minus values something's wrong with how we deal with nil/early responses)
            
            % create variable for correct response
            if p.stim_mat(i,3) == 1 % if trial matches blue
                d.correct_resp(block,i) = p.resp_keys{1};
                d.incorrect_resp(block,i) = p.resp_keys{2};
            elseif p.stim_mat(i,3) == 2 % if trial matches orange
                d.correct_resp(block,i) = p.resp_keys{2};
                d.incorrect_resp(block,i) = p.resp_keys{1};
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
            
            % display feedback
            DrawFormattedText(p.win,t.feedback, 'center', 'center', p.text_colour); %display feedback
            Screen('Flip', p.win);
            WaitSecs(p.feedback_time);
            Screen('Flip', p.win);
            
            %% post trial clean up
            
            % clear the response queue for the next trial and related floating variables
            KbQueueRelease();
            
        end % trial loop
        
    end % end block loop
    
    %% wrap up
    
    % save the trial data
    save(save_file); % save all data in '.mat' format
    
    % tell them it's over
    DrawFormattedText(p.win,'done!', 'center', 'center', p.text_colour); %display feedback
    Screen('Flip', p.win);
    WaitSecs(1);
    Screen('Flip', p.win);
    
    % close screen
    ShowCursor;
    KbQueueRelease(); %KbReleaseWait();
    clear block i ans; % clear specific indexes and stuff
    Screen('Close',p.win);
    
    %% run analysis
    if ~p.training_enabled
        [d.easy_threshold,d.hard_threshold,d.test_overview,d.test_summary] = match_thresholding(p,d,save_file);
    end
    
    % save the analysis results
    save(save_file); % save all data in '.mat' format
    
    %% finish
    
    fprintf('done running %s\n', mfilename);
    
catch err
    save(save_file);
    ShowCursor;
    KbQueueRelease(); %KbReleaseWait();
    sca; %Screen('Close',p.win);
    rethrow(err);
end

%% subfunctions

%% convert visual angles in degrees to pixels
% function pix = angle2pix(p,ang)
% pixSize = p.screen_width/p.resolution(1);   % cm/pix
% sz = 2*p.screen_distance*tan(pi*ang/(2*180));  %cm
% pix = round(sz/pixSize);   % pix
% return
% end
function f = keyswap_training(p,dots,d)
% keyswap_training(p,dots,d)
%
% will run some training on the EvAccum paradigm in a sandbox when called.
%
% requires KBQueueStart() called prior to function (so may want to
% call KBQueueFlush() prior to function).
% 
% will create KBQueue again before it wraps up, but be sure
% that queue is the same as the one you created outside the function.
%
% requires:
%   structure 'p' with experiment parameters from EvAccum
%       note: this utilises 'p.training' parameters
%   structure 'dots' with dot cloud parameters for moving_dots function
%   except the following which are specified in the function:
%       dots.direction
%       dots.coherence
%   from structure 'd':
%       d.stim_mat_all
%       d.easy_coherence
%       d.hard_coherence
%
% returns:
%   all output as a structure 'f'
%

%% NOV 19 - function by Dorian Minors
%% start function

p.stim_mat = d.stim_mat_all(:,:,1); % take the stimulus matrix from the first block to train on
f.max_trials = length(p.stim_mat(:,1));
block = 1; % couldn't be bothered removing all the references to block in the code, so I'll just do this

% open trial loop
i = 0;
while i < f.max_trials
    i = i + 1;
    fprintf('training trial %u of %u\n',i,f.max_trials); % report trial number to command window
    
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
        f.cue_onset(block,i) = Screen('Flip', p.win); % pull the time of the screen flip from the flip function while flipping
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
        t.waiting = []; % wait for user input
        while isempty(t.waiting)
            t.waiting = KbQueueWait();%([],3);
        end
        
        KbQueueFlush(); % flush the response queue from the continue screen presses
    end
    
    %% present dots
    
    % first need some info about the dots based on the trial order from 'p.stim_mat'
    dots.direction = p.stim_mat(i,4); % pulls direction in degrees from 'p.stim_mat' to give to moving_dots
    if p.stim_mat(i,5) == 1 % if this trial is supposed to be easy coherence (a '1' in 'p.stim_mat')
        dots.coherence = d.easy_coherence; % then give moving_dots value from 'd.easy_coherence'
    elseif p.stim_mat(i,5) == 2 % if the trial is supposed to be hard coherence (a '2' in 'p.stim_mat')
        dots.coherence = d.hard_coherence; % then instead give moving_dots value from 'd.hard_coherence'
    end
    
    % training protocol (reduce dots presentation time from 'p.training_dots_duration' to 'p.dots_duration' by one 'p.training_reduction' every 'p.training_interval')
    if block == 1 && i == 1
        fprintf('training - will reduce dots presentation from p.training_dots_duration (%u secs) to p.dots_duration (%u secs) trial by trial\n', p.training_dots_duration, p.dots_duration);
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
    
    % now run moving_dots
    KbQueueFlush(); % flush the response queue so any accidental presses recorded in the cue period won't affect responses in the dots period
    [f.dots_onset(block,i), t.pressed, t.firstPress] = moving_dots(p,dots); % pull time of first flip for dots, as well as information from KBQueueCheck from moving_dots
    
    %% deal with response and provide feedback (or abort if 'p.quitkey' pressed)
    
    % code response info
    f.resp_key_name{block,i} = KbName(t.firstPress); % get the name of the key used to respond - needs to be squiggly brackets or it wont work for no response
    f.resp_key_time(block,i) = sum(t.firstPress); % get the timing info of the key used to respond
    f.rt(block,i) = f.resp_key_time(block,i) - f.dots_onset(block,i); % rt is the timing of key info - time of dots onset (if you get minus values something's wrong with how we deal with nil/early responses)
    
    % create variable for correct response
    if p.stim_mat(i,7) == 1 % if trial matches blue
        f.correct_resp(block,i) = p.resp_keys{1};
        f.incorrect_resp(block,i) = p.resp_keys{2};
    elseif p.stim_mat(i,7) == 2 % if trial matches orange
        f.correct_resp(block,i) = p.resp_keys{2};
        f.incorrect_resp(block,i) = p.resp_keys{1};
    end
    
    % score response
    if strcmp(f.resp_key_name(block,i), f.correct_resp(block,i))
        f.correct(block,i) = 1; %correct trial
        t.feedback = 'correct';
    elseif strcmp(f.resp_key_name(block,i), f.incorrect_resp(block,i))
        f.correct(block,i) = 0; %incorrect trial
        t.feedback = 'incorrect';
    elseif strcmp(f.resp_key_name(block,i),p.quitkey)
        fclose('all');
        error('%s quit by user (p.quitkey pressed)\n', mfilename);
    else
        f.correct(block,i) = -1; % nil response
        f.rt(block,i) = 0;
        t.feedback = 'no valid input';
    end % end check correct
    
    % display some feedback
    DrawFormattedText(p.win,t.feedback, 'center', 'center', p.text_colour); %display feedback
    Screen('Flip', p.win);
    WaitSecs(p.feedback_time);
    Screen('Flip', p.win);
    
    %% post trial clean up
    KbQueueRelease();
    
    %set up a queue to collect response info
    t.queuekeys = [KbName(p.resp_keys{1}), KbName(p.resp_keys{2}), KbName(p.quitkey)]; % define the keys the queue cares about
    t.queuekeylist = zeros(1,256); % create a list of all possible keys (all 'turned off' i.e. zeroes)
    t.queuekeylist(t.queuekeys) = 1; % 'turn on' the keys we care about in the list (make them ones)
    KbQueueCreate([], t.queuekeylist); % initialises queue to collect response information from the list we made (not listening for response yet)
    KbQueueStart(); % starts delivering keypress info to the queue
    
end

% tell them it's over
DrawFormattedText(p.win,'\n all done training!\n ', 'center', 'center', p.text_colour);
Screen('Flip', p.win);
WaitSecs(1);
DrawFormattedText(p.win,'\n all done training!\n\n\n press either button to continue\n', 'center', 'center', p.text_colour);
Screen('Flip', p.win);
t.waiting = []; % wait for user input
while isempty(t.waiting)
    t.waiting = KbQueueWait();%([],3);
end
KbQueueFlush(); % flush the response queue from the continue screen presses

return
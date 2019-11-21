function [pressed, firstPress] = megtest_function(p,MEG)
% [dots_onset_time, pressed, firstPress] = megtest_function(p,dots,MEG,exp_trial)

if p.MEG_enabled == 1
    MEG.SendTrigger(1); % send a trigger for trial onset
    MEG.ResetClock; % reset the timer
    button_pressed = 0; % a counter to make sure we catch the first time a button was pressed
    MEG.WaitForButtonPress(1.5); % listen for button press
    pause(0.005); % quick pause before we reset triggers
    MEG.SendTrigger(0); % reset triggers
end

for n = 1:10
    
    while ~strcmp(MEG.LastButtonPress,p.resp_key{1}) || ~strcmp(MEG.LastButtonPress,p.resp_key{2})% isempty(MEG.LastButtonPress)
        pressed = KbQueueCheck();
        if pressed == 1; break; end % assumes you are listening for a quit key on experimenter keyboard
    end
    
    if ~isempty(MEG.LastButtonPress) && ~button_pressed % check for a keypress in the MEG key wait function every frame, if a key hasn't been pressed yet
        button_pressed = 1; % record that a key has been pressed this trial
        MEG.SendTrigger(2); % send a trigger
        firstPress{1} = MEG.LastButtonPress; % record the key pressed
        firstPress{2} = MEG.TimeOfLastButtonPress; % record the time of the key pressed
        pressed = NaN; % just put something in here
        pause(0.005); % quick pause before resetting
        MEG.SendTrigger(0); % reset the triggers
        if p.fix_trial_time == 0
            break % break the frame-loop which will end the function
        end
        fpidx = 0;
    elseif ~isempty(MEG.LastButtonPress) && button_pressed > 0 % if it's not the first time a key has been pressed
        button_pressed = 2; % record that another key has been pressed this trial
        fpidx = fpidx+1;
        firstPress{3+fpidx} = MEG.LastButtonPress; % record the key pressed
        firstPress{4+fpidx} = MEG.TimeOfLastButtonPress; % record the time of the key pressed
    end
    
end


return
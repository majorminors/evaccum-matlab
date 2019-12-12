function [pressed, firstPress] = megtest_function(p,MEG)
% [dots_onset_time, pressed, firstPress] = megtest_function(p,dots,MEG,exp_trial)

if p.MEG_enabled == 1
    MEG.SendTrigger(1); % send a trigger for trial onset
    button_pressed = 0; % a counter to make sure we catch the first time a button was pressed
    pause(0.005); % quick pause before we reset triggers
    MEG.ResetClock; % reset the timer
    MEG.SendTrigger(0); % reset triggers
end

fprintf('entering function response loop\n');

for n = 1:10
n
    

    
    fprintf('button coding block begins\n');
    if ~strcmp(MEG.LastButtonPress,p.continue_key) && ~button_pressed % check for a keypress in the MEG key wait function every frame, if a key hasn't been pressed yet
        button_pressed = 1; % record that a key has been pressed this trial
        MEG.SendTrigger(2); % send a trigger
        firstPress{1} = MEG.LastButtonPress; % record the key pressed
        firstPress{2} = MEG.TimeOfLastButtonPress; % record the time of the key pressed
        fprintf('recorded response as\n'); firstPress{1}
        pressed = NaN; % just put something in here
        pause(0.005); % quick pause before resetting
        MEG.SendTrigger(0); % reset the triggers
        fpidx = 0;
    elseif ~strcmp(MEG.LastButtonPress,firstPress{1}) && button_pressed > 0 % if it's not the first time a key has been pressed
        button_pressed = 2; % record that another key has been pressed this trial
        fpidx = fpidx+1;
        firstPress{3+fpidx} = MEG.LastButtonPress; % record the key pressed
        firstPress{4+fpidx} = MEG.TimeOfLastButtonPress; % record the time of the key pressed
        fprintf('recorded a subsequent response as\n'); firstPress{3+fpidx}
    end
    
end


return
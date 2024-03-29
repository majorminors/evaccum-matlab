function response_waiter(p,MEG)
% response_waiter(p,MEG)
% updated 18NOV19 - check function at bottom of script is same version
%
% will wait for button responses before continuing
%
% assumes you have KBQueueStart() called and listening for keys prior to function. may want to
% KbQueueFlush() following this function.
%
% requires from p:
%   p.MEG_enabled (0 or 1)
%   p.MEG_emulator_enabled (0 or 1)
% requires from MEG - a class of functions for the MEG interface with the National Instruments PCI 6503 card (MRC CBU)
%   MEG.LastButtonPress
%   MEG.WaitForButtonPress
%
% can use p.no_participant_response (1,0) to control whether participant
%   can respond, or if response limited to keyboard

%% last edit D. Minors 18 November 2019
%% start function

if p.MEG_enabled == 0 || p.no_participant_response == 1
    waiting = []; % wait for user input
    while isempty(waiting)
        waiting = KbQueueWait();%([],3) - this commented code is for later versions of MATLAB;
    end    
elseif p.MEG_enabled == 1
    MEG.WaitForButtonPress(0); % reset MEG button press to empty
    while isempty(MEG.LastButtonPress)
        MEG.WaitForButtonPress; % wait for user input
        if ~p.MEG_emulator_enabled
            [~, firstPress] = KbQueueCheck(); % listen on experimenter keyboard
            if strcmp(KbName(firstPress),p.quitkey)
                fclose('all');
                error('%s quit by user (p.quitkey pressed)\n', mfilename);
            elseif strcmp(KbName(firstPress),p.continuekey)
                break
            end 
        end
    end
end

return
%% MEG trigger test
% 
% sends triggers!

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);
p = struct(); % est structure for parameter values

% settings
p.pre_reset_pause = 0.005; % in secs
p.post_reset_pause = 1; % in secs
p.trigger_values = [1:1:255]; %[1,2,4,8];

% make sure pause is on
pause('on');

% enable MEG
MEG = MEGSynchClass;
pause(1); % I don't know if this is necessary, but seems sensible to wait a sec for it to load the class
MEG.SendTrigger(0); % reset triggers

% begin test
fprintf('beginning test: %s\n', mfilename);
for i = 1:length(p.trigger_values)
    MEG.SendTrigger(p.trigger_values(i));
	pause(p.pre_reset_pause)
	MEG.SendTrigger(0); % reset triggers
    pause(p.post_reset_pause);
end


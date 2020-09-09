%% MEG trigger test
% 
% sends triggers!

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);
p = struct(); % est structure for parameter values

% directory mapping
% addpath(); % the MEGSynchClass is in this folder, so you don't need this, but if you move things around, add the directory the class is here

% settings
p.pausetime = 1; % in seconds
p.reset_location = 1; % 1 if reset triggers before pause, 2 if reset triggers after pause
p.trigger_values = [1,2,4,8]; % [1:1:255];

% make sure pause is on
pause('on');

% enable MEG
MEG = MEGSynchClass;
pause(p.pausetime); % I don't know if this is necessary, but seems sensible to wait a sec for it to load the class
MEG.SendTrigger(0); % reset triggers

% begin test
fprintf('beginning test: %s\n', mfilename);
for i = 1:length(p.trigger_values)
    MEG.SendTrigger(p.trigger_values(i));
    pauser(p);
end



% function to order pause and reset
function pauser(p)
    if p.reset_location == 1
        MEG.SendTrigger(0); % reset triggers
        pause(p.pausetime);
    elseif p.reset_location == 2
        pause(p.pausetime);
        MEG.SendTrigger(0); % reset triggers
    end
end
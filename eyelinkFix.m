%% if the eyetracker gets stuck, run this twice (will error the first time)

clc;
clear all;
close all;
useOfEyelink=1;
dummymode=0;
edfFile='test';
trigger=1;
% connection with eyetracker, opening file
if ~EyelinkInit(dummymode)
    fprintf('Eyelink Init aborted.\n');
    cleanup(useOfEyelink, edfFile, trigger);
    return;
end
i = Eyelink('Openfile', edfFile);
if i~=0
    fprintf('Cannot create EDF file ''%s'' ', edfFile);
    cleanup(useOfEyelink, edfFile, trigger);
    return;
end
if Eyelink('IsConnected')~=1 && ~dummymode
    cleanup(useOfEyelink, edfFile, trigger);
    return;
end


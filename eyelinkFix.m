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

function cleanup(useOfEyelink, edfFile)

if useOfEyelink==1
    Eyelink('Command' , 'set_idle_mode');
    WaitSecs(0.5);
    Eyelink('CloseFile');
    %download data file
    try
        fprintf('Receiving data file ''%s''\n', edfFile );
        status=Eyelink('ReceiveFile');
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
        end
        if 2==exist(edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
        end
    catch
        fprintf('Problem receiving data file ''%s''\n', edfFile );
    end
    Eyelink('Shutdown');
end

Screen('CloseAll');

ListenChar(0);
ShowCursor;
end

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

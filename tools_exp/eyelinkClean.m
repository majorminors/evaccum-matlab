function eyelinkClean(useOfEyelink, edfFile, datadir)

if useOfEyelink==1
    Eyelink('Command' , 'set_idle_mode');
    WaitSecs(0.5);
    Eyelink('CloseFile');
    %download data file
    try
        fprintf('Receiving data file ''%s''\n', edfFile );
        status=Eyelink('ReceiveFile',edfFile,datadir,1);
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, datadir );
        elseif status == 0
            fprintf('Data file transfer cancelled\n');
        elseif status < 0
            fprintf('error in data file transfer\n');
        end
        if 2==exist(fullfile(datadir,edfFile), 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, datadir );
        end
    catch
        fprintf('Problem receiving data file ''%s''\n', edfFile );
    end
    Eyelink('Shutdown');
end

end

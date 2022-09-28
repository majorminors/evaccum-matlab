function spherefit = doMaxFilter(settings)
% modified version of Rik Henson's MaxFilter2pt2_CC


addpath('/neuro/meg_pd_1.2/'); % FIFACCESS toolbox if use some meg_misc functions below
addpath(genpath('/imaging/local/software/neuromag/'));
addpath(genpath('/hpc-software/neuromag/'));

cd(settings.maxfld); % you want to be wherever you want your files saved lest you lose them forever

%% 1 do SSS filtering

% aligning to the central run
nfiles = 1:length(settings.infname);
nfiles(1) = nfiles(ceil(numel(nfiles)/2));
nfiles(ceil(numel(nfiles)/2)) = 1;

save('settings.mat','settings');

count = 0; % since we're doing this from the central run, lets get a counter to report where we're up to
for i = nfiles
    
    count = count+1;
    fprintf('maxxxxxxing %.0f of %.0f raw files\n',count,numel(nfiles))
    
    infname  = settings.infname{i};
    outfname = settings.outfname{i};

    
    [fd, fstem, ~] = fileparts(outfname);
    
    if settings.overwrite %Do some cleaning otherwise Maxfilter will abort
        delf = dir([fd,'/',fstem,'*']);delf = {delf.name};
        delf = strcat([fd,'/'],delf);
        eval(sprintf('!rm -rf %s',strjoin(delf,' ')));
    end
    
    if exist([fstem,'_sphere_fit.txt'],'file') 
                eval(sprintf('!rm -rf %s',[fstem,'_sphere_fit.txt']));
    end
    
    spherefit = meg_fit_sphere_rik(infname,outfname);%might want to add EEG points

    posfname = [fstem '_headposition.pos']; % specifying head position filename
    
    % Initialise
    orgcmd = '';%#ok
    hpicmd = '';%#ok
    stcmd  = '';%#ok
    trcmd_def = '';
    
    % Origin and frame
    orgcmd = sprintf(' -frame head -origin %g %g %g',spherefit(1),spherefit(2),spherefit(3));
    
    % Autobad throughout (assumes 15mins=900s)
   
    %badstr  = sprintf(' -autobad %d -badlimit %d',900,7);
    
    badstr = settings.badchans;
   
    
    % SSS with ST:normal spatial filtering is applied in blocks of few secs
    % (buffers)
    
    stwin    = 10;%default data buffering = 10 secs; tSSS acts as a high-pass filter suppressing slow frequency components.
    %10 secs buffering corresponds to  a cut-off freq of 0.1Hz (1/10s). Longer
    %buffers can account for slower background variations but also increase memory
    %usage.
    
    stcorr   = 0.980;%the insight of SSSt is that if there are correlating waveforms in the internal and external signal, they must be artefact.
    stcmd    = sprintf(' -st %d -corr %g ',stwin,stcorr);
    %stcmd = ''; % To avoid jumps between end of 10s buffer, as warned in manual
    
    % HPI estimation and movement compensation
    hpistep = 10; hpisubt='amp';
    hpicmd = sprintf(' -linefreq 50 -hpistep %d -hpisubt %s -hpicons -movecomp inter -hp %s',hpistep,hpisubt,posfname);
    
    % If want to trans def in one go
    %outpfx    = 'transdef_';
    %trcmd_def = sprintf(' -trans default -frame head -origin %g %g %g',spherefit(1),spherefit(2)-13,spherefit(3)+6);
    
    % Preparing names of output and log files
    outsfx    = '_sss';
    outfname  = sprintf('%s%s.fif',fstem,outsfx);
    logfname  = sprintf('%s%s.log',fstem,outsfx);
    
    %store bad channels
    badfile  = ['bad',fstem,'.txt'];
    if exist(badfile,'file');delete(badfile);end
    %needs to be changed
    eval(sprintf('!cat %s | sed -n -e ''/Detected/p'' -e ''/Static/p'' | cut -f 5- -d '' '' > %s',logfname,badfile))
    % if exist(outfname)
    %      delete(outfname)
    % end
    
    % Skip if MF crashes if cHPI off for too long???? (consequences for later???)
    % if skipamount > 0
    %     skipstr = sprintf('  -skip 0 %d ',skipamount);
    % else
    skipstr = '';
    %end
    
    %mfcall = '! /neuro/bin/util/maxfilter-2.2'; % for eval
    %mfcall = '/neuro/bin/util/maxfilter-2.2'; % for unix
    %mfcall = '/imaging/local/software/neuromag/bin/util';
    mfcall = '/hpc-software/neuromag/bin/util/maxfilter-2.2.12';
    
    % Assembling MF command
    mfcmd_rest=[
        mfcall ' -f ' infname ' -o ' outfname,...
        '	 -ctc /neuro/databases/ctc/ct_sparse.fif' ' ',...
        '	 -cal /neuro/databases/sss/sss_cal.dat' ' ',...
        skipstr, badstr, orgcmd, stcmd, hpicmd, trcmd_def ' -v | tee ' logfname
        ];
    disp(mfcmd_rest);
    
    % Executing MF
    if ~exist(outfname,'file') %note that MaxFilter will return an error if a file already exists 
        %    eval(mfcmd_rest);
        %    unix(mfcmd_rest)
        [status, maxres] = unix(mfcmd_rest); % this stops screen-dumping?
        if status ~= 0
            disp(maxres)
            error('MaxFilter failed! FYI: this could happen because you are on an unlicenced node (only f, g, or h have licences) or sometimes if you run in desktop mode (reason unknown).')
        else
            disp(maxres);
            %plot head movements & fits
            savename = fstem;
            check_headmov({logfname})
            savefig(savename);
            close(gcf);
            
        end
        
        
        
    end
    
    
    
    %% 2. Trans mid - align coordinates to a chosen run (refname)
    
    if i ==nfiles(1); refname = outfname; end %align to middle run
    
    trcmd_def = sprintf(' -trans %s',refname);%refname needed to re-align within subject across sessions (usually the middle session is chosen)
    outsfx    = '_transmid';
    infname   = outfname;
    outfname  = sprintf('%s%s.fif',fstem,outsfx);
    logfname  = sprintf('%s%s.log',fstem,outsfx);
    
    % Assembling MF command
    mfcmd_rest=[
        mfcall ' -f ' infname ' -o ' outfname,...
        trcmd_def ' -force -v | tee ' logfname
        ];
    disp(mfcmd_rest);
    
    % Executing MF
    if ~exist(outfname,'file')
        %    eval(mfcmd_rest);
        %    unix(mfcmd_rest)
        [status, maxres] = unix(mfcmd_rest); % this stops screen-dumping?
        if status ~= 0
            error(sprintf('MaxFilter failed! %s',maxres))
        else
            disp(maxres);
        end
    end
    
    %% 3. Trans default (so that origin not same as SSS expansion origin above)
    
    trcmd_def = sprintf(' -trans default -origin %g %g %g -frame head ',spherefit(1),spherefit(2)-13,spherefit(3)+6);
    outsfx    = '_trans';
    infname   = outfname;
    outfname  = sprintf('%s%s.fif',fstem,outsfx);
    logfname  = sprintf('%s%s.log',fstem,outsfx);
    
    % Assembling MF command
    mfcmd_rest=[
        mfcall ' -f ' infname ' -o ' outfname,...
        trcmd_def ' -force -v | tee ' logfname
        ];
    disp(mfcmd_rest);
    
    % Executing MF
    if ~exist(outfname,'file')
        [status, maxres] = unix(mfcmd_rest); % this stops screen-dumping?
        if status ~= 0
            error('MaxFilter failed!');
        else
            disp(maxres);
            
        end
    end
    
    
end

return


function spherefit = meg_fit_sphere_rik(infname,outfname)

if nargin < 2
    outfname = infname;
end

% Assumes FIFACCESS tool hpipoints and Neuromag utility
% fit_sphere_to_points are on path

%% Read digitized points from FIF file
[pth, fstem, ~] = fileparts(outfname);
    
hptxtfname = fullfile(pth,[fstem '_headpoints.txt']); % specifying head points text file
if exist(hptxtfname,'file')~=2
    [co, ki] = hpipoints(infname); % HPI points  
    if size(co,2)~=3; co=co'; end % In case linux command above transposes on different machines (only fail if 3 headpoints too!)
    headpoints = co(ki>1,:);      % don't include the fiducial points
    headpoints = headpoints(find(~(headpoints(:,2)>0 & headpoints(:,3)<0)),:);      %#ok Remove nose points:
    save(hptxtfname,'-ASCII','headpoints');
end

%% Fitting sphere to headshape points
sphtxtfname = fullfile(pth,[fstem '_sphere_fit.txt']);
if exist(sphtxtfname,'file')~=2
    cmd_fit = ['/neuro/bin/util/fit_sphere_to_points ' hptxtfname];
    [status, spherefit] = unix(cmd_fit);
    if status ~= 0 || length(spherefit)<1
        error('Spherefit failed!')
    end
    spherefit = sscanf(spherefit,'%f')*1000;
    save(sphtxtfname,'-ASCII','spherefit');
else
    spherefit = textread(sphtxtfname,'%f');
end

return

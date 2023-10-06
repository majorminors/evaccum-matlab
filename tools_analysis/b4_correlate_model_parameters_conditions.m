%%%%%%%%%%%%
%% set up %%
%%%%%%%%%%%%

clear all

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';
datadir = fullfile(rootdir,'data','meg_pilot_4'); addpath(datadir);
preProcDir = fullfile(datadir,'%s','Preprocess'); % wants the id for the subject
toolsdir = fullfile(rootdir,'tools_analysis');
saveDir = fullfile(datadir,'model_correlations'); if ~exist(saveDir,'dir'); mkdir(saveDir); end
jobdir = fullfile(rootdir,'job_logging');
addpath(genpath(fullfile(toolsdir,'lib')))
ftDir = fullfile(toolsdir,'..','..','..','Toolboxes','fieldtrip');
addpath(ftDir); ft_defaults;

toolbox = fullfile(rootdir,'..','..','Toolboxes','gramm'); addpath(toolbox); clear toolbox
toolbox = fullfile(rootdir,'..','..','Toolboxes','BFF_repo'); addpath(genpath(toolbox)); clear toolbox

slopeTimeWindow = 10;

%%%%%%%%%%%%%%%%%%%%%%%%
%% load and prep data %%
%%%%%%%%%%%%%%%%%%%%%%%%

disp('loading data')

subjectParams = readtable(fullfile(rootdir,'dmc_lba','subjectParams.csv'));
meanParams = readtable(fullfile(rootdir,'dmc_lba','meanParams.csv'));

% create averages across manipulations and add them to the table
meanParams.b_ec = mean([meanParams.b_ecer, meanParams.b_echr],2);
meanParams.b_hc = mean([meanParams.b_hcer, meanParams.b_hchr],2);
meanParams.b_er = mean([meanParams.b_ecer, meanParams.b_hcer],2);
meanParams.b_hr = mean([meanParams.b_echr, meanParams.b_hchr],2);
meanParams.v1_ec = mean([meanParams.v1_ecer, meanParams.v1_echr],2);
meanParams.v1_hc = mean([meanParams.v1_hcer, meanParams.v1_hchr],2);
meanParams.v1_er = mean([meanParams.v1_ecer, meanParams.v1_hcer],2);
meanParams.v1_hr = mean([meanParams.v1_echr, meanParams.v1_hchr],2);
subjectParams.b_ec = mean([subjectParams.b_ecer, subjectParams.b_echr],2);
subjectParams.b_hc = mean([subjectParams.b_hcer, subjectParams.b_hchr],2);
subjectParams.b_er = mean([subjectParams.b_ecer, subjectParams.b_hcer],2);
subjectParams.b_hr = mean([subjectParams.b_echr, subjectParams.b_hchr],2);
subjectParams.v1_ec = mean([subjectParams.v1_ecer, subjectParams.v1_echr],2);
subjectParams.v1_hc = mean([subjectParams.v1_hcer, subjectParams.v1_hchr],2);
subjectParams.v1_er = mean([subjectParams.v1_ecer, subjectParams.v1_hcer],2);
subjectParams.v1_hr = mean([subjectParams.v1_echr, subjectParams.v1_hchr],2);

% load meeg data
varsToLoad = {...
    'ecer_responseLockedAverage','echr_responseLockedAverage','hcer_responseLockedAverage','hchr_responseLockedAverage',...
    'ec_responseLockedAverage','hc_responseLockedAverage','er_responseLockedAverage','hr_responseLockedAverage',...
    'ecer_coherenceLockedAverage','echr_coherenceLockedAverage','hcer_coherenceLockedAverage','hchr_coherenceLockedAverage',...
    'ec_coherenceLockedAverage','hc_coherenceLockedAverage','er_coherenceLockedAverage','hr_coherenceLockedAverage'...
    };

% alright, now we need to check which subjects we have data for, and which
% subject that dataset corresponds to in the DMC_table we fed into the LBA

% first we'll basically do a subset of what we did in b3_aggregate_data.m to
% determine where it skipped participants in producing timelocked averages

% load the subject information
allSubjects = importParticipants();

% loop through the subjects that have valid timelocked files
validSubjectNums = []; validSubjectIds = [];
dataTable = readtable(fullfile(datadir,'behavioural','DMC_table.csv'));
count = 0;
for subjectidx = 1:numel(allSubjects)
    thisSubject = allSubjects{subjectidx};
    %disp(thisSubject.id)
    timelockedFilename = fullfile(sprintf(preProcDir,thisSubject.id),'timelocked_averages.mat');
    if ~exist(timelockedFilename,'file')
        %warning('no file exists for this participant, skipping')
        continue
    elseif exist(timelockedFilename,'file')
        count = count+1;
        validSubjectIds = [validSubjectIds;{thisSubject.id}];
        % so we get the subject numbers which will correspond to the
        % subject numbers saved in the DMC file, but we want them in the
        % same format for comparison so convert to string and put in cells
        validSubjectNums = [validSubjectNums;{num2str(thisSubject.num)}];
        %         data{count} = load(timelockedFilename, varsToLoad{:});
    end
    clear timelockedFilename thisSubject
end; clear thisSubject allSubjects subjectidx

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    disp('no parallel pool is currently initialised---initialising');
    workers = numel(validSubjectIds);
    P=cbupool(workers);
    P.JobStorageLocation = jobdir;
    tempPath = fullfile(jobdir,'tmp');
    if ~exist(tempPath,'dir');mkdir(tempPath);end
    P.JobStorageLocation = tempPath;
    parpool(P,workers);
else
    disp('parallel pool has already been initialized---skipping');
end
parfor subject = 1:numel(validSubjectIds)
    thisSubject = validSubjectIds{subject};
    disp(['loading subject: ' thisSubject])
    timelockedFilename = fullfile(sprintf(preProcDir,thisSubject),'timelocked_averages.mat');
    ftData{subject} = load(timelockedFilename, varsToLoad{:});
end
disp('subject fieldtrip data loaded')
delete(gcp('nocreate')); clear workers;
rmdir(tempPath,'s');

% ok, let's clear out the invalid subjects so we can just loop through
validSubjectParams = subjectParams(ismember(subjectParams.Var1,validSubjectNums),:);

% set our channels
CPP = {'EEG040' 'EEG041' 'EEG042'};
%%
disp('collating data')

for subjectidx = 1:size(validSubjectParams,1)
    
    % first we'll get the parameter values
    paramTemplates = {'v1_%s' 'b_%s'};
    count = 0;
    splitStr = split(varsToLoad, '_');
    conds = splitStr(:,~contains(splitStr(:,:,2),'response'),1); % just need one of the two of these, for parameters
    clear splitStr
    param = table('Size', [numel(conds) 3], 'VariableTypes', {'string', 'string', 'double'},...
        'VariableNames', {'ParamName', 'Condition', 'ParamValue'});
    for thisCond = conds
        for thisParam = paramTemplates
            count = count+1;
            thisParam = sprintf(thisParam{:},thisCond{:});
            titleBits = split(thisParam, '_');
            
            param(count,:) = {titleBits{1}, titleBits{2}, validSubjectParams.(thisParam)(subjectidx)};
            
        end; clear thisParam titleBits
    end; clear thisCond conds count
    
    % now grab the amplitudes
    amplitude = table('Size', [numel(varsToLoad) 4], 'VariableTypes', {'string', 'string', 'cell', 'double'},...
        'VariableNames', {'LockedTo', 'Condition', 'Values', 'TimePoints'});
    count = 0;
    for thisVar = varsToLoad
        count = count+1;
        splitStr = split(thisVar, '_');
        condition = splitStr(1); condition = condition{:};
        lockedTo = splitStr(2);
        lockedTo = strrep(lockedTo, 'Average', ''); lockedTo = lockedTo{:};
        clear splitStr
        
        amplitude(count,:) = {lockedTo, condition, {getRawAmplitude(ftData{subjectidx}.(thisVar{:}),CPP)}, numel(getRawAmplitude(ftData{subjectidx}.(thisVar{:}),CPP))};
    end; clear count condition lockedTo thisVar
    
    % then create the local changes in slope
    slope = table('Size', [size(amplitude,1) 4], 'VariableTypes', {'string', 'string', 'cell', 'double'}, 'VariableNames', {'LockedTo', 'Condition', 'Slope', 'TimePoints'});
    for i = 1:size(amplitude,1)
        timepoints = amplitude.TimePoints(i);
        amplitudeData = amplitude.Values{i};
        lockedTo = amplitude.LockedTo{i};
        condition = amplitude.Condition{i};
        slopeData = zeros(1, timepoints);
        
        for t = 1:timepoints
            if t<slopeTimeWindow+1
                slopeData(t) = amplitudeData(1,t) - mean(amplitudeData(1,:)); % get something in a sensible range for visualisation
            else
                slopeData(t) = amplitudeData(1,t) - amplitudeData(1,t-slopeTimeWindow); % get local change in slope as time - time-10
            end
        end; clear t
        
        slope(i,:) = {lockedTo, condition, {slopeData}, timepoints};
        clear slopeData
    end; clear i
    
    data(subjectidx).param = param;
    data(subjectidx).amplitude = amplitude;
    data(subjectidx).slope = slope;
    clear param amplitude slope
    
end; clear subjectidx

disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run the correlations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('getting correlations')

correlationsSaveName = [saveDir filesep 'model_correlations_conditions_01.mat'];

if exist(correlationsSaveName,'file')
    disp('file exists: loading')
    load(correlationsSaveName)
else
    disp('file does not exist: running analysis')
    
    varNames = {...
        'Timepoint',...
        'ParameterName', 'ParameterCondition',...
        'Condition', 'LockedTo',...
        'Amplitudes', 'Slopes', 'ParamVals'...
        'R_Amplitude', 'R_Slope'};
    correlations = table('Size', [0 10], ...
        'VariableTypes', {...
        'double',...
        'string', 'string',...
        'string', 'string',...
        'cell', 'cell', 'cell',...
        'double', 'double'}, ...
        'VariableNames', varNames);
    
    parameters = unique([data(1).param.ParamName])';
    parameterConds = unique([data(1).param.Condition])';
    lockedTos = unique([data(1).amplitude.LockedTo])';
    conditions = unique([data(1).amplitude.Condition])';
    
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        disp('no parallel pool is currently initialised---initialising');
        workers = 100;
        P=cbupool(workers);
        P.JobStorageLocation = jobdir;
        tempPath = fullfile(jobdir,'tmp');
        if ~exist(tempPath,'dir');mkdir(tempPath);end
        P.JobStorageLocation = tempPath;
        parpool(P,workers);
    else
        disp('parallel pool has already been initialized---skipping');
    end
    count = 0;
    for thisParameter = parameters
        for thisLockedTo = lockedTos
            for thisCondition = conditions
                thisParameterCond = thisCondition; % in case you want to run across difference parameter conds later, you can just add another for loop
                count=count+1; fprintf('%.0f of %.0f\n',count,numel(parameters)*numel(lockedTos)*numel(conditions))
                fprintf('correlating param %s_%s for condition %s_%s\n',thisParameter,thisParameterCond,thisCondition,thisLockedTo)
                
                % init these to fetch the amplitudes and slopes
                amplitudes = [];
                slopes = [];
                params = [];
                for thisSubject = 1:numel(data)
                    
                    params(thisSubject,1) = data(thisSubject).param.ParamValue(strcmp(data(thisSubject).param.ParamName, thisParameter) & strcmp(data(thisSubject).param.Condition, thisParameterCond));
                    
                    slopeIdx = strcmp(data(thisSubject).slope.LockedTo,thisLockedTo) & strcmp(data(thisSubject).slope.Condition,thisCondition);
                    ampIdx = strcmp(data(thisSubject).amplitude.LockedTo,thisLockedTo) & strcmp(data(thisSubject).amplitude.Condition,thisCondition);
                    
                    theseTimepoints = data(thisSubject).slope.TimePoints(slopeIdx);
                    if theseTimepoints ~= data(thisSubject).amplitude.TimePoints(ampIdx)
                        error('timepoints for slope and amplitude not equal---you are using the same timepoints for both, so correct this')
                    end
                    parfor t = 1:theseTimepoints
                        timepoints(thisSubject,t) = t;
                        slopes(thisSubject,t) = data(thisSubject).slope.Slope{slopeIdx}(t);
                        amplitudes(thisSubject,t) = data(thisSubject).amplitude.Values{ampIdx}(t);
                    end; clear t
                    
                end; clear thisSubject
                
                if any(mean(timepoints) ~= 1:theseTimepoints); error('timepoints not equal across participants?'); end
                timepoints = mean(timepoints);
                parfor t = timepoints
                    % calculate correlations
                    timepointAmps{t} = amplitudes(:,t);
                    timepointSlopes{t} = slopes(:,t);
                    timepointParams{t} = params;
                    r_amplitude(t) = corr(amplitudes(:,t), params, 'type', 'Pearson');
                    r_slope(t) = corr(slopes(:,t), params, 'type', 'Pearson');
                end; clear t
                
                correlations = [correlations;...
                    table(timepoints',repmat(thisParameter,numel(timepoints),1),repmat(thisParameterCond,numel(timepoints),1),...
                    repmat(thisCondition,numel(timepoints),1),repmat(thisLockedTo,numel(timepoints),1),...
                    timepointAmps',timepointSlopes',timepointParams',r_amplitude',r_slope','VariableNames', varNames)
                    ];
                clear timepoints timepointAmps timepointSlopes timepointParams r_amplitude r_slope
                
            end; clear thisCondition
        end; clear thisLockedTo
    end; clear thisParameter
    clear count
    %     delete(gcp('nocreate')); clear workers;
    rmdir(tempPath,'s');
    
    
    save(correlationsSaveName,'correlations')
    
end

disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run the bayes analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('getting bayes factors')

% nullInterval = '-0.2,0.2';
nullInterval = '0.2,1';
bfSlopeSavename = [saveDir filesep 'model_correlations_conditionsOnly_slope_bfs_null_%s.mat'];
bfAmplitudeSavename = [saveDir filesep 'model_correlations_conditionsOnly_amplitude_bfs_null_%s.mat'];
if exist(sprintf(bfSlopeSavename,nullInterval),'file')
    disp('file exists for slope bfs: loading')
    load(sprintf(bfSlopeSavename,nullInterval));
    runSlope = 0;
else
    runSlope = 1;
end
if exist(sprintf(bfAmplitudeSavename,nullInterval),'file')
    disp('file exists for amplitude bfs: loading')
    load(sprintf(bfAmplitudeSavename,nullInterval));
    runAmp = 0;
else
    runAmp = 1;
end

if runSlope || runAmp
    
    disp(['not all bf files exist (amplitude: ' num2str(runAmp) ', slope: ' num2str(runSlope) '): calculating'])
    
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        disp('no parallel pool is currently initialised---initialising');
        workers = 100;
        P = cbupool(workers);
        P.JobStorageLocation = jobdir;
        tempPath = fullfile(jobdir,'tmp');
        if ~exist(tempPath,'dir');mkdir(tempPath);end
        P.JobStorageLocation = tempPath;
        parpool(P,workers);
    else
        disp('parallel pool has already been initialized---skipping');
    end
    
    parfor i = 1:numel(correlations.Timepoint)
        fprintf('%.0f of %.0f\n',i,numel(correlations.Timepoint))
        
        % add the R module and get the path to Rscript
        [status, result] = system('module add R && which Rscript');
        if status == 0
            %         disp('R module added successfully');
            RscriptPath = strtrim(result);
            %         disp(['Rscript path: ' RscriptPath]);
        else
            error('Failed to add R and/or locate Rscript');
        end
        
        theseArgs = ['rscale="medium",nullInterval=c(' nullInterval ')'];
        
        if runSlope
            % now run the rscript version of the bayes analysis
            sbf1(i) = bayesfactor_R_wrapper_corr(correlations.ParamVals{i}', correlations.Slopes{i}','Rpath',RscriptPath,'returnindex',1,...
                'args',theseArgs,'verbose',0);
            sbf2(i) = bayesfactor_R_wrapper_corr(correlations.ParamVals{i}', correlations.Slopes{i}','Rpath',RscriptPath,'returnindex',2,...
                'args',theseArgs,'verbose',0);
        end
        if runAmp
            % now run the rscript version of the bayes analysis
            abf1(i) = bayesfactor_R_wrapper_corr(correlations.ParamVals{i}', correlations.Amplitudes{i}','Rpath',RscriptPath,'returnindex',1,...
                'args',theseArgs,'verbose',0);
            abf2(i) = bayesfactor_R_wrapper_corr(correlations.ParamVals{i}', correlations.Amplitudes{i}','Rpath',RscriptPath,'returnindex',2,...
                'args',theseArgs,'verbose',0);
        end
        
    end
    if runSlope; slopeBFs.bf1 = sbf1; slopeBFs.bf2 = sbf2; save(sprintf(bfSlopeSavename,nullInterval),'slopeBFs'); end
    if runAmp; amplitudeBFs.bf1 = abf1; amplitudeBFs.bf2 = abf2; save(sprintf(bfAmplitudeSavename,nullInterval),'amplitudeBFs'); end
    delete(gcp('nocreate')); clear workers;
    rmdir(tempPath,'s');
    clear runSlope runAmp
end

correlations.SlopeBf1 = slopeBFs.bf1';
correlations.SlopeBf2 = slopeBFs.bf2';
correlations.AmplitudeBf1 = amplitudeBFs.bf1';
correlations.AmplitudeBf2 = amplitudeBFs.bf2';

disp('done')

%%%%%%%%%%
%% plot %%
%%%%%%%%%%

disp(nullInterval)
plotInside = 1; % plot bf for inside the null interval (1) or outside (0)?
plotSlope = 0; % plot slope (1) or amplitude (0)

disp('plotting')

for thisParam = unique(correlations.ParameterName)'
    for thisCond = unique(correlations.Condition)'
        for thisLockedTo = unique(correlations.LockedTo)'
            paramCond = thisCond; % if you want to look cross cond, you can add this back as a for loop
            
            figure;
            
            dataidx = ...
                strcmp(correlations.ParameterCondition,paramCond)...
                & strcmp(correlations.ParameterName,thisParam)...
                & strcmp(correlations.Condition,thisCond)...
                & strcmp(correlations.LockedTo,thisLockedTo);
            
            if plotSlope
                if plotInside
                    thesebfs = correlations.SlopeBf1(dataidx);
                else
                    thesebfs = correlations.SlopeBf2(dataidx);
                end
                theseRs = correlations.R_Slope(dataidx);
                allRs = correlations.R_Slope;
                tipo = 'Slope';
            elseif ~plotSlope
                if plotInside
                    thesebfs = correlations.AmplitudeBf1(dataidx);
                else
                    thesebfs = correlations.AmplitudeBf2(dataidx);
                end
                theseRs = correlations.R_Amplitude(dataidx);
                allRs = correlations.R_Amplitude;
                tipo = 'Amplitude';
            end
            
            slopeInvalidColour = [1, 0.6, 0.6];
            
            timepoints = correlations.Timepoint(dataidx);
            
            if contains(thisLockedTo,'response')
                onset = 600; % to subtract from timepoints
                xlims = [-600 200];
            elseif contains(thisLockedTo,'coherence')
                onset = 500; % to subtract from timepoints
                xlims = [-500 1500];
            end
            timepoints = timepoints-onset;
            
            % plot correlations
            subplot(2,1,1)
            color_map = summer(length(theseRs));
            color_map = color_map(:, [2, 1, 3]);  % Rearrange color map channels
            abs_Rs = abs(theseRs);
            scatter_colours = 1 - abs_Rs / max(abs_Rs);  % Scaling the color values to make darker as value moves away from 0
            scatter(timepoints, theseRs, [], scatter_colours, 'filled');
            colormap(color_map);
            %         cbh = colorbar;
            xlim(xlims)
            ylim([min(allRs)-0.1 max(allRs)+0.1])
            line(get(gca,'XLim'), [0 0], 'Color', [0.5 0.5 0.5], 'LineStyle', '-')
            line([0 0], get(gca,'YLim'), 'Color', [0.1 0.1 0.1], 'LineStyle', '--')
            if strcmp(tipo,'Slope')
                % plot a red bar over invalid timepoints
                y = get(gca,'YLim');
                x = get(gca,'XLim');
                patch([x(1) x(1) x(1)+slopeTimeWindow x(1)+slopeTimeWindow],...
                    [y(1) y(2) y(2) y(1)],...
                    slopeInvalidColour,'FaceAlpha',0.5);
                clear y x
            end
            xlabel('Timepoints');
            ylabel('Correlations (r)');
            title(['Correlation: ' thisParam{:} ' (' paramCond{:} ') and ' tipo ' of ' thisLockedTo{:} ' (' thisCond{:} ')'])
            
            % plot bayes factors
            subplot(2,1,2)
            load('bayes_colourmap.mat'); % in BFF repo
            exponential_minmax=2;
            val_col_map = logspace(-exponential_minmax,exponential_minmax,size(colours,1));
            scatter_colours = zeros(length(timepoints), 3);  % preallocate for efficiency
            for t = 1:length(timepoints)
                [~,idx] = min(abs(val_col_map-thesebfs(t)));
                scatter_colours(t, :) = colours(idx,1:3);
            end
            scatter(timepoints, thesebfs, 30, scatter_colours, 'filled');
            line(get(gca,'XLim'),[1 1], 'Color', [0.5 0.5 0.5], 'LineStyle', '--')
            ax = gca;
            set(ax,'YScale','log','XLim',[timepoints(1),timepoints(end)], ...
                'YLim',[1e-2 1e2],'YTick',10.^(-3:1:3))
            if strcmp(tipo,'Slope')
                % plot a red bar over invalid timepoints
                y = get(gca,'YLim');
                x = get(gca,'XLim');
                patch([x(1) x(1) x(1)+slopeTimeWindow x(1)+slopeTimeWindow],...
                    [y(1) y(2) y(2) y(1)],...
                    slopeInvalidColour,'FaceAlpha',0.5);
                clear y x
            end
            xlabel('Time (s)')
            ylabel('BF (log scale)')
            colormap(ax,colours)
            cbh = colorbar;
            caxis([-exponential_minmax,exponential_minmax])
            cbh.Units = 'normalized';
            cbh.Limits = [-exponential_minmax,exponential_minmax];
            cbh.Position(1) = 0.92;cbh.Position(3) = 0.01;cbh.Position(4) = ax.Position(4);cbh.Position(2) = ax.Position(2);
            cbh.Label.String = 'Bayes Factor';
            f = gcf; f.Position = [10 10 1600 1600];
            pause(1)
            cbh.Ticks = [-2, -1, -0.5, 0, 0.5, 1, 2];
            cbh.TickLabels=arrayfun(@(x) ['10^{' num2str(x) '}'], cbh.Ticks, 'UniformOutput', false);
            cbh.TickLabels(strcmp(cbh.TickLabels,'10^{0}')) = {'Inconclusive'};
            cbh.TickLabels(strcmp(cbh.TickLabels,'10^{0.5}') | strcmp(cbh.TickLabels,'10^{-0.5}')) = {'Moderate'};
            cbh.TickLabels(strcmp(cbh.TickLabels,'10^{1}') | strcmp(cbh.TickLabels,'10^{-1}')) = {'Strong'};
            
            print([saveDir filesep tipo '_' thisParam{:} '_' paramCond{:} '_' thisCond{:} '_' thisLockedTo{:} '.png'],'-dpng')
            
        end; clear thisLockedTo
    end; clear thisCond
end; clear thisParam
close all

disp('done')

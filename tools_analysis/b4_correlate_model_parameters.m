%%%%%%%%%%%%
%% set up %%
%%%%%%%%%%%%

clear all

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';
datadir = fullfile(rootdir,'data','meg_pilot_4'); addpath(datadir);
preProcDir = fullfile(datadir,'%s','Preprocess'); % wants the id for the subject
toolsdir = fullfile(rootdir,'tools_analysis');
saveDir = fullfile(datadir,'model_correlations'); if ~exist(saveDir,'dir'); mkdir(saveDir); end
addpath(genpath(fullfile(toolsdir,'lib')))
ftDir = fullfile(toolsdir,'..','..','..','Toolboxes','fieldtrip');
addpath(ftDir); ft_defaults;

toolbox = fullfile(rootdir,'..','..','Toolboxes','gramm'); addpath(toolbox); clear toolbox
toolbox = fullfile(rootdir,'..','..','Toolboxes','BFF_repo'); addpath(genpath(toolbox)); clear toolbox


%%%%%%%%%%%%%%%%%%%%%%%%
%% load and prep data %%
%%%%%%%%%%%%%%%%%%%%%%%%

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

averages = {'cOnsDiffAve', 'cRespDiffAve', 'rOnsDiffAve', 'rRespDiffAve'};
subjectwise = {'cOnsDiffAll', 'cRespDiffAll', 'rOnsDiffAll', 'rRespDiffAll'};

load(fullfile(datadir,'differences.mat'), subjectwise{:})
%load(fullfile(datadir,'differences.mat'), averages{:})

clear averages subjectwise

% alright, now we need to check which subjects we have data for, and which
% subject that dataset corresponds to in the DMC_table we fed into the LBA

% first we'll basically do a subset of what we did in b3_aggregate_data.m to
% determine where it skipped participants in producing timelocked averages

% load the subject information
allSubjects = importParticipants();

% loop through the subjects that have valid timelocked files
validSubjects = [];
dataTable = readtable(fullfile(datadir,'behavioural','DMC_table.csv'));
for subjectidx = 1:numel(allSubjects)
    thisSubject = allSubjects{subjectidx};
    %disp(thisSubject.id)
    timelockedFilename = fullfile(sprintf(preProcDir,thisSubject.id),'timelocked_averages.mat');
    if ~exist(timelockedFilename,'file')
        %warning('no file exists for this participant, skipping')
        continue
    elseif exist(timelockedFilename,'file')
        % so we get the subject numbers which will correspond to the
        % subject numbers saved in the DMC file, but we want them in the
        % same format for comparison so convert to string and put in cells
        validSubjects = [validSubjects;{num2str(thisSubject.num)}];
    end
    clear timelockedFilename thisSubject
end; clear thisSubject allSubjects subjectidx

% ok, let's clear out the invalid subjects so we can just loop through
validSubjectParams = subjectParams(ismember(subjectParams.Var1,validSubjects),:);

% set our channels
CPP = {'EEG040' 'EEG041' 'EEG042'};

for subjectidx = 1:size(validSubjectParams,1)
    param = table('Size', [4 3], 'VariableTypes', {'string', 'string', 'double'},...
        'VariableNames', {'ParamName', 'Condition', 'ParamValue'});
    param(1,:) = {'Decision Bdry', 'Coh', validSubjectParams.b_ec(subjectidx)-validSubjectParams.b_hc(subjectidx)};
    param(2,:) = {'Decision Bdry', 'Cat', validSubjectParams.b_er(subjectidx)-validSubjectParams.b_hr(subjectidx)};
    param(3,:) = {'Drift Rate', 'Coh', validSubjectParams.v1_ec(subjectidx)-validSubjectParams.v1_hc(subjectidx)};
    param(4,:) = {'Drift Rate', 'Cat', validSubjectParams.v1_er(subjectidx)-validSubjectParams.v1_hr(subjectidx)};

    amplitude = table('Size', [4 4], 'VariableTypes', {'string', 'string', 'cell', 'double'},...
        'VariableNames', {'LockedTo', 'Condition', 'Values', 'TimePoints'});
    amplitude(1,:) = {'Onset', 'Coherence', {getRawAmplitude(cOnsDiffAll{subjectidx},CPP)}, numel(getRawAmplitude(cOnsDiffAll{subjectidx},CPP))};
    amplitude(2,:) = {'Onset', 'Categorisation', {getRawAmplitude(rOnsDiffAll{subjectidx},CPP)}, numel(getRawAmplitude(rOnsDiffAll{subjectidx},CPP))};
    amplitude(3,:) = {'Response', 'Coherence', {getRawAmplitude(cRespDiffAll{subjectidx},CPP)}, numel(getRawAmplitude(cRespDiffAll{subjectidx},CPP))};
    amplitude(4,:) = {'Response', 'Categorisation', {getRawAmplitude(rRespDiffAll{subjectidx},CPP)}, numel(getRawAmplitude(rRespDiffAll{subjectidx},CPP))};

    slope = table('Size', [size(amplitude,1) 4], 'VariableTypes', {'string', 'string', 'cell', 'double'}, 'VariableNames', {'LockedTo', 'Condition', 'Slope', 'TimePoints'});

    for i = 1:size(amplitude,1)
        timepoints = amplitude.TimePoints(i);
        amplitudeData = amplitude.Values{i};
        lockedTo = amplitude.LockedTo{i};
        condition = amplitude.Condition{i};
        slopeData = zeros(1, timepoints);

        for t = 1:timepoints
            if t<11
                slopeData(t) = amplitudeData(1,t) - mean(amplitudeData(1,:)); % get something in a sensible range for visualisation
            else
                slopeData(t) = amplitudeData(1,t) - amplitudeData(1,t-10); % get local change in slope as time - time-10
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

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run the correlations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

correlationsSaveName = [saveDir filesep 'model_correlations_02.mat'];

if exist(correlationsSaveName,'file')
    load(correlationsSaveName)
else
    
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
    parpool(P,workers);
else
    disp('parallel pool has already been initialized---skipping');
end
    for thisParameter = parameters
        for thisParameterCond = parameterConds
            for thisLockedTo = lockedTos
                for thisCondition = conditions
                    
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
        end; clear thisParameterCond
    end; clear thisParameter
    
    save(correlationsSaveName,'correlations')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run the bayes analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    disp('no parallel pool is currently initialised---initialising');
    workers = 100;
    P=cbupool(workers);
    parpool(P,workers);
else
    disp('parallel pool has already been initialized---skipping');
end

nullInterval = '-0.3,0.3';
bfSlopeSavename = [saveDir filesep 'model_correlations_slope_bfs_null_%s.mat'];
bfAmplitudeSavename = [saveDir filesep 'model_correlations_amplitude_bfs_null_%s.mat'];
if exist(sprintf(bfSlopeSavename,nullInterval),'file')
    load(sprintf(bfSlopeSavename,nullInterval));
    runSlope = 0;
else
    runSlope = 1;
end
if exist(sprintf(bfAmplitudeSavename,nullInterval),'file')
    load(sprintf(bfAmplitudeSavename,nullInterval));
    runAmp = 0;
else
    runAmp = 1;
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
% % delete(gcp('nocreate')); clear workers
clear runSlope runAmp

correlations.SlopeBf1 = slopeBFs.bf1';
correlations.SlopeBf2 = slopeBFs.bf2';
correlations.AmplitudeBf1 = amplitudeBFs.bf1';
correlations.AmplitudeBf2 = amplitudeBFs.bf2';

%%%%%%%%%%
%% plot %%
%%%%%%%%%%

for paramCond = unique(correlations.ParameterCondition)'
    for thisParam = unique(correlations.ParameterName)'
        for thisCond = unique(correlations.Condition)'
            for thisLockedTo = unique(correlations.LockedTo)'

            figure;

            dataidx = ...
                strcmp(correlations.ParameterCondition,paramCond)...
                & strcmp(correlations.ParameterName,thisParam)...
                & strcmp(correlations.Condition,thisCond)...
                & strcmp(correlations.LockedTo,thisLockedTo);

            
% %             thesebfs = correlations.SlopeBf1(dataidx);
%             thesebfs = correlations.SlopeBf2(dataidx);
%             theseRs = correlations.R_Slope(dataidx);
%             allRs = correlations.R_Slope;
%             tipo = 'Slope';
            
% %             thesebfs = correlations.AmplitudeBf1(dataidx);
            thesebfs = correlations.AmplitudeBf2(dataidx);
            theseRs = correlations.R_Amplitude(dataidx);
            allRs = correlations.R_Amplitude;
            tipo = 'Amplitude';
            
            timepoints = correlations.Timepoint(dataidx);
            
            if contains(thisLockedTo,'Response')
                onset = 600; % to subtract from timepoints
                xlims = [-600 200];
            elseif contains(thisLockedTo,'Onset')
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
            xlabel('Timepoints');
            ylabel('Correlations (r)');
            title(['Correlation: ' thisParam{:} ' (' paramCond{:} ') and ' tipo ' of '  thisCond{:} ' ' thisLockedTo{:}])

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
end; clear paramCond
close all

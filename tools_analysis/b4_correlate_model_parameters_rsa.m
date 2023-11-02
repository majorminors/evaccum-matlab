%%%%%%%%%%%%
%% set up %%
%%%%%%%%%%%%

clear all

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';
datadir = fullfile(rootdir,'data','meg_pilot_4'); addpath(datadir);
subDir = fullfile(datadir,'%s'); % wants the id for the subject
toolsdir = fullfile(rootdir,'tools_analysis');
saveDir = fullfile(datadir,'model_correlations'); if ~exist(saveDir,'dir'); mkdir(saveDir); end
jobdir = fullfile(rootdir,'job_logging');
addpath(genpath(fullfile(toolsdir,'lib')))
% ftDir = fullfile(toolsdir,'..','..','..','Toolboxes','fieldtrip');
% addpath(ftDir); ft_defaults;

toolbox = fullfile(rootdir,'..','..','Toolboxes','gramm'); addpath(toolbox); clear toolbox
toolbox = fullfile(rootdir,'..','..','Toolboxes','BFF_repo'); addpath(genpath(toolbox)); clear toolbox


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
    %disp(thisSubject.id)s
    thisFile = fullfile(sprintf(subDir,thisSubject.id),'rsa.mat');
    if ~exist(thisFile,'file')
        %warning('no file exists for this participant, skipping')
        continue
    elseif exist(thisFile,'file')
        count = count+1;
        validSubjectIds = [validSubjectIds;{thisSubject.id}];
        % so we get the subject numbers which will correspond to the
        % subject numbers saved in the DMC file, but we want them in the
        % same format for comparison so convert to string and put in cells
        validSubjectNums = [validSubjectNums;{num2str(thisSubject.num)}];
        %         data{count} = load(timelockedFilename, varsToLoad{:});
    end
    clear thisFile thisSubject
end; clear thisSubject allSubjects subjectidx

for subject = 1:numel(validSubjectIds)
    thisSubject = validSubjectIds{subject};
    disp(['loading subject: ' thisSubject])
    thisFile = fullfile(sprintf(subDir,thisSubject),'rsa.mat');
    tmp = load(thisFile, 'rsa');
    rsa{subject} = tmp.rsa;
end
disp('subject data loaded')

% ok, let's clear out the invalid subjects so we can just loop through
validSubjectParams = subjectParams(ismember(subjectParams.Var1,validSubjectNums),:);

%%
disp('collating data')

for subjectidx = 1:size(validSubjectParams,1)
    
    % first we'll get the parameter values
    paramTemplates = {'v1_%s' 'b_%s'};
    count = 0;
    splitStr = split(fieldnames(rsa{subjectidx}), '_');
    conds = splitStr(~contains(splitStr(:,2),'response'),1); % just need one of the two of these, for parameters
    clear splitStr
    param = table('Size', [numel(conds) 3], 'VariableTypes', {'string', 'string', 'double'},...
        'VariableNames', {'ParamName', 'Condition', 'ParamValue'});
    for thisCond = conds'
        for thisParam = paramTemplates
            count = count+1;
            thisParam = sprintf(thisParam{:},thisCond{:});
            titleBits = split(thisParam, '_');
            
            param(count,:) = {titleBits{1}, titleBits{2}, validSubjectParams.(thisParam)(subjectidx)};
            
        end; clear thisParam titleBits
    end; clear thisCond conds count
    
    % now grab the rsaCorrs
    rsaCorrs = table('Size', [numel(fieldnames(rsa{subjectidx})) 5], 'VariableTypes', {'string', 'string', 'string', 'cell', 'double'},...
        'VariableNames', {'LockedTo', 'Condition', 'Model', 'Values', 'TimePoints'});
    count = 0;
    for thisVar = fieldnames(rsa{subjectidx})'
        for thisModel = fieldnames(rsa{subjectidx}.(thisVar{:}))'
            count = count+1;
            splitStr = split(thisVar, '_');
            condition = splitStr(1); condition = condition{:};
            lockedTo = splitStr(2);
            lockedTo = strrep(lockedTo, 'Average', ''); lockedTo = lockedTo{:};
            clear splitStr
            
            whichRsaCorrs = 'partial_samples'; % samples partial_samples fisher_transformed_samples fisher_transformed_partial_samples
            
            rsaCorrs(count,:) = {lockedTo, condition, thisModel {rsa{subjectidx}.(thisVar{:}).(thisModel{:})}, numel(rsa{subjectidx}.(thisVar{:}).(thisModel{:}).(whichRsaCorrs))};
        end
    end; clear count condition lockedTo thisVar thisModel
    
    
    data(subjectidx).param = param;
    data(subjectidx).rsaCorrs = rsaCorrs;
    clear param rsaCorrs slope
    
end; clear subjectidx

disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run the correlations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('getting correlations')

correlationsSaveName = [saveDir filesep 'model_correlations_rsa_01.mat'];

if exist(correlationsSaveName,'file')
    warning('file exists');
    response = input('Do you want to overwrite? (y/n): ','s');
    if strcmpi(response,'y') || strcmpi(response,'yes')
        disp('removing file...');
        system(['rm -rf ' correlationsSaveName]);
        if exist(correlationsSaveName,'file'); error('I couldnt delete it :('); end
        error('ok, good to go, run this section again');
    else
        disp('loading')
        load(correlationsSaveName)
    end
else
    disp('file does not exist: running analysis')
    
    varNames = {...
        'Timepoint',...
        'ParameterName', 'ParameterCondition',...
        'Condition', 'LockedTo',...
        'RsaCorrelations', 'ParamVals'...
        'R_Rsa', 'Model'};
    correlations = table('Size', [0 9], ...
        'VariableTypes', {...
        'double',...
        'string', 'string',...
        'string', 'string',...
        'cell', 'cell',...
        'double', 'string'}, ...
        'VariableNames', varNames);
    
    parameters = unique([data(1).param.ParamName])';
    parameterConds = unique([data(1).param.Condition])';
    lockedTos = unique([data(1).rsaCorrs.LockedTo])';
    conditions = unique([data(1).rsaCorrs.Condition])';
    
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
                clear rsaCorrs
                thisParameterCond = thisCondition; % in case you want to run across difference parameter conds later, you can just add another for loop
                count=count+1; fprintf('%.0f of %.0f\n',count,numel(parameters)*numel(lockedTos)*numel(conditions))
                fprintf('correlating param %s_%s for condition %s_%s\n',thisParameter,thisParameterCond,thisCondition,thisLockedTo)
                
                % init these to fetch the rsaCorrs
                params = [];
                for thisSubject = 1:numel(data)
                    
                    params(thisSubject,1) = data(thisSubject).param.ParamValue(strcmp(data(thisSubject).param.ParamName, thisParameter) & strcmp(data(thisSubject).param.Condition, thisParameterCond));
                    
                    theseTimepoints = data(thisSubject).rsaCorrs.TimePoints(strcmp(data(thisSubject).rsaCorrs.LockedTo,thisLockedTo) & strcmp(data(thisSubject).rsaCorrs.Condition,thisCondition));
                    if numel(unique(theseTimepoints)) > 1; error('you have different timepoints'); else; theseTimepoints = theseTimepoints(1); end
                    
                    models = unique(data(thisSubject).rsaCorrs.Model)';
                    
                    modelCount = 0;
                    for thisModel = models
                        modelCount = modelCount+1;
                        rsaIdx = strcmp(data(thisSubject).rsaCorrs.LockedTo,thisLockedTo) & strcmp(data(thisSubject).rsaCorrs.Condition,thisCondition) & strcmp(data(thisSubject).rsaCorrs.Model,thisModel);
                        rsaCorrs(thisSubject,:,modelCount) = data(thisSubject).rsaCorrs.Values{rsaIdx}.(whichRsaCorrs);
                    end; clear thisModel modelCount
                    
                end; clear thisSubject
                
                modelCount = 0;
                for thisModel = models
                    modelCount = modelCount+1;
                    timepoints = 1:theseTimepoints;
                    parfor t = timepoints
                        % calculate correlations
                        timepointRsa{t} = rsaCorrs(:,t,modelCount);
                        timepointParams{t} = params;
                        r_rsaCorrs(t) = corr(rsaCorrs(:,t,modelCount), params, 'type', 'Pearson');
                    end; clear t
                                        
                    correlations = [correlations;...
                        table(timepoints',repmat(thisParameter,numel(timepoints),1),repmat(thisParameterCond,numel(timepoints),1),...
                        repmat(thisCondition,numel(timepoints),1),repmat(thisLockedTo,numel(timepoints),1),...
                        timepointRsa',timepointParams',r_rsaCorrs',repmat(thisModel,numel(timepoints),1),'VariableNames', varNames)
                        ];
                    clear timepoints timepointRsa timepointParams r_rsaCorrs
                end
                
            end; clear thisCondition
        end; clear thisLockedTo
    end; clear thisParameter
    clear count
    response = input('Do you want to delete the parallel pool? (y/n): ','s');
    if strcmpi(response,'y') || strcmpi(response,'yes')
        disp('deleting...');
        delete(gcp('nocreate')); clear workers;
        system(['rm -rf ' tempPath]);
        if exist(tempPath,'file'); warning('I couldnt delete the job directory :('); end
    else
        disp('not deleting...')
    end
    
    
    save(correlationsSaveName,'correlations')
    
end

disp('done')

% tmp1 = load([saveDir filesep 'model_correlations_conditions_01.mat']);
% tmp2 = load([saveDir filesep 'model_correlations_conditionsOnly_amplitude_bfs_null_-0.05,0.05.mat']);
% tmp1.correlations.bfs = tmp2.amplitudeBFs.bf2';
%% get the rsa correlation bfs
tmp = load([datadir filesep 'rsa_0.5,Inf.mat']);
parameters = unique([data(1).param.ParamName])';
parameterConds = unique([data(1).param.Condition])';
lockedTos = unique([data(1).rsaCorrs.LockedTo])';
conditions = unique([data(1).rsaCorrs.Condition])';
models = unique([data(1).rsaCorrs.Model])';
count = 0;
for thisLockedTo = lockedTos
    rsaLock = strrep(thisLockedTo,'Locked','');
    for thisCondition = conditions
        for thisModel = models
            
            tmpbf = tmp.rsa.([thisCondition{:} '_' rsaLock{:}]).(thisModel{:}).pbfs;
                            
            for thisParameter = parameters

                dataidx = ...
                    strcmp(correlations.ParameterCondition,thisCondition)...
                    & strcmp(correlations.Condition,thisCondition)...
                    & strcmp(correlations.ParameterName,thisParameter)...
                    & strcmp(correlations.LockedTo,thisLockedTo)...
                    & strcmp(correlations.Model,thisModel);
                correlations.RsaBfs(dataidx) = tmpbf;
                
            end
        end
    end
end

disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run the bayes analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('getting bayes factors')

nullInterval = '-0.1,0.1';
% nullInterval = '0.2,1';
bfRsaSavename = [saveDir filesep 'model_correlations_conditionsOnly_rsaCorrs_bfs_null_%s.mat'];
if exist(sprintf(bfRsaSavename,nullInterval),'file')
    warning('rsaCorrs bfs for this null interval exist');
    response = input('Do you want to overwrite? (y/n): ','s');
    if strcmpi(response,'y') || strcmpi(response,'yes')
        disp('removing file...');
        system(['rm -rf ' sprintf(bfRsaSavename,nullInterval)]);
        if exist(sprintf(bfRsaSavename,nullInterval),'file'); error('I couldnt delete it :('); end
        runRsa = 1;
    else
        disp('loading')
        load(sprintf(bfRsaSavename,nullInterval));
        runRsa = 0;
    end
else
    runRsa = 1;
end

if runRsa
    
    disp('calculating bfs')
    
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
        
        if runRsa
            % now run the rscript version of the bayes analysis
            abf1(i) = bayesfactor_R_wrapper_corr(correlations.ParamVals{i}', correlations.RsaCorrelations{i}','Rpath',RscriptPath,'returnindex',1,...
                'args',theseArgs,'verbose',0);
            abf2(i) = bayesfactor_R_wrapper_corr(correlations.ParamVals{i}', correlations.RsaCorrelations{i}','Rpath',RscriptPath,'returnindex',2,...
                'args',theseArgs,'verbose',0);
        end
        
    end
    if runRsa; rsaCorrsBFs.bf1 = abf1; rsaCorrsBFs.bf2 = abf2; save(sprintf(bfRsaSavename,nullInterval),'rsaCorrsBFs'); end
    response = input('Do you want to delete the parallel pool? (y/n): ','s');
    if strcmpi(response,'y') || strcmpi(response,'yes')
        disp('deleting...');
        delete(gcp('nocreate')); clear workers;
        system(['rm -rf ' tempPath]);
        if exist(tempPath,'file'); warning('I couldnt delete the job directory :('); end
    else
        disp('not deleting...')
    end
    clear runRsa
end

correlations.bf1 = rsaCorrsBFs.bf1';
correlations.bf2 = rsaCorrsBFs.bf2';

disp('done')

%%%%%%%%%%
%% plot %%
%%%%%%%%%%

disp(nullInterval)
plotInside = 0; % plot bf for inside the null interval (1) or outside (0)?
plotMovie = 0;

disp('plotting')


for thisParam = unique(correlations.ParameterName)'
    for thisCond = unique(correlations.Condition)'
        for thisLockedTo = unique(correlations.LockedTo)'
            paramCond = thisCond; % if you want to look cross cond, you can add this back as a for loop
            for thisModel = unique(correlations.Model)'
                
                dataidx = ...
                    strcmp(correlations.ParameterCondition,paramCond)...
                    & strcmp(correlations.ParameterName,thisParam)...
                    & strcmp(correlations.Condition,thisCond)...
                    & strcmp(correlations.LockedTo,thisLockedTo)...
                    & strcmp(correlations.Model,thisModel);
                
                theseParams = correlations.ParamVals(dataidx);
                
                theseData = correlations.RsaCorrelations(dataidx);
                if plotInside
                    thesebfs = correlations.bf1(dataidx);
                else
                    thesebfs = correlations.bf2(dataidx);
                end
                theseRs = correlations.R_Rsa(dataidx);
                allRs = correlations.R_Rsa;
                theseRsaBfs = correlations.RsaBfs(dataidx);
                theseRsaBfs(theseRsaBfs >= 10^2) = 10^2;
                theseRsaBfs(theseRsaBfs >= 10^1 & theseRsaBfs < 10^2) = 10^1;
                theseRsaBfs(theseRsaBfs >= 10^0.5 & theseRsaBfs < 10^1) = 10^0.5;
                theseRsaBfs(theseRsaBfs > 10^-0.5 & theseRsaBfs < 10^0.5) = NaN;
                theseRsaBfs(theseRsaBfs <= 10^-0.5 & theseRsaBfs > 10^-1) = 10^-0.5;
                theseRsaBfs(theseRsaBfs <= 10^-1 & theseRsaBfs > 10^-2) = 10^-1;
                theseRsaBfs(theseRsaBfs <= 10^-2) = 10^-2;
                
                
                slopeInvalidColour = [1, 0.6, 0.6];
                
                timepoints = correlations.Timepoint(dataidx);
                
                if contains(thisLockedTo,'response')
                    onset = 600; % to subtract from timepoints
                    xlims = [-600 200];
                elseif contains(thisLockedTo,'coherence')
                    onset = 500; % to subtract from timepoints
                    xlims = [-200 1500];
                end
                timepoints = timepoints-onset;
                
                if plotMovie
                    % little movie of the correlation
                    for i = 100+onset:numel(theseParams)
                        % scatterplot data
                        scatterHandles(i) = scatter(theseParams{i}, theseData{i});
                        hold on
                        % plot the slope line
                        p = polyfit(theseParams{i}, theseData{i}, 1);
                        slope = p(1);
                        yFit = polyval(p, xlim);
                        plot(xlim, yFit, 'r--');
                        % add some labels
                        title([thisCond ' ' thisParam ' ' thisLockedTo])
                        xlabel('Parameter Values')
                        ylabel('Amplitude')
                        annotation('textbox', [0.7, 0.8, 0.1, 0.1], 'String', ['Time: ' num2str(i-onset)], 'FitBoxToText', 'on');
                        if slope > 0;slope_str = 'Positive';elseif slope < 0;slope_str = 'Negative';else;slope_str = 'Zero';end
                        annotation('textbox', [0.7, 0.7, 0.1, 0.1], 'String', ['Slope: ' slope_str], 'FitBoxToText', 'on');
                        % set the axis limits
                        xlim([min(min([theseParams{:}])), max(max([theseParams{:}]))]);
                        ylim([min(min([theseData{:}])), max(max([theseData{:}]))]);
                        % do a quick pause for a short duration
                        pause(0.001);
                        % clear the figure
                        clf;
                    end; clear i
                    close all
                end
                
                figure;
                
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
                modelTitle = strrep(thisModel{:}, '_', '-');
                title(['Correlation: ' upper(thisCond{:}) ' ' thisLockedTo{:} ' (' thisParam{:} ' and ' modelTitle ')'])
                clear modelTitle
                
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
                hold on
                scatter(timepoints, theseRsaBfs, 20, 'k', 'filled')
                hold off
                ax = gca;
                set(ax,'YScale','log','XLim',xlims, ...
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
                cbh.Ticks = [-2, -1, -0.5, 0, 0.5, 1, 2];
                cbh.TickLabels=arrayfun(@(x) ['10^{' num2str(x) '}'], cbh.Ticks, 'UniformOutput', false);
                cbh.TickLabels(strcmp(cbh.TickLabels,'10^{0}')) = {'Inconclusive'};
                cbh.TickLabels(strcmp(cbh.TickLabels,'10^{0.5}') | strcmp(cbh.TickLabels,'10^{-0.5}')) = {'Moderate'};
                cbh.TickLabels(strcmp(cbh.TickLabels,'10^{1}') | strcmp(cbh.TickLabels,'10^{-1}')) = {'Strong'};
                        
                pause(1)
                print([saveDir filesep 'RSA_' thisParam{:} '_' paramCond{:} '_' thisCond{:} '_' thisLockedTo{:} '_' thisModel{:} '.png'],'-dpng')
                
            end; clear thisModel
        end; clear thisLockedTo
    end; clear thisCond
end; clear thisParam
close all

disp('done')

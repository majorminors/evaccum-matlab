%%%%%%%%%%%%
%% set up %%
%%%%%%%%%%%%

disp('setting up')

clear all

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';
datadir = fullfile(rootdir,'data','meg_pilot_4'); addpath(datadir);
toolsdir = fullfile(rootdir,'tools_analysis');
jobdir = fullfile(rootdir,'job_logging');
addpath(genpath(fullfile(toolsdir,'lib')))
ftDir = fullfile(toolsdir,'..','..','..','Toolboxes','fieldtrip');
addpath(ftDir); ft_defaults;

figDir = fullfile(datadir,'rsaFigs'); if ~exist(figDir,'dir');mkdir(figDir);end

toolbox = fullfile(rootdir,'..','..','Toolboxes','gramm'); addpath(toolbox); clear toolbox
toolbox = fullfile(rootdir,'..','..','Toolboxes','BFF_repo'); addpath(genpath(toolbox)); clear toolbox

disp('done setting up')

%%%%%%%%%%%%%%%
%% load data %%
%%%%%%%%%%%%%%%

disp('loading data')

% since our folders are all e.g. S01, S12 etc, let's load information about
% the directories that have that pattern
tmp = dir([datadir filesep 'S*']);
subjectFolders = tmp([tmp.isdir]); clear tmp
% and loop through them
count = 0;
for subjectNum = 1:numel(subjectFolders)
    % clear thisFile % can't do this in a parfor loop
    % full path to the input file
    thisFile = fullfile(subjectFolders(subjectNum).folder,subjectFolders(subjectNum).name,'rsa.mat');
    % skip this loop if we can't find one
    if ~exist(thisFile,'file'); continue; end
    % let's check the subject directory part of thisFile matches the
    % subject number
    pathParts = strsplit(thisFile, '/');
    index = find(contains(pathParts, 'S'));
    subjectCode = pathParts{index};
    if ~strcmp(subjectCode,subjectFolders(subjectNum).name); error('file doesnt match subject'); end
    
    fprintf('this is subject: %s...',subjectFolders(subjectNum).name)
    
    % load it
    fprintf('loading...')
    count = count+1;
    data{count} = load(thisFile);
    fprintf('loaded\n')
end; clear subjectNum index pathParts subjectCode subjectFolders thisFile

disp('done loading')

fprintf('you have %.0f subjects\n', count); clear count


% let's reorganise things a bit
subjects = numel(data);
rsa = struct();
for subject = 1:subjects % loop through subjects
    thisSubj = data{subject};
    thisSubj = thisSubj.rsa;
    conditions = fieldnames(thisSubj); % get the conditions (fieldnames)
    
    for condition = 1:numel(conditions) % loop though the conditions
        thisCond = conditions{condition};
        newCondName = strrep(thisCond,'LockedAverage',''); % we don't need this bit, so let's rename it
        
        % create the field if it doesn't exist
        if ~isfield(rsa, newCondName)
            rsa.(newCondName) = struct();
        end
        
        % now let's loop through the models
        models = fieldnames(thisSubj.(thisCond));
        for model = 1:numel(models)
            thisModel = models{model};
            % grab the values
            samples = thisSubj.(thisCond).(thisModel).samples;
            trans_samples = thisSubj.(thisCond).(thisModel).fisher_transformed_samples;
            
            % concatenate the values row-wise
            if isfield(rsa.(newCondName),thisModel)
                rsa.(newCondName).(thisModel).vals = cat(1, rsa.(newCondName).(thisModel).vals, samples);
                rsa.(newCondName).(thisModel).tvals = cat(1, rsa.(newCondName).(thisModel).tvals, trans_samples);
            else
                rsa.(newCondName).(thisModel).vals = samples;
                rsa.(newCondName).(thisModel).tvals = trans_samples;
            end
        end; clear model models thisModel samples trans_samples
    end; clear condition newCondName thisCond
end; clear subject conditions thisSubj subjects

clear data % don't need this anymore

disp('done prepping data')

%%%%%%%%%%%%%%%%%%%
%% calculate bfs %%
%%%%%%%%%%%%%%%%%%%

disp('collecting bayes factors')

nullInterval = '0.5,Inf';
insideNull = 1;
bfSaveName = [datadir filesep sprintf('rsa_%s.mat',nullInterval)];
saveFigs = 0;

if exist(bfSaveName,'file')
    warning('bfs for this null interval exist');
    response = input('Do you want to overwrite? (y/n): ','s');
    if strcmpi(response,'y') || strcmpi(response,'yes')
        disp('removing file...');
        system(['rm -rf ' bfSaveName]);
        if exist(bfSaveName,'file'); error('I couldnt delete it :('); end
        doBfs = 1;
    else
        disp('loading')
        load(bfSaveName);
        doBfs = 0;
    end
else
    doBfs = 1;
end

if doBfs
    conditions = fieldnames(rsa);
    numConditions = numel(conditions);
    conditionsBfs = cell(1,numConditions);
    
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        disp('no parallel pool is currently initialised---initialising');
        workers = numConditions;
        P = cbupool(workers);
        P.JobStorageLocation = jobdir;
        tempPath = fullfile(jobdir,'tmp');
        if ~exist(tempPath,'dir');mkdir(tempPath);end
        P.JobStorageLocation = tempPath;
        parpool(P,workers);
    else
        disp('parallel pool has already been initialized---skipping');
    end
    parfor condition = 1:numConditions
        thisCond = conditions{condition};
        models = fieldnames(rsa.(thisCond));
        numModels = numel(models);
        bfs = cell(1,numModels);
        tbfs = cell(1,numModels);
        
        % add the R module and get the path to Rscript
        [status, result] = system('module add R && which Rscript');
        if status == 0
            disp('R module added successfully');
            RscriptPath = strtrim(result);
            disp(['Rscript path: ' RscriptPath]);
        else
            error('Failed to add R and/or locate Rscript');
        end
        
        for model = 1:numModels
            thisModel = models{model};
            
            % now run the rscript version of the bayes analysis
            %   we can also get the bf for the complementary interval
            %   by specifying insideNull = 2 (default is 1)
            bfs{model} = bayesfactor_R_wrapper(rsa.(thisCond).(thisModel).vals','Rpath',RscriptPath,'returnindex',insideNull,...
                'args',['mu=0,rscale="medium",nullInterval=c(' nullInterval ')']);
            tbfs{model} = bayesfactor_R_wrapper(rsa.(thisCond).(thisModel).tvals','Rpath',RscriptPath,'returnindex',insideNull,...
                'args',['mu=0,rscale="medium",nullInterval=c(' nullInterval ')']);
        end
        conditionsBfs{condition} = bfs;
        conditionstBfs{condition} = tbfs;
    end; clear bfs tbfs
    
    % put the results back into the rsa struct
    for condition = 1:numConditions
        thisCond = conditions{condition};
        models = fieldnames(rsa.(thisCond));
        for model = 1:numel(models)
            thisModel = models{model};
            rsa.(thisCond).(thisModel).bfs = conditionsBfs{condition}{model};
            rsa.(thisCond).(thisModel).tbfs = conditionstBfs{condition}{model};
        end; clear models model thisModel
    end; clear conditions condition thisCond
    save(bfSaveName,'rsa')
    clear conditionsBfs conditionstBfs
    
    response = input('Do you want to delete the parallel pool? (y/n): ','s');
    if strcmpi(response,'y') || strcmpi(response,'yes')
        disp('deleting...');
        delete(gcp('nocreate')); clear workers;
        system(['rm -rf ' tempPath]);
        if exist(tempPath,'file'); warning('I couldnt delete the job directory :('); end
    else
        disp('not deleting...')
    end

end; clear doBfs

disp('done collecting bayes factors')

%%%%%%%%%%%%%%%%
%% plot them! %%
%%%%%%%%%%%%%%%%

disp('doing plots')

errTerm = 'error'; % std 'error' or std 'deviation'
plotTransformed = 0;

if plotTransformed
    whichVals = 'tvals';
else
    whichVals = 'vals';
end

% get ylims for the rsa
count = 0;
conditions = fieldnames(rsa);
for condition = 1:numel(conditions)
    thisCond = conditions{condition};
    models = fieldnames(rsa.(thisCond));
    for model = 1:numel(models)
        thisModel = models{model};
        count = count+1;
        vals = mean(rsa.(thisCond).(thisModel).(whichVals));
        switch errTerm
            case 'deviation'
                thisErr = std(rsa.(thisCond).(thisModel).(whichVals));
            case 'error'
                thisErr = std(rsa.(thisCond).(thisModel).(whichVals)) / sqrt(size(rsa.(thisCond).(thisModel).(whichVals), 1));
        end
        ymax(count) = max([max(thisErr) max(vals)]);
        ymin(count) = min([min(thisErr) min(vals)]);
    end; clear models model thisModel
end; clear conditions condition thisCond count
ymax = max(ymax)+max(ymax)/5;
ymin = min(ymin)-min(ymin)/5;

for lock = {'coherence' 'response'}
    
    figure;
    conditions = fieldnames(rsa);
    conditions = conditions(contains(conditions,lock));
    for condition = 1:numel(conditions)
        thisCond = conditions{condition};
                
        models = fieldnames(rsa.(thisCond));
        for model = 1:numel(models)
            thisModel = models{model};
            
            timepoints = 1:numel(rsa.(thisCond).(thisModel).bfs);
            bfs = rsa.(thisCond).(thisModel).bfs;
            vals = mean(rsa.(thisCond).(thisModel).(whichVals));
            std_dev = std(rsa.(thisCond).(thisModel).(whichVals));
            std_err = std(rsa.(thisCond).(thisModel).(whichVals)) / sqrt(size(rsa.(thisCond).(thisModel).(whichVals), 1));
            switch errTerm
                case 'deviation'
                    thisErr = std_dev;
                case 'error'
                    thisErr = std_err;
            end
            
            if contains(thisCond,'response')
                onset = 600; % to subtract from timepoints
                xlims = [-600 200];
                condName = thisCond;
            elseif contains(thisCond,'coherence')
                onset = 500; % to subtract from timepoints
                xlims = [-500 1500];
                condName = strrep(thisCond,'coherence','onset');
            end
            timepoints = timepoints-onset;
            
            % first the correlation
            subplot(2*numel(conditions),numel(models),(condition-1)*numel(models)*2 + model);
            corrCol = [0.4 0.8 0.6];
            fillCol = [1.0 0.7 0.7];
            x_fill = [timepoints, fliplr(timepoints)]; % x values for the fill
            y_upper = vals + thisErr; % upper bound of the fill
            y_lower = vals - thisErr; % lower bound of the fill
            fill(x_fill, [y_upper, fliplr(y_lower)], fillCol, 'FaceAlpha', 0.5, 'EdgeColor', 'none')
            hold on
            plot(timepoints,vals,'color',corrCol)
            hold off
            xlabel('Time (s)')
            ylabel(['Decoding (std' errTerm(1:3) ')'])
            ylim([ymin ymax])
            xlim([timepoints(1),timepoints(end)])
            title([strrep(condName,'_',' ') ' ' thisModel])
            clear corrCol fillCol
            
            % then the bayesfactors
            subplot(2*numel(conditions),numel(models),condition*numel(models)*2 - numel(models) + model);
            load('bayes_colourmap.mat'); % in BFF repo
            exponential_minmax=6;
            val_col_map = logspace(-exponential_minmax,exponential_minmax,size(colours,1));
            scatter_colours = zeros(length(timepoints), 3);  % preallocate for efficiency
            for t = 1:length(timepoints)
                [~,idx] = min(abs(val_col_map-bfs(t)));
                scatter_colours(t, :) = colours(idx,1:3);
            end
            scatter(timepoints, bfs, 30, scatter_colours, 'filled');
            line(get(gca,'XLim'),[1 1], 'Color', [0.7 0.7 0.7], 'LineStyle', '--')
            ax = gca;
            set(ax,'YScale','log','XLim',[timepoints(1),timepoints(end)], ...
                'YLim',[1*10^(-exponential_minmax) 1*10^exponential_minmax],'YTick',10.^(-exponential_minmax:2:exponential_minmax))
            xlabel('Time (s)')
            ylabel('BF (log scale)')
            colormap(colours)
            cbh = colorbar;
            caxis([-exponential_minmax,exponential_minmax])
            cbh.Units = 'normalized';
            cbh.Limits = [-exponential_minmax,exponential_minmax];
            cbh.Position(1) = 0.92;cbh.Position(3) = 0.01;cbh.Position(4) = ax.Position(4);cbh.Position(2) = ax.Position(2);
            cbh.Label.String = 'Bayes Factor';
            f = gcf; f.Position = [10 10 1600 1600];
            if exponential_minmax<=4
                cbh.Ticks = [-exponential_minmax, -1, -0.5, 0, 0.5, 1, exponential_minmax];
            else
                cbh.Ticks = [-exponential_minmax, -1, 0, 1, exponential_minmax];
            end
            cbh.TickLabels=arrayfun(@(x) ['10^{' num2str(x) '}'], cbh.Ticks, 'UniformOutput', false);
            cbh.TickLabels(strcmp(cbh.TickLabels,'10^{0}')) = {'Inconclusive'};
            cbh.TickLabels(strcmp(cbh.TickLabels,'10^{0.5}') | strcmp(cbh.TickLabels,'10^{-0.5}')) = {'Moderate'};
            cbh.TickLabels(strcmp(cbh.TickLabels,'10^{1}') | strcmp(cbh.TickLabels,'10^{-1}')) = {'Strong'};
            
        end; clear models model thisModel
    end; clear conditions condition thisCond
    
    if saveFigs
        tmp = split(condName,'_');
        print([figDir filesep tmp{2} '_rsa.png'],'-dpng')
        clear tmp
    end

end; clear lock

disp('done doing plots')


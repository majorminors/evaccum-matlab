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
            part_samples = thisSubj.(thisCond).(thisModel).partial_samples;
            part_trans_samples = thisSubj.(thisCond).(thisModel).fisher_transformed_partial_samples;
            timepoints = thisSubj.(thisCond).(thisModel).a.fdim.values{1};
            
            % concatenate the values row-wise
            if isfield(rsa.(newCondName),thisModel)
                rsa.(newCondName).(thisModel).vals = cat(1, rsa.(newCondName).(thisModel).vals, samples);
                rsa.(newCondName).(thisModel).tvals = cat(1, rsa.(newCondName).(thisModel).tvals, trans_samples);
                rsa.(newCondName).(thisModel).pvals = cat(1, rsa.(newCondName).(thisModel).vals, part_samples);
                rsa.(newCondName).(thisModel).ptvals = cat(1, rsa.(newCondName).(thisModel).tvals, part_trans_samples);
                if any(rsa.(newCondName).(thisModel).timepoints ~= timepoints); error('you have different timepoints'); end
            else
                rsa.(newCondName).(thisModel).vals = samples;
                rsa.(newCondName).(thisModel).tvals = trans_samples;
                rsa.(newCondName).(thisModel).pvals = part_samples;
                rsa.(newCondName).(thisModel).ptvals = part_trans_samples;
                rsa.(newCondName).(thisModel).timepoints = timepoints;
            end
            clear timepoints
        end; clear model models thisModel samples trans_samples part_samples part_trans_samples
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
        
        % initialise these so the parfor loop can fill them up
        bfs = cell(1,numModels);
        tbfs = cell(1,numModels);
        pbfs = cell(1,numModels);
        ptbfs = cell(1,numModels);
        
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
            pbfs{model} = bayesfactor_R_wrapper(rsa.(thisCond).(thisModel).pvals','Rpath',RscriptPath,'returnindex',insideNull,...
                'args',['mu=0,rscale="medium",nullInterval=c(' nullInterval ')']);
            ptbfs{model} = bayesfactor_R_wrapper(rsa.(thisCond).(thisModel).ptvals','Rpath',RscriptPath,'returnindex',insideNull,...
                'args',['mu=0,rscale="medium",nullInterval=c(' nullInterval ')']);
        end
        conditionsBfs{condition} = bfs;
        conditionstBfs{condition} = tbfs;
        conditionspBfs{condition} = pbfs;
        conditionsptBfs{condition} = ptbfs;
    end; clear bfs tbfs
    
    % put the results back into the rsa struct
    for condition = 1:numConditions
        thisCond = conditions{condition};
        models = fieldnames(rsa.(thisCond));
        for model = 1:numel(models)
            thisModel = models{model};
            rsa.(thisCond).(thisModel).bfs = conditionsBfs{condition}{model};
            rsa.(thisCond).(thisModel).tbfs = conditionstBfs{condition}{model};
            rsa.(thisCond).(thisModel).pbfs = conditionspBfs{condition}{model};
            rsa.(thisCond).(thisModel).ptbfs = conditionsptBfs{condition}{model};
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% analyse rsa differences %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nullInterval = '0.5,Inf';
insideNull = 1;
bfDiffSaveName = [datadir filesep sprintf('rsa_diffs_%s.mat',nullInterval)];
saveFigs = 1;

if exist(bfDiffSaveName,'file')
    warning('bfs for this null interval exist');
    response = input('Do you want to overwrite? (y/n): ','s');
    if strcmpi(response,'y') || strcmpi(response,'yes')
        disp('removing file...');
        system(['rm -rf ' bfDiffSaveName]);
        if exist(bfDiffSaveName,'file'); error('I couldnt delete it :('); end
        doBfDiffs = 1;
    else
        disp('loading')
        load(bfDiffSaveName);
        doBfDiffs = 0;
    end
else
    doBfDiffs = 1;
end

models = {'stim' 'decbdry' 'dec_simple' 'resp'};
modelNames = {'Motion Direction' 'Decision Boundary' 'Motion Classification' 'Response'};
coh = {'ec_' 'hc_'};
cat = {'er_' 'hr_'};
conds = {'ecer_' 'echr_' 'hcer_' 'hchr_'};

if doBfDiffs
    
    [status, result] = system('module add R && which Rscript');
    if status == 0
        disp('R module added successfully');
        RscriptPath = strtrim(result);
        disp(['Rscript path: ' RscriptPath]);
    else
        error('Failed to add R and/or locate Rscript');
    end
    
    for model = models
        fprintf('model %s...',model{:});
        
        for lockedTo = {'response' 'coherence'}
            fprintf('locked to %s\n',lockedTo{:});
            
            if ~exist('bfs','var')
                diffBfs = struct();
            end
            if ~isfield(diffBfs,'coh')
                diffBfs.coh = [];
            end
            if ~isfield(diffBfs,'cat')
                diffBfs.cat = [];
            end
            if ~isfield(diffBfs,'int')
                diffBfs.int = [];
            end
            
            diffBfs.coh = [diffBfs.coh {model{:}; lockedTo{:}; bayesfactor_R_wrapper(rsa.([coh{1} lockedTo{:}]).(model{:}).vals'-rsa.([coh{2} lockedTo{:}]).(model{:}).vals',...
                'Rpath',RscriptPath,'returnindex',insideNull,...
                'args',['mu=0,rscale="medium",nullInterval=c(' nullInterval ')'])}];
            diffBfs.cat = [diffBfs.cat {model{:}; lockedTo{:}; bayesfactor_R_wrapper(rsa.([cat{1} lockedTo{:}]).(model{:}).vals'-rsa.([cat{2} lockedTo{:}]).(model{:}).vals',...
                'Rpath',RscriptPath,'returnindex',insideNull,...
                'args',['mu=0,rscale="medium",nullInterval=c(' nullInterval ')'])}];
            diffBfs.int = [diffBfs.int {model{:}; lockedTo{:}; bayesfactor_R_wrapper(...
                (rsa.([conds{1} lockedTo{:}]).(model{:}).vals'-rsa.([conds{2} lockedTo{:}]).(model{:}).vals')-...
                (rsa.([conds{3} lockedTo{:}]).(model{:}).vals'-rsa.([conds{4} lockedTo{:}]).(model{:}).vals'),...
                'Rpath',RscriptPath,'returnindex',insideNull,...
                'args',['mu=0,rscale="medium",nullInterval=c(' nullInterval ')'])}];
            
        end; clear lockedTo
        
    end; clear model
    
    save(bfDiffSaveName,'diffBfs')
end; clear doBfDiffs

types = {'coh', 'cat'};

for h = {'response' 'coherence'}
    figure;
    for i = 1:numel(types)
        type = types{i};
        for j = 1:numel(models)
            model_name = models{j};
            theseBfs = diffBfs.(type){3,find(strcmp(model_name, diffBfs.(type)(1,:)) & strcmp(h, diffBfs.(type)(2,:)))};
            subplot(numel(types), numel(models), (i-1)*numel(models)+j);
            
            timepoints = rsa.(['ecer_' h{:}]).(model_name).timepoints';
            
            if contains(h,'response')
                xlims = [-0.6 0.2];
            elseif contains(h,'coherence')
                xlims = [-0.2 1.5];
            end
            
            load('bayes_colourmap.mat'); % in BFF repo
            exponential_minmax=6;
            val_col_map = logspace(-exponential_minmax,exponential_minmax,size(colours,1));
            scatter_colours = zeros(length(timepoints), 3);  % preallocate for efficiency
            for t = 1:length(timepoints)
                [~,idx] = min(abs(val_col_map-theseBfs(t)));
                scatter_colours(t, :) = colours(idx,1:3);
            end
            scatter(timepoints, theseBfs, 30, scatter_colours, 'filled');
            line(get(gca,'XLim'),[1 1], 'Color', [0.7 0.7 0.7], 'LineStyle', '--') % plot a line along inconclusive evidence
            ax = gca;
            set(ax,'YScale','log','XLim',xlims, ...
                'YLim',[1*10^(-exponential_minmax) 1*10^exponential_minmax],'YTick',10.^(-exponential_minmax:2:exponential_minmax))
            xlabel('Time (s)')
            ylabel('BF (log scale)')
            if strcmp(type,'cat')
                tmp = 'Cat Diff';
            elseif strcmp(type,'coh')
                tmp = 'Coh Diff';
            elseif strcmp(type,'int')
                tmp = 'Cat Interaction';
            end
            title([tmp ' in ' modelNames{j}]);
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
            
            
        end; clear j
    end; clear i
    
    tmp = h{:}; tmp(1) = upper(tmp(1));
    suptitle([tmp ' Locked'])
    clear tmp
    
    if saveFigs
        print([figDir filesep h{:} '_locked_rsa_diffs.png'],'-dpng')
    end
end; clear h


%%%%%%%%%%%%%%%%
%% plot them! %%
%%%%%%%%%%%%%%%%

disp('doing plots')
saveFigs = 1;
plotConds = {'stim' 'decbdry' 'dec_simple' 'resp'};
plotNames = {'Motion Direction' 'Decision Boundary' 'Motion Classification' 'Response'};

errTerm = 'error'; % std 'error' or std 'deviation'
plotTransformed = 1;
plotPartials = 0;



whichVals = 'vals';
whichBfs = 'bfs';
if plotTransformed
    whichVals = ['t' whichVals];
    whichBfs = ['t' whichBfs];
end
if plotPartials
    whichVals = ['p' whichVals];
    whichBfs = ['p' whichBfs];
end

plotGroups = {{'ecer_' 'echr_' 'hcer_' 'hchr_'} {'ec_' 'hc_' 'er_' 'hr_'}};
groupNames = {'conditions' 'manipulations'};
for thisGroup = 1:numel(plotGroups)
    tmp = fieldnames(rsa);
    conditions = tmp(startsWith(fieldnames(rsa),plotGroups{thisGroup})); clear tmp
    
    % get ylims for the rsa
    count = 0;
    for condition = 1:numel(conditions)
        thisCond = conditions{condition};
        % models = fieldnames(rsa.(thisCond));
        models = plotConds;
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
    end; clear condition thisCond count
    ymax = max(ymax)+max(ymax)/5;
    ymin = min(ymin)-min(ymin)/5;
    
    for lock = {'coherence' 'response'}
        
        figure;
        lockConds = conditions(contains(conditions,lock));
        for condition = 1:numel(lockConds)
            thisCond = lockConds{condition};
            
            % models = fieldnames(rsa.(thisCond));
            models = plotConds;
            for model = 1:numel(models)
                thisModel = models{model};
                
                timepoints = rsa.(thisCond).(thisModel).timepoints';
                bfs = rsa.(thisCond).(thisModel).(whichBfs);
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
                    xlims = [-0.6 0.2];
                    condName = thisCond;
                elseif contains(thisCond,'coherence')
                    xlims = [-0.2 1.5];
                    condName = strrep(thisCond,'coherence','onset');
                end
                tmp = abs(xlims-timepoints); % get absolute difference between desired xlims and timepoints
                [~,tmp] = min(tmp); % get the indices of the minimum difference
                xlims = timepoints(tmp)'; % use the timepoints closest to the desired xlims (by minimum absolute difference)
                
                % first the correlation
                subplot(2*numel(lockConds),numel(models),(condition-1)*numel(models)*2 + model);
                corrCol = [0.4 0.8 0.6];
                fillCol = [1.0 0.7 0.7];
                x_fill = [timepoints', fliplr(timepoints')]; % x values for the fill
                y_upper = vals + thisErr; % upper bound of the fill
                y_lower = vals - thisErr; % lower bound of the fill
                fill(x_fill, [y_upper, fliplr(y_lower)], fillCol, 'FaceAlpha', 0.5, 'EdgeColor', 'none') % the fill
                hold on
                plot(timepoints',vals,'color',corrCol) % the actual rsa correlation
                line(get(gca,'XLim'), [0 0], 'Color', [0.1 0.1 0.1], 'LineStyle', '-') % plot a line along zero on the y
                xlabel('Time (s)')
                ylabel(['Decoding (std' errTerm(1:3) ')'])
                ylim([ymin ymax])
                line([0 0], get(gca,'YLim'), 'Color', [0.7 0.7 0.7], 'LineStyle', '--') % plot a line to mark the timepoint of interest
                hold off
                xlim(xlims)
                title([regexprep(regexprep(strrep(condName,'_',' '), '(?<!\S)(\S)(\S*)', '${upper($1)}${lower($2)}'), '(\S+\s)(\S+)', '$1$2-Locked')...
                    ' ' plotNames{strcmp(plotConds,thisModel)}])
                clear corrCol fillCol
                
                % then the bayesfactors
                subplot(2*numel(lockConds),numel(models),condition*numel(models)*2 - numel(models) + model);
                load('bayes_colourmap.mat'); % in BFF repo
                exponential_minmax=6;
                val_col_map = logspace(-exponential_minmax,exponential_minmax,size(colours,1));
                scatter_colours = zeros(length(timepoints), 3);  % preallocate for efficiency
                for t = 1:length(timepoints)
                    [~,idx] = min(abs(val_col_map-bfs(t)));
                    scatter_colours(t, :) = colours(idx,1:3);
                end
                scatter(timepoints, bfs, 30, scatter_colours, 'filled');
                line(get(gca,'XLim'),[1 1], 'Color', [0.7 0.7 0.7], 'LineStyle', '--') % plot a line along inconclusive evidence
                ax = gca;
                set(ax,'YScale','log','XLim',xlims, ...
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
        end; clear condition thisCond lockConds
        
        if saveFigs
            tmp = split(condName,'_');
            print([figDir filesep groupNames{thisGroup} '_' tmp{2} '_rsa.png'],'-dpng')
            clear tmp
        end
        
    end; clear lock
end; clear thisGroup conditions

disp('done doing plots')


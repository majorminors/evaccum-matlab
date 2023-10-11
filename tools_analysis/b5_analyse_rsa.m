%%%%%%%%%%%%
%% set up %%
%%%%%%%%%%%%

clear all

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';
datadir = fullfile(rootdir,'data','meg_pilot_4'); addpath(datadir);
toolsdir = fullfile(rootdir,'tools_analysis');
jobdir = fullfile(rootdir,'job_logging');
addpath(genpath(fullfile(toolsdir,'lib')))
ftDir = fullfile(toolsdir,'..','..','..','Toolboxes','fieldtrip');
addpath(ftDir); ft_defaults;

toolbox = fullfile(rootdir,'..','..','Toolboxes','gramm'); addpath(toolbox); clear toolbox
toolbox = fullfile(rootdir,'..','..','Toolboxes','BFF_repo'); addpath(genpath(toolbox)); clear toolbox

%%%%%%%%%%%%%%%
%% load data %%
%%%%%%%%%%%%%%%

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
    
    fprintf('this is subject: %s\n',subjectFolders(subjectNum).name)
    
    % load it
    disp('loading')
    count = count+1;
    data{count} = load(thisFile);
    disp('loaded')
end; clear subjectNum index pathParts subjectCode subjectFolders thisFile count

disp('done loading')

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
            fieldValue = thisSubj.(thisCond).(thisModel).samples;
            
            % concatenate the values row-wise
            if isfield(rsa.(newCondName),thisModel)
                rsa.(newCondName).(thisModel).vals = cat(1, rsa.(newCondName).(thisModel).vals, fieldValue);
            else
                rsa.(newCondName).(thisModel).vals = fieldValue;
            end
        end; clear model models thisModel fieldValue
    end; clear condition newCondName thisCond
end; clear subject conditions thisSubj subjects

clear data % don't need this anymore

%%%%%%%%%%%%%%%%%%%
%% calculate bfs %%
%%%%%%%%%%%%%%%%%%%

nullInterval = '0.5,Inf';
bfSaveName = [datadir filesep sprintf('rsa_%s.mat',nullInterval)];

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
    for condition = 1:numel(conditions)
        thisCond = conditions{condition};
        
        models = fieldnames(rsa.(thisCond));
        for model = 1:numel(models)
            thisModel = models{model};
            
            % add the R module and get the path to Rscript
            [status, result] = system('module add R && which Rscript');
            if status == 0
                disp('R module added successfully');
                RscriptPath = strtrim(result);
                disp(['Rscript path: ' RscriptPath]);
            else
                error('Failed to add R and/or locate Rscript');
            end
            
            % now run the rscript version of the bayes analysis
            %   we can also get the bf for the complementary interval
            %   by specifying complementary = 2. Let's set a default:
            if ~exist('complementary','var'); complementary = 1; end
            rsa.(thisCond).(thisModel).bfs = bayesfactor_R_wrapper(rsa.(thisCond).(thisModel).vals','Rpath',RscriptPath,'returnindex',complementary,...
                'args',['mu=0,rscale="medium",nullInterval=c(' nullInterval ')']);
            
        end; clear models model thisModel
    end; clear conditions condition thisCond
    
    save(bfSaveName,'rsa')
end; clear doBfs

%%%%%%%%%%%%%%%%
%% plot them! %%
%%%%%%%%%%%%%%%%

errTerm = 'error'; % std 'error' or std 'deviation'

% get ylims for the rsa
count = 0;
conditions = fieldnames(rsa);
for condition = 1:numel(conditions)
    thisCond = conditions{condition};
    models = fieldnames(rsa.(thisCond));
    for model = 1:numel(models)
        thisModel = models{model};
        count = count+1;
        vals = mean(rsa.(thisCond).(thisModel).vals);
        switch errTerm
            case 'deviation'
                thisErr = std(rsa.(thisCond).(thisModel).vals);
            case 'error'
                thisErr = std(rsa.(thisCond).(thisModel).vals) / sqrt(size(rsa.(thisCond).(thisModel).vals, 1));
        end
        ymax(count) = max([max(thisErr) max(vals)]);
        ymin(count) = min([min(thisErr) min(vals)]);
    end; clear models model thisModel
end; clear conditions condition thisCond count
ymax = max(ymax);
ymin = min(ymin);

for lock = {'coherence' 'response'}
    
    figure;
    subplotCounter = 0;
    conditions = fieldnames(rsa);
    conditions = conditions(contains(conditions,lock));
    for condition = 1:numel(conditions)
        thisCond = conditions{condition};
                
        models = fieldnames(rsa.(thisCond));
        for model = 1:numel(models)
            thisModel = models{model};
            subplotCounter = subplotCounter+1;
            
            timepoints = 1:numel(rsa.(thisCond).(thisModel).bfs);
            bfs = rsa.(thisCond).(thisModel).bfs;
            vals = mean(rsa.(thisCond).(thisModel).vals);
            std_dev = std(rsa.(thisCond).(thisModel).vals);
            std_err = std(rsa.(thisCond).(thisModel).vals) / sqrt(size(rsa.(thisCond).(thisModel).vals, 1));
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
            ylabel('RSA Correlation')
            ylim([ymin ymax])
            xlim([timepoints(1),timepoints(end)])
            title([strrep(condName,'_',' ') ' ' thisModel])
            clear corrCol fillCol
            
            % then the bayesfactors
            subplot(2*numel(conditions),numel(models),condition*numel(models)*2 - numel(models) + model);
            load('bayes_colourmap.mat'); % in BFF repo
            exponential_minmax=4;
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
                'YLim',[1e-4 1e4],'YTick',10.^(-4:1:4))
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
            %         cbh.Ticks = [-6, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 6];
            cbh.TickLabels=arrayfun(@(x) ['10^{' num2str(x) '}'], cbh.Ticks, 'UniformOutput', false);
            cbh.TickLabels(strcmp(cbh.TickLabels,'10^{0}')) = {'Inconclusive'};
            cbh.TickLabels(strcmp(cbh.TickLabels,'10^{0.5}') | strcmp(cbh.TickLabels,'10^{-0.5}')) = {'Moderate'};
            cbh.TickLabels(strcmp(cbh.TickLabels,'10^{1}') | strcmp(cbh.TickLabels,'10^{-1}')) = {'Strong'};
            
        end; clear models model thisModel
    end; clear conditions condition thisCond
end; clear lock

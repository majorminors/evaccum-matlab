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

% gather the data up
for subjectidx = 1:size(validSubjectParams,1)
    %r = b4_collateRs(r,timepoint,{'Decision Bdry (coh)'},{'Onset Locked Coherence Amplitude'},b_coh_diff, coh_onset(timepoint, :));
    data.params(subjectidx) = {'Decision Bdry' 'Coh' validSubjectParams.b_ec(subjectidx)-validSubjectParams.b_hc(subjectidx)};
    data.params(subjectidx) = {'Decision Bdry' 'Cat' validSubjectParams.b_er(subjectidx)-validSubjectParams.b_hr(subjectidx)};
    
    b_coh_diff(subjectidx) = ;
    b_cat_diff(subjectidx) = validSubjectParams.b_er(subjectidx)-validSubjectParams.b_hr(subjectidx);
    v_coh_diff(subjectidx) = validSubjectParams.v1_ec(subjectidx)-validSubjectParams.v1_hc(subjectidx);
    v_cat_diff(subjectidx) = validSubjectParams.v1_er(subjectidx)-validSubjectParams.v1_hr(subjectidx);

    coh_onset(:,subjectidx) = getRawAmplitude(cOnsDiffAll{subjectidx},CPP);
    coh_response(:,subjectidx) = getRawAmplitude(cRespDiffAll{subjectidx},CPP);
    cat_onset(:,subjectidx) = getRawAmplitude(rOnsDiffAll{subjectidx},CPP);
    cat_response(:,subjectidx) = getRawAmplitude(rRespDiffAll{subjectidx},CPP);

    % get the change in the local slope
    for t = 1:length(coh_onset(:,subjectidx))
        if t<11
            coh_onset_slope(t,subjectidx) = coh_onset(t,1)-mean(coh_onset(:,1));
            cat_onset_slope(t,subjectidx) = cat_onset(t,1)-mean(cat_onset(:,1));
        else
            coh_onset_slope(t,subjectidx) = coh_onset(t,1)-coh_onset(t-10,1);
            cat_onset_slope(t,subjectidx) = cat_onset(t,1)-cat_onset(t-10,1);
        end
    end; clear t
    for t = 1:length(coh_response(:,subjectidx))
        if t<11
            coh_response_slope(t,subjectidx) = coh_response(t,1)-mean(coh_response(:,1));
            cat_response_slope(t,subjectidx) = cat_response(t,1)-mean(cat_response(:,1));
        else
            coh_response_slope(t,subjectidx) = coh_response(t,1)-coh_response(t-10,1);
            cat_response_slope(t,subjectidx) = cat_response(t,1)-cat_response(t-10,1);
        end
    end; clear t

end; clear subjectidx


% now correlate

r = struct; 

for timepoint = 1:size(coh_onset,1)

    % correlation between decision boundary for e/h coherence and evoked response for e/h coherence onset
    r = b4_collateRs(r,timepoint,{'Decision Bdry (coh)'},{'Onset Locked Coherence Amplitude'},b_coh_diff, coh_onset(timepoint, :));
    % correlation between decision boundary for e/h coherence and evoked response for e/h categorisation onset
    r = b4_collateRs(r,timepoint,{'Decision Bdry (coh)'},{'Onset Locked Categorisation Amplitude'},b_coh_diff,cat_onset(timepoint, :));
    % correlation between decision boundary for e/h categorisation and evoked response for e/h coherence onset
    r = b4_collateRs(r,timepoint,{'Decision Bdry (cat)'},{'Onset Locked Coherence Amplitude'},b_cat_diff,coh_onset(timepoint, :));
    % correlation between decision boundary for e/h categorisation and evoked response for e/h categorisation onset
    r = b4_collateRs(r,timepoint,{'Decision Bdry (cat)'},{'Onset Locked Categorisation Amplitude'},b_cat_diff,cat_onset(timepoint, :));
    % correlation between drift rate for e/h coherence and evoked response for e/h coherence onset
    r = b4_collateRs(r,timepoint,{'Drift Rate (coh)'},{'Onset Locked Coherence Amplitude'},v_coh_diff,coh_onset(timepoint, :));
    % correlation between drift rate for e/h coherence and evoked response for e/h categorisation onset
    r = b4_collateRs(r,timepoint,{'Drift Rate (coh)'},{'Onset Locked Categorisation Amplitude'},v_coh_diff,cat_onset(timepoint, :));
    % correlation between drift rate for e/h categorisation and evoked response for e/h coherence onset
    r = b4_collateRs(r,timepoint,{'Drift Rate (cat)'},{'Onset Locked Coherence Amplitude'},v_cat_diff,coh_onset(timepoint, :));
    % correlation between drift rate for e/h categorisation and evoked response for e/h categorisation onset
    r = b4_collateRs(r,timepoint,{'Drift Rate (cat)'},{'Onset Locked Categorisation Amplitude'},v_cat_diff,cat_onset(timepoint, :));
    
end

for timepoint = 1:size(coh_response,1)

    % correlation between decision boundary for e/h coherence and evoked response for e/h coherence response
    r = b4_collateRs(r,timepoint,{'Decision Bdry (coh)'},{'Response Locked Coherence Amplitude'},b_coh_diff,coh_response(timepoint, :));
    % correlation between decision boundary for e/h coherence and evoked response for e/h categorisation response
    r = b4_collateRs(r,timepoint,{'Decision Bdry (coh)'},{'Response Locked Categorisation Amplitude'},b_coh_diff,cat_response(timepoint, :));
    % correlation between decision boundary for e/h categorisation and evoked response for e/h coherence response
    r = b4_collateRs(r,timepoint,{'Decision Bdry (cat)'},{'Response Locked Coherence Amplitude'},b_cat_diff,coh_response(timepoint, :));
    % correlation between decision boundary for e/h categorisation and evoked response for e/h categorisation response
    r = b4_collateRs(r,timepoint,{'Decision Bdry (cat)'},{'Response Locked Categorisation Amplitude'},b_cat_diff,cat_response(timepoint, :));
    % correlation between drift rate for e/h coherence and evoked response for e/h coherence response
    r = b4_collateRs(r,timepoint,{'Drift Rate (coh)'},{'Response Locked Coherence Amplitude'},v_coh_diff,coh_response(timepoint, :));
    % correlation between drift rate for e/h coherence and evoked response for e/h categorisation response
    r = b4_collateRs(r,timepoint,{'Drift Rate (coh)'},{'Response Locked Categorisation Amplitude'},v_coh_diff,cat_response(timepoint, :));
    % correlation between drift rate for e/h categorisation and evoked response for e/h coherence response
    r = b4_collateRs(r,timepoint,{'Drift Rate (cat)'},{'Response Locked Coherence Amplitude'},v_cat_diff,coh_response(timepoint, :));
    % correlation between drift rate for e/h categorisation and evoked response for e/h categorisation response
    r = b4_collateRs(r,timepoint,{'Drift Rate (cat)'},{'Response Locked Categorisation Amplitude'},v_cat_diff,cat_response(timepoint, :));
    
end

% workers = 100;
% P=cbupool(workers);
% parpool(P,workers);
% nullInterval = '0.2,1';
% parfor i = 1:numel(r.timepoint)
%     fprintf('%.0f of %.0f\n',i,numel(r.timepoint))
%     
%     % add the R module and get the path to Rscript
%     [status, result] = system('module add R && which Rscript');
%     if status == 0
% %         disp('R module added successfully');
%         RscriptPath = strtrim(result);
% %         disp(['Rscript path: ' RscriptPath]);
%     else
%         error('Failed to add R and/or locate Rscript');
%     end
%     
%     theseArgs = ['rscale="medium",nullInterval=c(' nullInterval ')'];
%     
%     % now run the rscript version of the bayes analysis
%     bf1(i) = bayesfactor_R_wrapper_corr(r.paramVals{i}, r.rawAmplitude{i},'Rpath',RscriptPath,'returnindex',1,...
%         'args',theseArgs,'verbose',0);
%     bf2(i) = bayesfactor_R_wrapper_corr(r.paramVals{i}, r.rawAmplitude{i},'Rpath',RscriptPath,'returnindex',2,...
%         'args',theseArgs,'verbose',0);
%     
% end
% % delete(gcp('nocreate')); clear workers
bfSavename = [saveDir filesep 'model_correlation_bfs_null_%s.mat'];
% save(sprintf(bfSavename,nullInterval),'bf1','bf2')

load(sprintf(bfSavename,'-0.2,0.2'))

r.bf1 = bf1;
r.bf2 = bf2;

save([saveDir filesep 'model_correlations.mat'],'r')

subset = strcmp(r.dataType,'coh_ons') | strcmp(r.dataType,'cat_ons')
g = gramm('x',r.timepoint,'y',r.r,'color',r.param,'subset',subset)
g.facet_grid(r.dataType,[])
% Plot raw data as points
g.geom_line()
% % Plot linear fits of the data with associated confidence intervals
% g.stat_glm()
% Set appropriate names for legends
g.set_names('row','','x','Timepoint','y','Correlation (r)','color','Parameter')
% %Set figure title
% g.set_title('Fuel economy of new cars between 1970 and 1982')
% Do the actual drawing
g.draw()


subset = strcmp(r.dataType,'coh_resp') | strcmp(r.dataType,'cat_resp')
g = gramm('x',r.timepoint,'y',r.r,'color',r.param,'subset',subset)
g.facet_grid(r.dataType,[])
% Plot raw data as points
g.geom_line()
% % Plot linear fits of the data with associated confidence intervals
% g.stat_glm()
% Set appropriate names for legends
g.set_names('row','','x','Timepoint','y','Correlation (r)','color','Parameter')
% %Set figure title
% g.set_title('Fuel economy of new cars between 1970 and 1982')
% Do the actual drawing
g.draw()


for dataType = unique(r.dataType)
    for param = unique(r.param)
        figure;

        dataidx = strcmp(r.dataType,dataType) & strcmp(r.param,param);
%         thesebfs = r.bf1(dataidx);
        thesebfs = r.bf2(dataidx);
        rs = r.r(dataidx);
        timepoints = r.timepoint(dataidx);
        if contains(dataType,'Response')
            onset = 600; % to subtract from timepoints
            xlims = [-600 200];
        elseif contains(dataType,'Onset')
            onset = 500; % to subtract from timepoints
            xlims = [-500 1500];
        end
        timepoints = timepoints-onset;
        
        % plot correlations
        subplot(2,1,1)
        color_map = summer(length(rs));
        color_map = color_map(:, [2, 1, 3]);  % Rearrange color map channels
        abs_rs = abs(rs);
        scatter_colours = 1 - abs_rs / max(abs_rs);  % Scaling the color values to make darker as value moves away from 0
        scatter(timepoints, rs, [], scatter_colours, 'filled');
        colormap(color_map);
%         cbh = colorbar;
        xlim(xlims)
        ylim([min(r.r)-0.1 max(r.r)+0.1])
        line(get(gca,'XLim'), [0 0], 'Color', [0.5 0.5 0.5], 'LineStyle', '-')
        line([0 0], get(gca,'YLim'), 'Color', [0.1 0.1 0.1], 'LineStyle', '--')
        xlabel('Timepoints');
        ylabel('Correlations (r)');
        title(['Correlation between ' param{:} ' and ' dataType])

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

        print([saveDir filesep param{:} '_' dataType{:} '.png'],'-dpng')
        
    end
end
% close all

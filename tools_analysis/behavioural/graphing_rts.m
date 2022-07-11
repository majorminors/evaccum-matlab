        %% Graph rts as histogram

% produces:
% d.allsubjs:
%   participant id | rts | accuracy | coherence condition (1 = easy, 2 = hard) | rule condition (1 = easy, 2 = hard)
% d.subjects: structure with above information separated by both subject and variable
% histograms for each subject and all subjects in following layout:
%   easy coh, easy rule | eash coh, hard rule
%   -----------------------------------------
%   hard coh, easy rule | hard coh, hard rule

% currently only saving rts - turn p.save_figs on or off to do this
% you can save variables using path 'p.save_file'

%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);
p = struct(); % keep some of our parameters tidy - these are things you might want to edit
d = struct(); % set up a structure for the data info - these are outputs we want to save
t = struct(); % set up a structure for temp data - these are things that are liable to change throughout the script

% set up variables
rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab'; % % root directory - used to inform directory mappings
datadir = fullfile(rootdir,'data','meg_pilot_4');
%p.subjectrange = [1:4,7:9,12,20:23,26];
p.subjectrange = [1:14,19:23,25,26,29:33,35:37];
p.save_figs = 0;
p.datafilepattern = '*_EvAccum.mat';
p.vars = {'d','block'}; % which variables do you need
p.savefilename = 'rt_dists';

% directory mapping
addpath(genpath(fullfile(rootdir, 'tools_analysis'))); % add tools folder to path (don't think we need this, but in case)
behavdatadir = fullfile(datadir,'behavioural'); % add matlab behavioural data
outdir = fullfile(datadir,'descriptives'); % find or make directory to output lba fit results
if ~exist(outdir,'dir')
    mkdir(outdir);
end
save_file = fullfile(outdir,p.savefilename);

%% get the data
d.fileinfo = dir(fullfile(behavdatadir, p.datafilepattern)); % find all the datafiles and get their info
d.allsubjs = {}; % create this so we can work with it
for i = p.subjectrange%length(d.fileinfo) % loop through each
    t.path = fullfile(behavdatadir, d.fileinfo(i).name); % get the full path to the file
    fprintf(1, 'working with %s\n', t.path); % print that so you can check
    
    t.alldata = load(t.path,p.vars{:}); % load in the data
    
    t.id = t.alldata.d.participant_id;
    t.blocks = t.alldata.block; % how many blocks?
    
    t.rts = []; t.accuracy = []; t.conds = [];
    for block = 1:t.blocks % go through each block of data
        rts = t.alldata.d.rt(block,:)';
        accuracy = t.alldata.d.correct(block,:)';
        conds = t.alldata.d.stim_mat_all(:,[5 8]);
        t.rts = [t.rts;rts];
        t.accuracy = [t.accuracy;accuracy];
        t.conds = [t.conds;conds];
    end
    clear block
    
    % strip invalid
    t.invalid = t.accuracy == -1;
    t.rts(t.invalid,:) = [];
    t.accuracy(t.invalid,:) = [];
    t.conds(t.invalid,:) = [];
    
    
%     count = 0;
%     figure;
%     for coh_idx = 1:2
%         for rule_idx = 1:2
%             count = count + 1;
%             idx = t.conds(:,1)==coh_idx & t.conds(:,2)==rule_idx;
%             idx = idx & t.accuracy;
%             subplot(2,2,count)
%             histogram(t.rts(idx),10)
%         end
%     end; clear idx coh_idx rule_idx count
%     if p.save_figs
%         export_fig(fullfile(outdir,[num2str(t.id) '_rt_hist.jpeg']),'-transparent');
%     end
    
    d.subjects(i).id = t.id;
    d.subjects(i).rts = t.rts;
    d.subjects(i).accuracy = t.accuracy;
    d.subjects(i).conditions = t.conds;
    
    t.idcolumn = ones(length(d.subjects(i).rts),1)*t.id;
    
    t.thissubj = [t.idcolumn d.subjects(i).rts, d.subjects(i).accuracy d.subjects(i).conditions];
    
    d.allsubjs = [d.allsubjs;t.thissubj];
    % d.allsubjs is participant id | rts | accuracy | coherence condition (1 = easy, 2 = hard) | rule condition (1 = easy, 2 = hard)
    
end

d.allsubjs = cell2mat(d.allsubjs); % convert it from a cell array
condlabs = {'ecer' 'echr' 'hcer' 'hchr'};
count = 0;
figure;
for coh_idx = 1:2
    for rule_idx = [2,1]
        count = count + 1;
        idx = d.allsubjs(:,4)==coh_idx & d.allsubjs(:,5)==rule_idx;
        idx = idx & d.allsubjs(:,3);
        subplot(2,2,count)
        xlabel(condlabs{count});
        h = histogram(d.allsubjs(idx,2),'FaceColor',[0.0 0.502 0.502]);
        hold on
        h.NumBins = 100;
        ylim([0, 300]);
        xlim([0,1.6])
    end
end; clear idx coh_idx rule_idx count
if p.save_figs
    export_fig(fullfile(outdir,'allsubjs_rt_hist.png'),'-transparent');
end
fprintf('done\n'); % print that so you can check


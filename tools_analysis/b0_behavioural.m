%% set up

clear all
% set up paths
addpath /hpc-software/matlab/cbu/

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';
datadir = fullfile(rootdir,'data','meg_pilot_4'); addpath(datadir);
behavfile = fullfile(datadir,'behavioural','%s_MEGRTs.mat');
toolsdir = fullfile(rootdir,'tools_analysis');
addpath(genpath(fullfile(toolsdir,'lib')))
brainDataFile = ['Preprocess' filesep 'timelocked_averages.mat'];


toolbox = fullfile(rootdir,'..','..','Toolboxes','gramm'); addpath(toolbox); clear toolbox
toolbox = fullfile(rootdir,'..','..','Toolboxes','bayesFactor'); addpath(toolbox); clear toolbox

%% load

% since our folders are all e.g. S01, S12 etc, let's load information about
% the directories that have that pattern
tmp = dir([datadir filesep 'S*']);
subjectFolders = tmp([tmp.isdir]); clear tmp
% and loop through them
for subjectNum = 1:numel(subjectFolders)
    clear thisDataFile thisBehavFile thisSubjData
    
    subjectId = subjectFolders(subjectNum).name;
    fprintf('\nthis is subject %s\n',subjectId)

    % full path to the behavioural file
    thisBehavFile = sprintf(behavfile,subjectId);
    % skip this loop if we can't find one
    if ~exist(thisBehavFile,'file'); continue; end
    % check there is brain data for this subject
    thisDataFile = sprintf(subjectFolders(subjectNum).folder,subjectFolders(subjectNum).name,brainDataFile);
    % skip this loop if we can't find one
    if ~exist(thisDataFile,'file'); continue; end
    
    % load it in
    tmp = load(thisBehavFile);
    thisSubjData = tmp.megBehaviouralData;
    clear tmp
    
    % re-composed matrix looks like:
    %  1)  cue direction (1-4)
    %  2)  dot motion direction condition (1-8)
    %  3)  coherence difficulty (1 = easy, 2 = hard)
    %  4)  matching difficulty (1 = easy, 2 = difficult)
    %  5)  unique number for each trial condition (1-64)
    %  6)  gives a number based on 9 for meg triggers
    %  7)  conditions (HcHr, HcLr, LcHr, LcLr)
    %  8)  accuracy (0, 1, or -1 for missed trials)
    %  9)  rts from behavioural data
    % 10)  rts from MEG
    % 11)  something to tag which run this sequence of data belongs to
    % 12)  we later add the fixation dots timings for each trial
    % 13)  we later add a row of unique numbers to tag each trial

    % stack em
    if ~exist('allTrials','var'); allTrials = []; end
    allTrials = [allTrials;thisSubjData];
    
    if ~exist('a','var'); a = struct(); end
    if ~isfield(a,'matrts'); a.matrts = []; end
    if ~isfield(a,'megrts'); a.megrts = []; end
    if ~isfield(a,'correct'); a.correct = []; end
    if ~isfield(a,'conditions'); a.conditions = []; end
    if ~isfield(a,'subject'); a.subject = []; end
    
    if ~exist('s','var'); s = struct(); end
    if ~isfield(s,'matrts'); s.matrts = []; end
    if ~isfield(s,'megrts'); s.megrts = []; end
    if ~isfield(s,'correct'); s.correct = []; end
    if ~isfield(s,'conditions'); s.conditions = []; end
    if ~isfield(s,'subject'); s.subject = []; end
    
    % average them
    for thisCond = unique(thisSubjData(:,7))'
        
        thisCond = thisCond{:};
        
        fprintf('this is condition %s\n',thisCond)

        idx = find(strcmp(thisSubjData(:,7),thisCond));
        lIdx = strcmp(thisSubjData(:,7),thisCond);
                
        a.matrts = [a.matrts;str2double(thisSubjData(idx,9))];
        a.megrts = [a.megrts;str2double(thisSubjData(idx,10))];
        a.correct = [a.correct;str2double(thisSubjData(idx,8))];
        a.conditions = [a.conditions;repmat({thisCond},numel(idx),1)];
        a.subject = [a.subject;repmat({subjectId},numel(idx),1)];

        s.matrts = [s.matrts;mean(str2double(thisSubjData(lIdx & str2double(thisSubjData(:,9)) > 0,9)))];
        s.megrts = [s.megrts;mean(str2double(thisSubjData(lIdx & str2double(thisSubjData(:,10)) > 0,10)))];
        cleanedCorr = str2double(thisSubjData(lIdx & str2double(thisSubjData(:,8)) > -1,8));
        s.correct = [s.correct;sum(cleanedCorr)/numel(cleanedCorr)*100];
        s.conditions = [s.conditions;{thisCond}];
        s.subject = [s.subject;{subjectId}];

    end
    
end


close all
clear g


g(1,1) = gramm('x',s.conditions,'y',s.megrts);
g(1,1).stat_summary('geom',{'bar' 'black_errorbar'},'type','sem','dodge',0.3,'width',0.3)
% g(1,1).axe_property('YLim',[500 850]);
% g(1,1).facet_grid([],r.fullConditions);
g(1,1).set_names('x','Conditions','y','Mean RT (ms) with SEM','column','','row','','color','Mode','lightness','Congruency');
g(1,1).set_title('a) Mean RTs');
% g(1,1).set_color_options('map','d3_20','legend','expand');
% figure('Position',[100 100 2000 1000]);
g.draw();

[a b c d] = ttest(s.megrts(strcmp(s.conditions,'EcEr') | strcmp(s.conditions,'EcHr')),s.megrts(strcmp(s.conditions,'HcEr') | strcmp(s.conditions,'HcHr')))
fprintf('%.3f\n',b)

[a b c d] = ttest(s.megrts(strcmp(s.conditions,'EcEr') | strcmp(s.conditions,'HcEr')),s.megrts(strcmp(s.conditions,'EcHr') | strcmp(s.conditions,'HcHr')))
fprintf('%.3f\n',b)



















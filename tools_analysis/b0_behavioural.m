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
figDir = fullfile(datadir,'behavFigs'); if ~exist(figDir,'dir'); mkdir(figDir); end


toolbox = fullfile(rootdir,'..','..','Toolboxes','gramm'); addpath(toolbox); clear toolbox
toolbox = fullfile(rootdir,'..','..','Toolboxes','bayesFactor'); addpath(toolbox); clear toolbox


teal = [0.2, 0.6, 0.7];
coral = [0.9, 0.4, 0.3];
lilac = [0.7, 0.5, 0.8];

n_colors = 3;
n_lightness = 10;

r = [linspace(teal(1), coral(1), n_lightness); linspace(coral(1), lilac(1), n_lightness); linspace(lilac(1), teal(1), n_lightness)];
g = [linspace(teal(2), coral(2), n_lightness); linspace(coral(2), lilac(2), n_lightness); linspace(lilac(2), teal(2), n_lightness)];
b = [linspace(teal(3), coral(3), n_lightness); linspace(coral(3), lilac(3), n_lightness); linspace(lilac(3), teal(3), n_lightness)];

colormap_matrix = reshape([r(:), g(:), b(:)], [], 3);

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
    if ~isfield(s,'coherence'); s.coherence = []; end
    if ~isfield(s,'categorisation'); s.categorisation = []; end
    
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
        if startsWith(thisCond,'Ec')
            s.coherence = [s.coherence;{'Easy Coh'}];
        elseif startsWith(thisCond,'Hc')
            s.coherence = [s.coherence;{'Hard Coh'}];
        else
            error('what the hell')
        end
        if endsWith(thisCond,'Er')
            s.categorisation = [s.categorisation;{'Easy Cat'}];
        elseif endsWith(thisCond,'Hr')
            s.categorisation = [s.categorisation;{'Hard Cat'}];
        end
        s.subject = [s.subject;{subjectId}];

    end
    
end

%%

close all
clear g
g(1,1) = gramm('x',s.conditions,'y',s.megrts);
g(1,1).stat_summary('geom',{'bar' 'black_errorbar'},'type','sem','dodge',0.3,'width',0.3)
g(1,1).axe_property('YLim',[650 950]);
g(1,1).set_names('x','Conditions','y','Mean RT (ms) with SEM');
g(1,1).set_title('Mean RTs for each condition');
g(1,1).set_color_options('map',teal);
% figure('Position',[100 100 2000 1000]);
g(1,2) = gramm('x',s.coherence,'y',s.megrts);
g(1,2).stat_summary('geom',{'bar' 'black_errorbar'},'type','sem','dodge',0.3,'width',0.3)
g(1,2).axe_property('YLim',[650 950]);
g(1,2).set_names('x','Conditions','y','Mean RT (ms) with SEM');
g(1,2).set_title('Mean RTs for Coherence Difficulty');
g(1,2).set_color_options('map',coral);
g(1,3) = gramm('x',s.categorisation,'y',s.megrts);
g(1,3).stat_summary('geom',{'bar' 'black_errorbar'},'type','sem','dodge',0.3,'width',0.3)
g(1,3).axe_property('YLim',[650 950]);
g(1,3).set_names('x','Conditions','y','Mean RT (ms) with SEM');
g(1,3).set_title('Mean RTs for Categorisation Difficulty');
g(1,3).set_color_options('map',lilac);
g.draw();
f = gcf; f.Position = [10 10 1600 1600];
print([figDir filesep 'rts.png'],'-dpng')

disp('are easy and hard coherence different?')
[a b c d] = ttest(s.megrts(strcmp(s.coherence,'Easy Coh')),s.megrts(strcmp(s.coherence,'Hard Coh')));
fprintf('%.3f\n',b)
disp(d)
disp('are easy and hard categorisation different?')
[a b c d] = ttest(s.megrts(strcmp(s.categorisation,'Easy Cat')),s.megrts(strcmp(s.categorisation,'Hard Cat')));
fprintf('%.3f\n',b)
disp(d)

close all
clear g
g(1,1) = gramm('x',s.conditions,'y',s.correct);
g(1,1).stat_summary('geom',{'bar' 'black_errorbar'},'type','sem','dodge',0.3,'width',0.3)
g(1,1).axe_property('YLim',[60 90]);
g(1,1).set_names('x','Conditions','y','Mean Accuracy (%) with SEM');
g(1,1).set_title('Mean Accuracy for each condition');
g(1,1).set_color_options('map',teal);
% figure('Position',[100 100 2000 1000]);
g(1,2) = gramm('x',s.coherence,'y',s.correct);
g(1,2).stat_summary('geom',{'bar' 'black_errorbar'},'type','sem','dodge',0.3,'width',0.3)
% g(1,2).axe_property('YLim',[650 950]);
g(1,2).set_names('x','Conditions','y','Mean Accuracy (%) with SEM');
g(1,2).set_title('Mean Accuracy for Coherence Difficulty');
g(1,2).set_color_options('map',coral);
g(1,3) = gramm('x',s.categorisation,'y',s.correct);
g(1,3).stat_summary('geom',{'bar' 'black_errorbar'},'type','sem','dodge',0.3,'width',0.3)
% g(1,3).axe_property('YLim',[650 950]);
g(1,3).set_names('x','Conditions','y','Mean Accuracy (%) with SEM');
g(1,3).set_title('Mean Accuracy for Categorisation Difficulty');
g(1,3).set_color_options('map',lilac);
g.draw();
f = gcf; f.Position = [10 10 1600 1600];
print([figDir filesep 'acc.png'],'-dpng')

disp('are easy and hard coherence different?')
[a b c d] = ttest(s.correct(strcmp(s.coherence,'Easy Coh')),s.correct(strcmp(s.coherence,'Hard Coh')));
fprintf('%.3f\n',b)
disp(d)
disp('are easy and hard categorisation different?')
[a b c d] = ttest(s.correct(strcmp(s.categorisation,'Easy Cat')),s.correct(strcmp(s.categorisation,'Hard Cat')));
fprintf('%.3f\n',b)
disp(d)

writetable(makeTableWithNans(...
    {'EcEr' 'EcHr' 'HcEr' 'HcHr'},...
    s.megrts(strcmp(s.conditions,'EcEr')),...
    s.megrts(strcmp(s.conditions,'EcHr')),...
    s.megrts(strcmp(s.conditions,'HcEr')),...
    s.megrts(strcmp(s.conditions,'HcHr'))),...
    [datadir filesep 'rts-for-jasp.csv'])

writetable(makeTableWithNans(...
    {'EcEr' 'EcHr' 'HcEr' 'HcHr'},...
    s.correct(strcmp(s.conditions,'EcEr')),...
    s.correct(strcmp(s.conditions,'EcHr')),...
    s.correct(strcmp(s.conditions,'HcEr')),...
    s.correct(strcmp(s.conditions,'HcHr'))),...
    [datadir filesep 'accuracy-for-jasp.csv'])

close all
clear g
g(1,1) = gramm('x',s.coherence,'y',s.megrts,'color',s.categorisation);
g(1,1).stat_summary('geom',{'line' 'errorbar'},'type','sem','dodge',0.3,'width',0.3)
g(1,1).axe_property('YLim',[650 950]);
g(1,1).set_names('x','Coherence','y','Mean RTs (ms) with SEM','column','','row','','color','Categorisation','lightness','');
g(1,1).set_title('Mean RTs by Condition');
g(1,1).set_color_options('map',colormap_matrix);
figure('Position',[100 100 2000 1000]);
g(1,2) = gramm('x',s.coherence,'y',s.correct,'color',s.categorisation);
g(1,2).stat_summary('geom',{'line' 'errorbar'},'type','sem','dodge',0.3,'width',0.3)
g(1,2).axe_property('YLim',[60 90]);
g(1,2).set_names('x','Coherence','y','Mean Accuracy (%) with SEM','column','','row','','color','Categorisation','lightness','');
g(1,2).set_title('Mean Accuracy by Condition');
g(1,2).set_color_options('map',colormap_matrix);
% figure('Position',[100 100 2000 1000]);
g.draw();

[a b c d] = ttest(s.megrts(strcmp(s.coherence,'Easy Coh') & strcmp(s.categorisation,'Easy Cat')),...
    s.megrts(strcmp(s.coherence,'Hard Coh') & strcmp(s.categorisation,'Easy Cat')),'alpha',0.05/4);
fprintf('%.3f\n',b)
[a b c d] = ttest(s.megrts(strcmp(s.coherence,'Easy Coh') & strcmp(s.categorisation,'Hard Cat')),...
    s.megrts(strcmp(s.coherence,'Hard Coh') & strcmp(s.categorisation,'Hard Cat')),'alpha',0.05/4);
fprintf('%.3f\n',b)
disp(d)
[a b c d] = ttest(s.megrts(strcmp(s.categorisation,'Easy Cat') & strcmp(s.coherence,'Easy Coh')),...
    s.megrts(strcmp(s.categorisation,'Hard Cat') & strcmp(s.coherence,'Easy Coh')),'alpha',0.05/4);
fprintf('%.3f\n',b)
disp(d)
[a b c d] = ttest(s.megrts(strcmp(s.categorisation,'Easy Cat') & strcmp(s.coherence,'Hard Coh')),...
    s.megrts(strcmp(s.categorisation,'Hard Cat') & strcmp(s.coherence,'Hard Coh')),'alpha',0.05/4);
fprintf('%.3f\n',b)
disp(d)
















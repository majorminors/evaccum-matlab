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

% since our folders are all e.g. S01, S12 etc, let's load information about
% the directories that have that pattern
tmp = dir([datadir filesep 'S*']);
subjectFolders = tmp([tmp.isdir]); clear tmp
% and loop through them
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
    data{subjectNum} = load(thisFile);
    disp('loaded')
end

disp('done loading')

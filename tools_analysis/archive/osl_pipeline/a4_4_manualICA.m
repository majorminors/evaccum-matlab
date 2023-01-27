function a4_4_manualICA(thisSubject,datadir,toolsdir,runLocal)

addpath /hpc-software/matlab/cbu/
addpath('/neuro/meg_pd_1.2/');       % FIFACCESS toolbox if use some meg_misc functions below
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/');
addpath(fullfile(toolsdir,'fieldtrip-20190410')); ft_defaults; % not needed if running locally, but needed if running on the cluster for some reason
% addpath(genpath(fullfile(toolsdir,'eeglab13_5_4b')))
addpath(genpath(fullfile(toolsdir,'meg_data')))
oslDir = fullfile('/imaging/woolgar/projects/Dorian/Toolboxes/osl'); % osl wants this if run in shared user mode
addpath(fullfile(oslDir,'osl-core')); % osl_startup can error without this
if runLocal; oslUserMode = 'user'; elseif ~runLocal; oslUserMode = 'shared';end

rootDir = fullfile(datadir,thisSubject.id);
inputFolder  = fullfile(rootDir,'Preprocess'); % where should we find the files?
outputFolder  = fullfile(rootDir,'Preprocess'); % where should we put the output?
behavDir   = fullfile(datadir,'behavioural'); % where is our behavioural data?
behavfile = fullfile(behavDir,[thisSubject.id '_Evaccum.mat']); % what is it called?
megrtfile = fullfile(behavDir,[thisSubject.id '_MEGRTs.mat']); % and since I have put my MEG triggers somewhere else, we will specify that too
trlfile  = fullfile(outputFolder, 'run%s_trl.mat'); % this contains our trial related info for epoching - should exist (generated with megTriggers script)
addpath(rootDir);

%% manual ICA

disp('now doing manual ica')
% WIP: this does a manual identification, but I haven't
% checked it works/saves properly and whatnot

inputFileInfo = dir(fullfile(inputFolder,'ica*.mat'));
inputFile = [inputFileInfo.folder filesep inputFileInfo.name];

D = spm_eeg_load(inputFile); % load this, in case we're just re-running for doing manual checks
cfg = [];
cfg.output = [outputFolder filesep 'eeg_layout.lay'];
cfg.elec = D.sensors('EEG');
eeg_layout = ft_prepare_layout(cfg);

% we then want to redo the artefact rejection to make changes to the assignment if we like
D = osl_africa(D,...
    'do_ica',false,...
    'do_ident','manual',...
    'used_maxfilter',true,...
    'artefact_channels',{'EOG','ECG'},...
    'montagename','ols africa ica denoised meg'...
    );
% now will look for EEG
D = osl_africa(D,...
    'do_ica',false,...
    'do_ident','manual',...
    'eeg_layout',eeg_layout.cfg.output,...
    'modality',{'EEG'},...
    'used_maxfilter',true,...
    'artefact_channels',{'EOG','ECG'},...
    'montagename','ols africa ica denoised eeg'...
    );

D.save(); % we do need to save the MEEG object to commit marked bad components to disk after the manual ica

% osl africa saves an 'online montage' (a linear combination of the original sensor data) that we can explore
D = spm_eeg_load(inputFile); % load our first dataset
has_montage(D) % should show available montages

D_pre_ica = D.montage('switch',0); % grab the original data
D_post_ica = D.montage('switch',2); % grab the data transformed by the ica

%oslview(D_raw); % if you'd like to compare these
%oslview(D_filtered); % if you'd like to compare these
%oslview(D_deartefacted); % if you'd like to compare these

% but now we can start to plot interesting stuff
artefactSensor = 2; % change this to an ECG or EOG sensor - something that is artefactual
sensorToCheck = 4; % change this to a sensor near your artefact sensor - we will compare it to the artefact before and after ICA
% note that the way we've done this means there are two ica montages,
% one for eeg and one for meg, so if we have selected the wrong
% montage+sensor combination, the plots will only show the artefact
% sensor
figure;
subplot(2,1,1);
plot(D_pre_ica.time(1:10000),D_pre_ica(artefactSensor,1:10000)); % takes first 10000 sample points
title('artefact channel')
xlim([0 10]);
xlabel('Time (s)');

subplot(2,1,2);
plot(D_pre_ica.time(1:10000),D_pre_ica(sensorToCheck,1:10000)); % takes first 10000 sample points
title('artefact contaminated channel')
xlim([0 10]);
hold on;
plot(D_post_ica.time(1:10000),D_post_ica(sensorToCheck,1:10000),'r');
xlim([0 10]);
xlabel('Time (s)');
legend({'pre ICA' 'post ICA'});

clear artefactSensor sensorToCheck
close all

return
end
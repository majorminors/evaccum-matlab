% here we will:
% 1. view the raw values of onsets and responses by subject: check it looks
% sensible
% made this while investigating the results, then re-wrote most of this into
% a subject-based function to run on the cluster (b2)

clear all

% rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';
rootdir = pwd;
datadir = fullfile(rootdir,'data','meg_pilot_4'); addpath(datadir);
toolsdir = fullfile(rootdir,'tools_analysis');
ftDir = fullfile(toolsdir,'..','..','..','Toolboxes','fieldtrip');

addpath(genpath(toolsdir))

data_inspection(datadir, toolsdir, ftDir, 1,{'if' '.fif' ''})

% sendToCluster(@data_inspection, {datadir,toolsdir,ftDir,1,{'if' '.fif' ''}}, {rootdir datadir toolsdir})


function data_inspection(datadir, toolsdir, ftDir, overwrite,fileMods)

if ~exist('overwrite','var'); overwrite = 0; end
if ~exist('fileMods','var'); fileMods = {'' '.mat' 'v'}; end


% set up paths
addpath /hpc-software/matlab/cbu/
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/');

addpath(ftDir); ft_defaults;
addpath(genpath(fullfile(toolsdir,'lib')))

% file definition
inputFileName = ['Preprocess' filesep 'run*' fileMods{1} '_transmid' fileMods{2}];
trlfile = '%s_trl.mat'; % this contains our trial related info for epoching - should exist (generated with megTriggers script)
behavfile = fullfile(datadir,'behavioural','%s_MEGRTs.mat'); % grab our behavioural data + output from a3_megtriggers
saveFile = fullfile(datadir,'timelockedData');

% first loop through subjects

% since our folders are all e.g. S01, S12 etc, let's load information about
% the directories that have that pattern
tmp = dir([datadir filesep 'S*']);
subjectFolders = tmp([tmp.isdir]); clear tmp
% and loop through them
for subjectNum = 1:numel(subjectFolders)
    disp('this is subject:')
    disp(subjectFolders(subjectNum).name);
    % grab info about all the meeg data files we care about
    theseFiles = dir([subjectFolders(subjectNum).folder filesep subjectFolders(subjectNum).name filesep inputFileName]);
    % skip this loop if we don't have any files to work with
    if isempty(theseFiles); continue; end
    % and the behavioural data for this subject
    thisBehavFile = sprintf(behavfile,subjectFolders(subjectNum).name);
    % and make a figure directory
    thisFigDir = fullfile(subjectFolders(subjectNum).folder, subjectFolders(subjectNum).name,'figs');
    if ~exist(thisFigDir,'dir'); mkdir(thisFigDir); end
    % loop through those files
    for fileNum = 1:numel(theseFiles)
        % let's make the current data file easy to reference
        thisFile = [theseFiles(fileNum).folder filesep theseFiles(fileNum).name];
        % let's find what the run label is (e.g. 'run1' etc) (looks for
        % 'run' followed by any amount of numerical digits, and returns
        % {'run[numbers]'}, so we pop it out of the cell array. This would
        % break if you had more than one 'run[numbers]' in your filename
        % because it would return more than one item in the cell array
        runLabel = regexp(thisFile,'run\d+','match'); runLabel = runLabel{:};
        % and get the run number in a similar way, since I later found out I need this too
        thisRun = regexp(thisFile,'run(\d*)_','tokens'); thisRun = thisRun{1}{1};
        % now we can find the trial information for this run
        thisTrlFile = [theseFiles(fileNum).folder filesep sprintf(trlfile,runLabel)];
        
        % ok, now let's load in the data file
        cfg = [];
        cfg.continuous = 'yes'; % this data is not epoched, it's continuous
        if strcmp(fileMods{2},'.fif')
            cfg.dataset = thisFile;
            rawData = convertFif(cfg,ftDir); % fieldtrip clears the mne toolbox after first use? so let's force this in a distinct environment, in case this is something they did on purpose
        elseif strcmp(fileMods{2},'.mat')
            rawData = ft_preprocessing(cfg,loadFtData(thisFile)); % load it
        end
        
        % let's make a layout using the data
        [combinedLayout, megLayout, eegLayout] = makeLayout(rawData);
        
        % and now lets break it into trials using our custom trial function
        % (this trial function equally works with older versions of
        % fieldtrip---just use cfg.trialfun = 'trlFromFile', instead of
        % feeding the output directly to cfg.trl (although maybe this would
        % work even for older versions of ft?)
        
        %% first define data by the onset of fixation dots
        cfg = [];
        cfg.trlfile = thisTrlFile;
        cfg.behavfile = thisBehavFile;
        cfg.run = thisRun;
        cfg.lockTo = 'fixation';
        cfg.trl = trlFromFile(cfg);
        fixOnsetsData{fileNum} = ft_redefinetrial(cfg,rawData);
        fixOnsetsData{fileNum} = rejectArtefacts(fixOnsetsData{fileNum}, combinedLayout,ftDir);
        
        % so with this, you can select e.g. MEG, then in the interactive
        % figure from multiplot, you can draw a box around the channels you
        % care about for this effect, click in the box to see the average,
        % then copy the title of the results to just grab those channels in
        % future;
        %         cfg = [];
        %         cfg.trials = find(fixOnsetsData{fileNum}.trialinfo(:,2) == 1); % grab just the good trials
        %         fixOnsets{subjectNum} = ft_timelockanalysis(cfg, fixOnsetsData{fileNum}); clear tmpData
        %         cfg = [];
        %         cfg.showlabels = 'yes';
        %         cfg.fontsize = 6;
        %         cfg.layout = megLayout;
        %         cfg.channel = 'MEG';
        %         ft_multiplotER(cfg, fixOnsets{subjectNum});
        %         cfg.layout = eegLayout;
        %         cfg.channel = 'EEG';
        %         ft_multiplotER(cfg, fixOnsets{subjectNum});
        
        %% now onset of coherent dots
        cfg = [];
        cfg.trlfile = thisTrlFile;
        cfg.behavfile = thisBehavFile;
        cfg.run = thisRun;
        cfg.trl = trlFromFile(cfg);
        cohOnsetsData{fileNum} = ft_redefinetrial(cfg,rawData);
        cohOnsetsData{fileNum} = rejectArtefacts(cohOnsetsData{fileNum}, combinedLayout,ftDir);
        
        
        %% now response locked
        cfg = [];
        cfg.trlfile = thisTrlFile;
        cfg.behavfile = thisBehavFile;
        cfg.run = thisRun;
        cfg.lockTo = 'response';
        cfg.pre = -600;
        cfg.post = 200;
        cfg.trl = trlFromFile(cfg);
        responseLockedData{fileNum} = ft_redefinetrial(cfg,rawData);
        responseLockedData{fileNum} = rejectArtefacts(responseLockedData{fileNum}, combinedLayout,ftDir);
        
        
        
    end
    
    %% now graph at subject level
    
    %% first fixations
    
    % append them all
    cfg = [];
    cfg.keepsampleinfo  = 'yes';
    tmpData = ft_appenddata(cfg,fixOnsetsData{:});
    clear fixOnsetsData
    
    cfg = [];
    cfg.trials = find(tmpData.trialinfo(:,2) == 1); % grab just the good trials
    fixOnsets{subjectNum} = ft_timelockanalysis(cfg, tmpData);
    clear tmpData
    
    cfg = [];
    cfg.showlabels = 'yes';
    cfg.fontsize = 6;
    
    % let's plot the mean of some occipital channels
    occiesMeg = {'MEG1711' 'MEG1731' 'MEG1732' 'MEG1733' 'MEG1741' 'MEG1742' 'MEG1743' 'MEG1922' 'MEG1923' 'MEG1931' 'MEG1932' 'MEG1933' 'MEG2112' 'MEG2113' 'MEG2121' 'MEG2122' 'MEG2123' 'MEG2131' 'MEG2132' 'MEG2133' 'MEG2141' 'MEG2142' 'MEG2143' 'MEG2332' 'MEG2333' 'MEG2342' 'MEG2343' 'MEG2543'};
    occiesEeg = {'EEG001' 'EEG003' 'EEG018' 'EEG028' 'EEG029' 'EEG039' 'EEG051' 'EEG052' 'EEG053' 'EEG054' 'EEG055' 'EEG056' 'EEG057' 'EEG058' 'EEG059' 'EEG060' 'EEG061' 'EEG062' 'EEG063' 'EEG064'};
    cfg.channel = occiesMeg;
    ft_singleplotER(cfg, fixOnsets{subjectNum});
    %         theseLims = xlim(); xlim([theseLims(1) 0.5]); clear theseLims
    print([thisFigDir filesep 'meg_fixation_occave.png'], '-dpng');
    cfg.channel = occiesEeg;
    ft_singleplotER(cfg, fixOnsets{subjectNum});
    %         theseLims = xlim(); xlim([theseLims(1) 0.5]); clear theseLims
    print([thisFigDir filesep 'eeg_fixation_occave.png'], '-dpng');
    
    close all
    
    %% now coherence
    
    % append them all
    cfg = [];
    cfg.keepsampleinfo  = 'yes';
    tmpData = ft_appenddata(cfg,cohOnsetsData{:});
    clear cohOnsetsData
    
    cfg = [];
    cfg.trials = find(tmpData.trialinfo(:,2) == 1); % grab just the good trials
    cohOnsets{subjectNum} = ft_timelockanalysis(cfg, tmpData);
    clear tmpData
    
    cfg = [];
    cfg.showlabels = 'yes';
    cfg.fontsize = 6;
    
    % let's plot the mean of some occipital channels
    cfg.channel = occiesMeg;
    ft_singleplotER(cfg, cohOnsets{subjectNum});
    print([thisFigDir filesep 'meg_onset_occave.png'], '-dpng');
    cfg.channel = occiesEeg;
    ft_singleplotER(cfg, cohOnsets{subjectNum});
    print([thisFigDir filesep 'eeg_onset_occave.png'], '-dpng');
    
    close all
    
    %% now response locked
    
    % append them all
    cfg = [];
    cfg.keepsampleinfo  = 'yes';
    tmpData = ft_appenddata(cfg,responseLockedData{:});
    clear responseLockedData
    
    cfg = [];
    cfg.trials = find(tmpData.trialinfo(:,2) == 1); % grab just the good trials
    respLocked{subjectNum} = ft_timelockanalysis(cfg, tmpData);
    
    cfg = [];
    cfg.showlabels = 'no';
    cfg.fontsize = 2;
    
    % let's plot the mean of some midline parietal channels
    midParMeg = {'MEG0431' 'MEG0632' 'MEG0633' 'MEG0711' 'MEG0712' 'MEG0713' 'MEG0721' 'MEG0722' 'MEG0723' 'MEG0732' 'MEG0733' 'MEG0741' 'MEG0742' 'MEG0743' 'MEG1042' 'MEG1043' 'MEG1113' 'MEG1143' 'MEG1821' 'MEG1823' 'MEG1831' 'MEG1832' 'MEG1833' 'MEG1841' 'MEG1842' 'MEG2011' 'MEG2021' 'MEG2241' 'MEG2242'};
    midParEeg = {'EEG022' 'EEG023' 'EEG024' 'EEG033' 'EEG034' 'EEG035' 'EEG044' 'EEG045' 'EEG046' 'EEG055' 'EEG056' 'EEG057'};
    cfg.channel = midParMeg;
    ft_singleplotER(cfg, respLocked{subjectNum});
    hold on
    plot([0 0], ylim, 'k--')
    hold off
    print([thisFigDir filesep 'meg_resp_midparave.png'], '-dpng');
    cfg.channel = midParEeg;
    ft_singleplotER(cfg, respLocked{subjectNum})
    hold on
    plot([0 0], ylim, 'k--')
    hold off
    print([thisFigDir filesep 'eeg_resp_midparave.png'], '-dpng');
    
    close all
    
    % now the same thing but broken into our conditions
    cfg = [];
    cfg.trials = find(tmpData.trialinfo(:,2) == 1 & tmpData.trialinfo(:,1) == 1); % good trials, 1st condition
    ecer_respLocked{subjectNum} = ft_timelockanalysis(cfg, tmpData);
    cfg.trials = find(tmpData.trialinfo(:,2) == 1 & tmpData.trialinfo(:,1) == 2); % good trials, 1st condition
    echr_respLocked{subjectNum} = ft_timelockanalysis(cfg, tmpData);
    cfg.trials = find(tmpData.trialinfo(:,2) == 1 & tmpData.trialinfo(:,1) == 3); % good trials, 1st condition
    hcer_respLocked{subjectNum} = ft_timelockanalysis(cfg, tmpData);
    cfg.trials = find(tmpData.trialinfo(:,2) == 1 & tmpData.trialinfo(:,1) == 4); % good trials, 1st condition
    hchr_respLocked{subjectNum} = ft_timelockanalysis(cfg, tmpData);
    
    cfg.trials = find(tmpData.trialinfo(:,2) == 1 & (tmpData.trialinfo(:,1) == 1 | tmpData.trialinfo(:,1) == 2)); % good trials, easy coherence
    ec_respLocked{subjectNum} = ft_timelockanalysis(cfg, tmpData);
    cfg.trials = find(tmpData.trialinfo(:,2) == 1 & (tmpData.trialinfo(:,1) == 3 | tmpData.trialinfo(:,1) == 4)); % good trials, hard coherence
    hc_respLocked{subjectNum} = ft_timelockanalysis(cfg, tmpData);
    
    cfg.trials = find(tmpData.trialinfo(:,2) == 1 & (tmpData.trialinfo(:,1) == 1 | tmpData.trialinfo(:,1) == 3)); % good trials, easy rule
    er_respLocked{subjectNum} = ft_timelockanalysis(cfg, tmpData);
    cfg.trials = find(tmpData.trialinfo(:,2) == 1 & (tmpData.trialinfo(:,1) == 2 | tmpData.trialinfo(:,1) == 4)); % good trials, hard rule
    hr_respLocked{subjectNum} = ft_timelockanalysis(cfg, tmpData);
    
    clear tmpData
    
    cfg = [];
    cfg.showlabels = 'no';
    cfg.fontsize = 2;
    
    cfg.channel = midParMeg;
    ft_singleplotER(cfg, ecer_respLocked{subjectNum},echr_respLocked{subjectNum},hcer_respLocked{subjectNum},hchr_respLocked{subjectNum});
    hold on
    plot([0 0], ylim, 'k--')
    legend({'ecer' 'echr' 'hcer' 'hchr'},'Location','southwest')
    hold off
    print([thisFigDir filesep 'meg_resp_midparCondave.png'], '-dpng');
    cfg.channel = midParEeg;
    ft_singleplotER(cfg, ecer_respLocked{subjectNum},echr_respLocked{subjectNum},hcer_respLocked{subjectNum},hchr_respLocked{subjectNum});
    hold on
    plot([0 0], ylim, 'k--')
    legend({'ecer' 'echr' 'hcer' 'hchr'},'Location','southwest')
    hold off
    print([thisFigDir filesep 'eeg_resp_midparCondave.png'], '-dpng');
    
    cfg.channel = midParMeg;
    ft_singleplotER(cfg, ec_respLocked{subjectNum},hc_respLocked{subjectNum},er_respLocked{subjectNum},hr_respLocked{subjectNum});
    hold on
    plot([0 0], ylim, 'k--')
    legend({'ec' 'hc' 'er' 'hr'},'Location','southwest')
    hold off
    print([thisFigDir filesep 'meg_resp_midparTypeave.png'], '-dpng');
    cfg.channel = midParEeg;
    ft_singleplotER(cfg, ec_respLocked{subjectNum},hc_respLocked{subjectNum},er_respLocked{subjectNum},hr_respLocked{subjectNum});
    hold on
    plot([0 0], ylim, 'k--')
    legend({'ec' 'hc' 'er' 'hr'},'Location','southwest')
    hold off
    print([thisFigDir filesep 'eeg_resp_midparTypeave.png'], '-dpng');
    
    close all
    
end

% first let's save that in case it crashes
% disp('saving')
% save(saveFile,'fixOnsets','cohOnsets','respLocked',...
%     'ecer_respLocked','echr_respLocked','hcer_respLocked','hchr_respLocked',...
%     '-v7.3');
% % save(saveFile,'fixOnsets','cohOnsets','respLocked',...
% %     'ecer_respLocked','echr_respLocked','hcer_respLocked','hchr_respLocked',...
% %     'ec_respLocked','hc_respLocked','er_respLocked','hr_respLocked',...
% %     '-v7.3');
% disp('saved')

getFull = @(x) find(~cellfun(@isempty,x));
grandFigdir = [saveFile '_figs'];
if ~exist(grandFigdir); mkdir(grandFigdir); end

%% now the same thing but group level!
cfg = [];
fixOnsetsGrand = ft_timelockgrandaverage(cfg, fixOnsets{getFull(fixOnsets)});
cfg = [];
cfg.showlabels = 'yes';
cfg.fontsize = 6;
% let's plot the mean of some occipital channels
occiesMeg = {'MEG1711' 'MEG1731' 'MEG1732' 'MEG1733' 'MEG1741' 'MEG1742' 'MEG1743' 'MEG1922' 'MEG1923' 'MEG1931' 'MEG1932' 'MEG1933' 'MEG2112' 'MEG2113' 'MEG2121' 'MEG2122' 'MEG2123' 'MEG2131' 'MEG2132' 'MEG2133' 'MEG2141' 'MEG2142' 'MEG2143' 'MEG2332' 'MEG2333' 'MEG2342' 'MEG2343' 'MEG2543'};
occiesEeg = {'EEG001' 'EEG003' 'EEG018' 'EEG028' 'EEG029' 'EEG039' 'EEG051' 'EEG052' 'EEG053' 'EEG054' 'EEG055' 'EEG056' 'EEG057' 'EEG058' 'EEG059' 'EEG060' 'EEG061' 'EEG062' 'EEG063' 'EEG064'};
cfg.channel = occiesMeg;
ft_singleplotER(cfg, fixOnsetsGrand);
print([grandFigdir filesep 'meg_fixation_occave.png'], '-dpng');
cfg.channel = occiesEeg;
ft_singleplotER(cfg, fixOnsetsGrand);
print([grandFigdir filesep 'eeg_fixation_occave.png'], '-dpng');
close all

cfg = [];
cohOnsetsGrand = ft_timelockgrandaverage(cfg, cohOnsets{getFull(cohOnsets)});
cfg = [];
cfg.showlabels = 'yes';
cfg.fontsize = 6;
% let's plot the mean of some occipital channels
cfg.channel = occiesMeg;
ft_singleplotER(cfg, cohOnsetsGrand);
print([grandFigdir filesep 'meg_onset_occave.png'], '-dpng');
cfg.channel = occiesEeg;
ft_singleplotER(cfg, cohOnsetsGrand);
print([grandFigdir filesep 'eeg_onset_occave.png'], '-dpng');

close all
    
cfg = [];
respLockedGrand = ft_timelockgrandaverage(cfg, respLocked{getFull(respLocked)});
cfg = [];
cfg.showlabels = 'no';
cfg.fontsize = 2;
% let's plot the mean of some midline parietal channels
midParMeg = {'MEG0431' 'MEG0632' 'MEG0633' 'MEG0711' 'MEG0712' 'MEG0713' 'MEG0721' 'MEG0722' 'MEG0723' 'MEG0732' 'MEG0733' 'MEG0741' 'MEG0742' 'MEG0743' 'MEG1042' 'MEG1043' 'MEG1113' 'MEG1143' 'MEG1821' 'MEG1823' 'MEG1831' 'MEG1832' 'MEG1833' 'MEG1841' 'MEG1842' 'MEG2011' 'MEG2021' 'MEG2241' 'MEG2242'};
midParEeg = {'EEG022' 'EEG023' 'EEG024' 'EEG033' 'EEG034' 'EEG035' 'EEG044' 'EEG045' 'EEG046' 'EEG055' 'EEG056' 'EEG057'};
cfg.channel = midParMeg;
ft_singleplotER(cfg, respLockedGrand);
hold on
plot([0 0], ylim, 'k--')
hold off
print([grandFigdir filesep 'meg_resp_midparave.png'], '-dpng');
cfg.channel = midParEeg;
ft_singleplotER(cfg, respLockedGrand)
hold on
plot([0 0], ylim, 'k--')
hold off
print([grandFigdir filesep 'eeg_resp_midparave.png'], '-dpng');
close all

cfg = [];
ecer_respLockedGrand = ft_timelockgrandaverage(cfg, ecer_respLocked{getFull(ecer_respLocked)});
cfg = [];
echr_respLockedGrand = ft_timelockgrandaverage(cfg, echr_respLocked{getFull(echr_respLocked)});
cfg = [];
hcer_respLockedGrand = ft_timelockgrandaverage(cfg, hcer_respLocked{getFull(hcer_respLocked)});
cfg = [];
hchr_respLockedGrand = ft_timelockgrandaverage(cfg, hchr_respLocked{getFull(hchr_respLocked)});
cfg = [];
cfg.showlabels = 'no';
cfg.fontsize = 2;
cfg.channel = midParMeg;
ft_singleplotER(cfg, ecer_respLockedGrand,echr_respLockedGrand,hcer_respLockedGrand,hchr_respLockedGrand);
hold on
plot([0 0], ylim, 'k--')
legend({'ecer' 'echr' 'hcer' 'hchr'},'Location','southwest')
hold off
print([grandFigdir filesep 'meg_resp_midparCondave.png'], '-dpng');
cfg.channel = midParEeg;
ft_singleplotER(cfg, ecer_respLockedGrand,echr_respLockedGrand,hcer_respLockedGrand,hchr_respLockedGrand);
hold on
plot([0 0], ylim, 'k--')
legend({'ecer' 'echr' 'hcer' 'hchr'},'Location','southwest')
hold off
print([grandFigdir filesep 'eeg_resp_midparCondave.png'], '-dpng');
close all

cfg = [];
ec_respLockedGrand = ft_timelockgrandaverage(cfg, ec_respLocked{getFull(ec_respLocked)});
cfg = [];
hc_respLockedGrand = ft_timelockgrandaverage(cfg, hc_respLocked{getFull(hc_respLocked)});
cfg = [];
er_respLockedGrand = ft_timelockgrandaverage(cfg, er_respLocked{getFull(er_respLocked)});
cfg = [];
hr_respLockedGrand = ft_timelockgrandaverage(cfg, hr_respLocked{getFull(hr_respLocked)});
cfg = [];
cfg.showlabels = 'no';
cfg.fontsize = 2;
cfg.channel = midParMeg;
ft_singleplotER(cfg, ec_respLockedGrand,hc_respLockedGrand,er_respLockedGrand,hr_respLockedGrand);
hold on
plot([0 0], ylim, 'k--')
legend({'ec' 'hc' 'er' 'hr'},'Location','southwest')
hold off
print([grandFigdir filesep 'meg_resp_midparTypeave.png'], '-dpng');
cfg.channel = midParEeg;
ft_singleplotER(cfg, ec_respLockedGrand,hc_respLockedGrand,er_respLockedGrand,hr_respLockedGrand);
hold on
plot([0 0], ylim, 'k--')
legend({'ec' 'hc' 'er' 'hr'},'Location','southwest')
hold off
print([grandFigdir filesep 'eeg_resp_midparTypeave.png'], '-dpng');
close all




return
end

function [combinedLayout, megLayout, eegLayout] = makeLayout(data)

%% prep layout
% make an eeg layout
cfg = [];
cfg.elec = data.elec;
eegLayout = ft_prepare_layout(cfg);
% and an meg layout
cfg = [];
cfg.grad = data.grad;
megLayout = ft_prepare_layout(cfg);
% combine them
cfg = [];
combinedLayout = ft_appendlayout(cfg, eegLayout, megLayout);

return
end

function cleanedData = rejectArtefacts(inputData, layout,ftDir)
addpath(ftDir); ft_defaults;

cfg = [];
% here's an example of badchannels (NO CHANNELS SPECIFIED!)
% but we won't worry about that for now
%         cfg.layout = layout;
%         cfg.method = 'distance';
%         cfg.neighbours = ft_prepare_neighbours(cfg, inputData);
%         cfg = ft_badchannel(cfg, inputData);
%         badMegChans = cfg.badchannel;
cfg.artfctdef.reject = 'complete'; % remove 'complete' trials
cfg.metric = 'zvalue';
cfg.threshold = 4;
cfg.channel = 'MEG';
cfg = ft_badsegment(cfg, inputData);
cleanedData = ft_rejectartifact(cfg, inputData);
cfg.channel = 'EEG';
cfg = ft_badsegment(cfg, cleanedData);
disp('dont worry about warnings about nans---this wont spread, its contained to the trial youre naning')
cleanedData = ft_rejectartifact(cfg, cleanedData);

return
end

function convertedData = convertFif(cfg, ftDir)
% fieldtrip clears the mne toolbox after first use? so let's force this in
% a distinct environment in case fieldtrip does this for a good reason

addpath(ftDir); ft_defaults; addpath(fullfile(ftDir, 'external','mne'));
convertedData = ft_preprocessing(cfg); % load it
            
return
end

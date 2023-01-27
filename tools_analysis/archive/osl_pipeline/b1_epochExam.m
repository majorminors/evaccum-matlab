%% epoch exam at group level

clear all


rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';
scriptdir = fullfile(rootdir,'tools_analysis'); cd(scriptdir)
datadir = fullfile(rootdir,'data','meg_pilot_4'); addpath(datadir);
behavDataDir = fullfile(datadir,'behavioural');
outputFolder = fullfile(datadir,'epochExam'); if ~exist(outputFolder,'dir'); mkdir(outputFolder); end

oslDir = fullfile('/imaging/woolgar/projects/Dorian/Toolboxes/osl'); % osl wants this if run in shared user mode
addpath(fullfile(oslDir,'osl-core')); % osl_startup can error without this

saveName = fullfile(outputFolder,'allData');
subjectRange = 0;%[-1 7]; % 0 does all; array like [-1 4] does 4 to end
runLocal = 1; % doesn't much now, but should be on for osl
disp('starting up osl')
if runLocal; oslUserMode = 'user'; elseif ~runLocal; oslUserMode = 'shared';end
osl_startup(oslDir,oslUserMode)

allSubjects = importParticipants();
if ~subjectRange; subjectRange = 1:numel(allSubjects); end
isneg = @(val) val < 0 ;
if any(isneg(subjectRange)); subjectRange = subjectRange(~isneg(subjectRange)):1:numel(allSubjects); end

counter = 0;
for subjectidx = subjectRange
    %     disp('starting with subject: '); disp(subjectidx);
    
    thisSubject = allSubjects{subjectidx}; % let's make this easier to reference
    
    % let's skip participants who are not useable
    if ~thisSubject.usable
        disp('marked not usable - skipping')
        continue
    end
    disp('this is subject')
    disp(thisSubject.id)
    counter = counter+1;
    
    subjectDir = fullfile(datadir,thisSubject.id);
    preprocFolder = fullfile(subjectDir,'Preprocess');
    behavData = fullfile(behavDataDir,[thisSubject.id '_EvAccum.mat']);
    megrtFile = fullfile(behavDataDir,[thisSubject.id '_MEGRTs.mat']);
    rawDataFile = sprintf([preprocFolder,filesep,'allRuns_trans.mat']);
    filterPrefix = 'f'; % most programs default to an f, but note we filter twice so will have a file with ff prepended
    artefactPrefix = 'b';
    icaPrefix = 'ica'; % we add this to the file ourselves
    epochPrefix = 'e'; % we also choose this
    filteredDataFile = addPrefix(rawDataFile,[filterPrefix filterPrefix]);
    deArtefactedFile = addPrefix(filteredDataFile,artefactPrefix);
    icaedFile = addPrefix(deArtefactedFile,icaPrefix);
    epochedFile = addPrefix(icaedFile,epochPrefix);
    deArtefactedEpochs = addPrefix(epochedFile,artefactPrefix);
    
    disp('loading epoched data')
    D_epoched_deArtefacted = spm_eeg_load(deArtefactedEpochs);
    
    data(counter).subject = thisSubject;
    
    tmp = load(megrtFile,'megBehaviouralData');
    data(counter).behav = tmp.megBehaviouralData;
    clear tmp
    tmp = load(behavData,'p');
    data(counter).fixationTimings = tmp.p.fixation_dots_duration_vector;
    clear tmp

    D_raw = D_epoched_deArtefacted.montage('switch',0); % grab the original data
    D_clean_meg = D_epoched_deArtefacted.montage('switch',1); % grab the data transformed by the ica
    D_clean_eeg = D_epoched_deArtefacted.montage('switch',2); % grab the data transformed by the ica
    
%     megData = spm2fieldtrip(D_clean_meg);
%     eegData = spm2fieldtrip(D_clean_eeg);
%     
%     for trialRecodeIdx = 1:size(megData.trialinfo,1)
%         if megData.trialinfo(trialRecodeIdx,1) == 1 % it coded this as EcHr
%             megData.trialinfo(trialRecodeIdx,1) = 2;
%             eegData.trialinfo(trialRecodeIdx,1) = 2;
%         elseif megData.trialinfo(trialRecodeIdx,1) == 2 % it coded this as HcEr
%             megData.trialinfo(trialRecodeIdx,1) = 3;
%             eegData.trialinfo(trialRecodeIdx,1) = 3;
%         elseif megData.trialinfo(trialRecodeIdx,1) == 3 % it coded this as EcEr
%             megData.trialinfo(trialRecodeIdx,1) = 1;
%             eegData.trialinfo(trialRecodeIdx,1) = 1;
%         elseif megData.trialinfo(trialRecodeIdx,1) == 4 % it coded this as HcHr
%             do nothing for recoding trial
%         else
%             error('codes not right');
%         end
%     end
%     
%     megData.trialinfo(indtrial(D_raw,'EcEr','good'),2) = 1;
%     eegData.trialinfo(indtrial(D_raw,'EcEr','good'),2) = 1;
%     megData.trialinfo(indtrial(D_raw,'EcHr','good'),2) = 1;
%     eegData.trialinfo(indtrial(D_raw,'EcHr','good'),2) = 1;
%     megData.trialinfo(indtrial(D_raw,'HcEr','good'),2) = 1;
%     eegData.trialinfo(indtrial(D_raw,'HcEr','good'),2) = 1;
%     megData.trialinfo(indtrial(D_raw,'HcHr','good'),2) = 1;
%     eegData.trialinfo(indtrial(D_raw,'HcHr','good'),2) = 1;

    data(counter).D_raw = D_raw; % grab the original data
    data(counter).D_clean_meg = D_clean_meg; % grab the data transformed by the ica
    data(counter).D_clean_eeg = D_clean_eeg; % grab the data transformed by the ica
    
    data(counter).EcEr.inds = indtrial(D_raw,'EcEr');
    data(counter).EcHr.inds = indtrial(D_raw,'EcHr');
    data(counter).HcEr.inds = indtrial(D_raw,'HcEr');
    data(counter).HcHr.inds = indtrial(D_raw,'HcHr');
    data(counter).EcEr_g.inds = indtrial(D_raw,'EcEr','good');
    data(counter).EcHr_g.inds = indtrial(D_raw,'EcHr','good');
    data(counter).HcEr_g.inds = indtrial(D_raw,'HcEr','good');
    data(counter).HcHr_g.inds = indtrial(D_raw,'HcHr','good');
    
    data(counter).eogs = D_raw.indchantype('EOG');
    data(counter).ecg = D_raw.indchantype('ECG');
    data(counter).eegs = D_raw.indchantype('EEG');
    data(counter).magnetos = D_raw.indchantype('MEGMAG');
    data(counter).planars = D_raw.indchantype('MEGPLANAR');
    
    
end



% megData = spm2fieldtrip(data(subject).D_clean_meg);
% cfg = [];
% cfg.trials = 1:size(data_meg.trialinfo,1);
% timeLockedData = ft_timelockanalysis(cfg,data_meg);
% 
% cfg = [];
% cfg.fontsize = 6;
% cfg.layout = 'neuromag306all.lay';
% figure
% ft_multiplotER(cfg,timeLockedData)
% 
% cfg = [];
% cfg.xlim = [0.5 1];
% ft_topoplotER(cfg,timeLockedData)

%% plot onsets from occipital sensors
disp('plotting onsets')

% 0 epoch start; 1500 dots onset; 3000 trial end; 4000 epoch end.
timePoints = [1000:3000];
channelLabels = {'MEG2041' 'MEG2111' 'MEG2141' 'MEG1911' ...
'MEG1921' 'MEG1931' 'MEG1941' 'MEG1731' 'MEG1741' 'MEG1641' 'MEG1721' ...
'MEG1711' 'MEG2031' 'MEG2121' 'MEG2131' ...
'MEG2341''MEG 2331' 'MEG2311' 'MEG2321' 'MEG2511' 'MEG2541' 'MEG2431' 'MEG2521' ...
'MEG2531'};
for subject = 1:numel(data)
    trialIndices = [data(subject).EcEr_g.inds,data(subject).EcHr_g.inds,data(subject).HcEr_g.inds,data(subject).HcHr_g.inds];
    dataObject = data(subject).D_clean_meg;
    [trialAverage, acrossSensors] = returnTrialData(dataObject,channelLabels,timePoints,trialIndices);
    plot(trialAverage)
    hold on
    markX([500 2000])
    hold off
end



%% plot motor channels leading to response (motor prep signal)

% 0 epoch start; 1500 dots onset; 3000 trial end; 4000 epoch end.
channelLabels = {'MEG0631' 'MEG0711' 'MEG0741' 'MEG1831' ...
'MEG2011' 'MEG0421' 'MEG0431' 'MEG1821' 'MEG1841' 'MEG0411' 'MEG0441' 'MEG1811' ...
'MEG1631' 'MEG1041' 'MEG0721' 'MEG0731' 'MEG2241' ...
'MEG2021' 'MEG1111' 'MEG1141' 'MEG2211' 'MEG2231' 'MEG1121' 'MEG1131' 'MEG2221' ...
'MEG2441'};
for subject = 1:numel(data)
      
    trialIndices = [data(subject).EcEr_g.inds,data(subject).EcHr_g.inds,data(subject).HcEr_g.inds,data(subject).HcHr_g.inds];
    
    % so we take a window around the response
    timeWindow = [-1000 500];
    % now we need to loop through the trials
    clear theseRts thisTimeWindow theseTrialData
    for trialNum = trialIndices
            % find the response time
            theseRts(trialNum) = str2num(data(subject).behav.megBehaviouralData(trialNum,10));
            % also need to add the variable random dots fixation
            theseRts(trialNum) = theseRts(trialNum) + (1000*data(subject).fixationTimings(trialNum));
            % and add the onset
            theseRts(trialNum) = theseRts(trialNum) + 1000;
            
            thisTimeWindow = timeWindow + theseRts(trialNum);
            
            timePoints = thisTimeWindow(1):thisTimeWindow(2);
            thisTrialIdx = trialNum;
            % get data
            theseTrialData(trialNum,:) = returnTrialData(data(subject).D_clean_meg,channelLabels,timePoints,thisTrialIdx);
            
    end
        % get rid of zeros - wherever this doesn't match the condition but also
        % when it does match the condition, and they didn't answer (properly)
        goodIdx = theseRts>0;
        theseTrialData = mean(theseTrialData(goodIdx,:));

    plot(theseTrialData)
    hold on
    markX(1000)
    hold off
end

%% topo from start of dots by 50ms
for subject = 1:numel(data)
    
    trialIndices = data(subject).EcEr_g.inds;
    %channelLabels = 'MEG';
    % now we need to loop through the trials
    clear theseRts thisOnset theseTrialData
    for trialNum = trialIndices
            % find  the variable random dots fixation
            thisOnset = 1000*data(subject).fixationTimings(trialNum);
            theseRts(trialNum) = thisOnset;
            % and add the onset
            theseRts(trialNum) = theseRts(trialNum) + 1000;
                       
            timePoints = theseRts(trialNum):theseRts(trialNum)+1000;
            thisTrialIdx = trialNum;
            % get data
            [~,theseTrialData(:,:,trialNum)] = returnTrialData(data(subject).D_clean_meg,[],timePoints,thisTrialIdx);
            
    end
        % get rid of zeros - wherever this doesn't match the condition but also
        % when it does match the condition, and they didn't answer (properly)
%         goodIdx = theseRts>0;
%         theseTrialData = theseTrialData(goodIdx,:);

        for thisTime = 2:50:numel(theseTrialData)
            topo=squeeze(mean(theseTrialData(:,thisTime,:),3));
            plotSensorTopographies(data(subject).D_clean_meg,topo,{'MEGMAG','MEGPLANAR'})
        end


end
%% take a baseline?

%% plot a timeseries with errorbars for the four conditions


% 0 epoch start; 1500 dots onset; 3000 trial end; 4000 epoch end.
channelLabels = {'MEG0631' 'MEG0711' 'MEG0741' 'MEG1831' ...
'MEG2011' 'MEG0421' 'MEG0431' 'MEG1821' 'MEG1841' 'MEG0411' 'MEG0441' 'MEG1811' ...
'MEG1631' 'MEG1041' 'MEG0721' 'MEG0731' 'MEG2241' ...
'MEG2021' 'MEG1111' 'MEG1141' 'MEG2211' 'MEG2231' 'MEG1121' 'MEG1131' 'MEG2221' ...
'MEG2441'};
for subject = 1:numel(data)
    
    trialIndices{1} = data(subject).EcEr_g.inds;
    trialIndices{2} = data(subject).EcHr_g.inds;
    trialIndices{3} = data(subject).HcEr_g.inds;
    trialIndices{4} = data(subject).HcHr_g.inds;
    
    % so we take a window around the response
    timeWindow = [-1000 500];
    for condNum = 1:numel(trialIndices)
        % now we need to loop through the trials
        clear theseRts thisTimeWindow theseTrialData
        for trialNum = trialIndices{condNum}
            % find the response time
            theseRts(trialNum) = str2num(data(subject).behav(trialNum,10));
            % also need to add the variable random dots fixation
            theseRts(trialNum) = theseRts(trialNum) + (1000*data(subject).fixationTimings(trialNum));
            % and add the onset
            theseRts(trialNum) = theseRts(trialNum) + 1000;
            
            thisTimeWindow = timeWindow + theseRts(trialNum);
            
            timePoints = thisTimeWindow(1):thisTimeWindow(2);
            thisTrialIdx = trialNum;
            % get data
            theseTrialData(trialNum,:) = returnTrialData(data(subject).D_clean_meg,channelLabels,timePoints,thisTrialIdx);
            
        end
        % get rid of zeros - wherever this doesn't match the condition but also
        % when it does match the condition, and they didn't answer (properly)
        goodIdx = theseRts>0;
        theseTrialData(~goodIdx,:) = [];
        allCondData{condNum} = theseTrialData;
        
    
    end
    
    for condNum = 1:numel(allCondData)
        [thisMean, thisSem] = returnStats(allCondData{condNum});
        plot(thisMean)
        hold on
        plot(thisMean+thisSem)
        plot(thisMean-thisSem)
        markX(1000)
        hold off
    end

end





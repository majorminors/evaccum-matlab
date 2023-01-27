function a5_epochExam(thisSubject,datadir,toolsdir,runLocal)

%% setup

rootDir = fullfile(datadir,thisSubject.id);
inputFolder = fullfile(rootDir,'Preprocess');
outputFolder = fullfile(inputFolder,'epochExam'); if ~exist(outputFolder,'dir'); mkdir(outputFolder); end
oslDir = fullfile('/imaging/woolgar/projects/Dorian/Toolboxes/osl'); % osl wants this if run in shared user mode
addpath(fullfile(oslDir,'osl-core')); % osl_startup can error without this
if runLocal; oslUserMode = 'user'; elseif ~runLocal; oslUserMode = 'shared';end

behavDataDir = fullfile(datadir,'behavioural');
behavData = fullfile(behavDataDir,[thisSubject.id '_EvAccum.mat']);
% and we'll output this revised version of it to play with
megrtFile = fullfile(behavDataDir,[thisSubject.id '_MEGRTs.mat']);
rawDataFile = sprintf([inputFolder,filesep,'allRuns_trans.mat']);

filterPrefix = 'f'; % most programs default to an f, but note we filter twice so will have a file with ff prepended
artefactPrefix = 'b';
icaPrefix = 'ica'; % we add this to the file ourselves
epochPrefix = 'e'; % we also choose this
filteredDataFile = addPrefix(rawDataFile,[filterPrefix filterPrefix]);
deArtefactedFile = addPrefix(filteredDataFile,artefactPrefix);
icaedFile = addPrefix(deArtefactedFile,icaPrefix);
epochedFile = addPrefix(icaedFile,epochPrefix);
deArtefactedEpochs = addPrefix(epochedFile,artefactPrefix);

disp('starting up osl')
osl_startup(oslDir,oslUserMode)

%% continuous data checks

disp('loading continuous data')
% can oslview any of the following
D_continuous_raw = spm_eeg_load(rawDataFile);
D_continuous_filtered = spm_eeg_load(filteredDataFile);
D_continuous_deArtefacted = spm_eeg_load(deArtefactedFile);
D_continuous_ica = spm_eeg_load(icaedFile);

% can check any of the above
% oslview(D_continuous_ica);

% let's grab some sensors
eogs = D_continuous_raw.indchantype('EOG');
ecg = D_continuous_raw.indchantype('ECG');
eegs = D_continuous_raw.indchantype('EEG');
magnetos = D_continuous_raw.indchantype('MEGMAG');
planars = D_continuous_raw.indchantype('MEGPLANAR');
D_epoched_deArtefacted = spm_eeg_load(deArtefactedEpochs); % grab this real quick, from the last preprocessing step
badChannels = D_epoched_deArtefacted.badchannels;
clear D_epoched_deArtefacted

% let's compare our continuous data against our artefacts---see that it's nice and cleaned up
has_montage(D_continuous_ica) % should show available montages
D_raw = D_continuous_ica.montage('switch',0); % grab the original data
D_clean_meg = D_continuous_ica.montage('switch',1); % grab the data transformed by the ica
D_clean_eeg = D_continuous_ica.montage('switch',2); % grab the data transformed by the ica

% should be sensor by time by trial
% ok, so this check shows me that the meg ica montage does not do anything
% to the eeg data, and the eeg ica montage does not do anything to the meg
% data - good to know!
% 
% rawEegData = mean(D_raw([eegs(:)],:));
% rawMegmagData = mean(D_raw([magnetos(:)],:));
% rawMegplanarData = mean(D_raw([planars(:)],:));
% 
% DmegEegData = mean(D_clean_meg([eegs(:)],:));
% DmegMegmagData = mean(D_clean_meg([magnetos(:)],:));
% DmegMegplanarData = mean(D_clean_meg([planars(:)],:));
% 
% DeegEegData = mean(D_clean_eeg([eegs(:)],:));
% DeegMegmagData = mean(D_clean_eeg([magnetos(:)],:));
% DeegMegplanarData = mean(D_clean_eeg([planars(:)],:));
% 
% disp(['raw eeg = MEG ica? ' num2str(sum(rawEegData == DmegEegData)/D_raw.nsamples)])
% disp(['raw eeg = EEG ica? ' num2str(sum(rawEegData == DeegEegData))])
% 
% disp(['raw megmag = MEG ica? ' num2str(sum(rawMegmagData == DeegMegmagData)/D_raw.nsamples)])
% disp(['raw megmag = EEG ica? ' num2str(sum(rawMegmagData == DmegMegmagData))])
% 
% disp(['raw megplanar = MEG ica? ' num2str(sum(rawMegplanarData == DeegMegplanarData)/D_raw.nsamples)])
% disp(['raw megplanar = EEG ica? ' num2str(sum(rawMegplanarData == DmegMegplanarData))])
% 
% we can also plot it
% 
% subplot(3,1,1);
% plot(rawEegData,'m');
% xlabel('time','FontSize',15);
% set(gca,'FontSize',15)
% ylabel('mean eeg sensor vals','FontSize',15);
% hold on
% plot(DmegEegData,'c');
% plot(DeegEegData,'b');
% 
% subplot(3,1,2);
% plot(rawMegmagData,'m');
% xlabel('time','FontSize',15);
% set(gca,'FontSize',15)
% ylabel('mean megmag sensor vals','FontSize',15);
% hold on
% plot(DmegMegmagData,'c');
% plot(DeegMegmagData,'b');
% 
% subplot(3,1,3);
% plot(rawMegplanarData,'m');
% xlabel('time','FontSize',15);
% set(gca,'FontSize',15)
% ylabel('mean megplanar sensor vals','FontSize',15);
% hold on
% plot(DmegMegplanarData,'c');
% plot(DeegMegplanarData,'b');

% now we can check that we've successfully cleaned the data
% make a little construct to loop through
first = @(x) x(1); getProblemSensors = @(allSensorInds,chanLabels) setdiff(allSensorInds(find(ismember(D_continuous_raw.chanlabels(allSensorInds),chanLabels))),badChannels);

artefactChecks(1).artefactSensor = eogs(2);
artefactChecks(1).sensorToCheck = first(getProblemSensors(eegs,{'EEG002' 'EEG004' 'EEG005' 'EEG006' 'EEG007' 'EEG008' 'EEG009' 'EEG010' 'EEG011' 'EEG012' 'EEG013' 'EEG014' 'EEG015' 'EEG016' 'EEG017' 'EEG018'}));
artefactChecks(1).cleanedData = D_clean_eeg;
artefactChecks(1).channelTitle = 'eeg';
artefactChecks(2).artefactSensor = eogs(2);
artefactChecks(2).sensorToCheck = first(getProblemSensors(magnetos,{'MEG0521' 'MEG0811' 'MEG0611' 'MEG0911' 'MEG0311' 'MEG0921' 'MEG0641' 'MEG531' 'MEG0821' 'MEG0941' 'MEG0931'}));
artefactChecks(2).cleanedData = D_clean_meg;
artefactChecks(2).channelTitle = 'megmag';
artefactChecks(2).sensorToCheck = first(getProblemSensors(planars,{'MEG0522' 'MEG0523' 'MEG0812' 'MEG0813' 'MEG0612' 'MEG0613' 'MEG0912' 'MEG0913' 'MEG0312' 'MEG0313' 'MEG0922' 'MEG0923' 'MEG0642' 'MEG0643' 'MEG532' 'MEG533' 'MEG0822' 'MEG0823' 'MEG0942' 'MEG0943' 'MEG0932' 'MEG0933'}));
artefactChecks(2).cleanedData = D_clean_meg;
artefactChecks(2).channelTitle = 'megplanars';

for figIdx = 1:numel(artefactChecks)
    thisArtefactSensor = artefactChecks(figIdx).artefactSensor;
    thisSensorToCheck = artefactChecks(figIdx).sensorToCheck;
    thisCleanedData = artefactChecks(figIdx).cleanedData;
    thisChannelTitle = artefactChecks(figIdx).channelTitle;
    
    figureHandle = plotAgainstArtefact(thisArtefactSensor,thisSensorToCheck,D_raw,thisCleanedData,thisChannelTitle);
    save_as_png(figureHandle, [outputFolder filesep 'artefact_check_' thisChannelTitle '.png'], 150, 1200, 800)
    close all
    clear figureHandle
end

clear D_raw D_clean_eeg D_clean_meg

%% epoched data checks

disp('loading epoched data')
% cannot oslview these
D_epoched_ica = spm_eeg_load(epochedFile);
D_epoched_deArtefacted = spm_eeg_load(deArtefactedEpochs);

has_montage(D_epoched_deArtefacted) % should show available montages

trl = load([inputFolder,filesep,'allRuns_trl.mat']); % if we want to check this, for some reason

D_raw = D_epoched_deArtefacted.montage('switch',0); % grab the original data
D_clean_meg = D_epoched_deArtefacted.montage('switch',1); % grab the data transformed by the ica
D_clean_eeg = D_epoched_deArtefacted.montage('switch',2); % grab the data transformed by the ica

% should be sensor by time by trial
% response = @(sensor,trial) plot(D.time,mean(D(sensor,:,trial),3));
% average_response = @(sensor) plot(D.time,mean(D(sensor,:,:),3));


EcEr = indtrial(D_raw,'EcEr');
EcHr = indtrial(D_raw,'EcHr');
HcEr = indtrial(D_raw,'HcEr');
HcHr = indtrial(D_raw,'HcHr');
EcEr_g = indtrial(D_raw,'EcEr','good');
EcHr_g = indtrial(D_raw,'EcHr','good');
HcEr_g = indtrial(D_raw,'HcEr','good');
HcHr_g = indtrial(D_raw,'HcHr','good');

% plotting average of all sensors for raw and cleaned as a 2D
% matrix for comparison and also to see what sensors are doing what
subplot(2,3,1);imagesc(D_raw.time,[],squeeze(mean(D_raw([planars(:)],:,EcEr_g),3)));
xlabel('Time (seconds)','FontSize',20);
ylabel('Sensors','FontSize',15);colorbar
title('raw megplanars','FontSize',15)
set(gca,'FontSize',15)

subplot(2,3,2); % plots magnetometers
imagesc(D_raw.time,[],squeeze(mean(D_raw([magnetos(:)],:,EcEr_g),3)));
xlabel('Time (seconds)','FontSize',15);
ylabel('Sensors','FontSize',15);colorbar
title('raw megmags','FontSize',15)
set(gca,'FontSize',15)

subplot(2,3,3); % plots eeg
imagesc(D_raw.time,[],squeeze(mean(D_raw([eegs(:)],:,EcEr_g),3)));
xlabel('Time (seconds)','FontSize',15);
ylabel('Sensors','FontSize',15);colorbar
title('raw eegs','FontSize',15)
set(gca,'FontSize',15)

subplot(2,3,4);imagesc(D_clean_meg.time,[],squeeze(mean(D_clean_meg([planars(:)],:,EcEr_g),3)));
xlabel('Time (seconds)','FontSize',20);
ylabel('Sensors','FontSize',15);colorbar
title('clean megplanars','FontSize',15)
set(gca,'FontSize',15)

subplot(2,3,5); % plots magnetometers
imagesc(D_clean_meg.time,[],squeeze(mean(D_clean_meg([magnetos(:)],:,EcEr_g),3)));
xlabel('Time (seconds)','FontSize',15);
ylabel('Sensors','FontSize',15);colorbar
title('clean megmags','FontSize',15)
set(gca,'FontSize',15)

subplot(2,3,6); % plots eeg
imagesc(D_clean_eeg.time,[],squeeze(mean(D_clean_eeg([eegs(:)],:,EcEr_g),3)));
xlabel('Time (seconds)','FontSize',15);
ylabel('Sensors','FontSize',15);colorbar
title('clean eegs','FontSize',15)
set(gca,'FontSize',15)

save_as_png(gcf, [outputFolder filesep 'epoch_average_sensors_clean_v_raw.png'], 150, 1200, 800)

close all

% epoch plotting of occipital/parietal sensors (not all, or signal will be
% gone
% 1. subject average (evoked response) locked to onset of random dots (so
% epoch middle minus p.fixation_dots_duration?) and see a massive onset!
% 2. do the same for motion onset
% 3. do the sa me for response (so add back RT) - should be parietal
% then do the average across participants

%% let's look at the onsets

for cond = [1 2 3 4]
    for sens = [1 2 3]
        
        switch sens
            case 1
                theseSensors = eegs;
                thisClnData = D_clean_eeg;
                sensorType = 'eegs';
            case 2
                theseSensors = magnetos;
                thisClnData = D_clean_meg;
                sensorType = 'megmags';
            case 3
                theseSensors = planars;
                thisClnData = D_clean_meg;
                sensorType = 'megmags';
        end
        
        switch cond
            case 1
                trialIndices = EcEr_g;
                condName = 'EcEr';
            case 2
                trialIndices = EcHr_g;
                condName = 'EcHr';
            case 3
                trialIndices = HcEr_g;
                condName = 'HcEr';
            case 4
                trialIndices = HcHr_g;
                condName = 'HcHr';
        end
        % should be sensor by time by trial
        rawTrls = mean(mean(D_raw([theseSensors(:)],:,trialIndices),3));
        clnTrls = mean(mean(thisClnData([theseSensors(:)],:,trialIndices),3));
        
        plot(rawTrls(1000:end),'r')
        title(['average ' sensorType ' onset for ' condName])
        hold on
        plot(clnTrls(1000:end))
        plot([1500 1500],[ylim],'g'); % mark the end of the trial
        save_as_png(gcf, [outputFolder filesep 'average_onset_' sensorType '_' condName '.png'], 150, 800, 800)
        close all
    end
    for timepoint = [1000 1500 1700 2000]
        topo=squeeze(mean(thisClnData(:,timepoint,trialIndices),3));
        cfg = [];
        cfg.output = [inputFolder filesep 'eeg_layout.lay'];
        cfg.elec = D_raw.sensors('EEG');
        eeg_layout = ft_prepare_layout(cfg);
        plotSensorTopographies(thisClnData,topo,{'MEGMAG','MEGPLANAR','EEG'},1,[inputFolder filesep 'eeg_layout.lay'],cfg.elec.label)
        save_as_png(gcf, [outputFolder filesep 'average_onset_topo_' num2str(timepoint-1000) '_' condName '.png'], 150, 1200, 800)
        close all
    end
end

%% now responses --- much harder

% so we take a window around the response
timeWindow = [-1000 500];
load(behavData,'p')
fixationDots = p.fixation_dots_duration_vector;
load(megrtFile,'megBehaviouralData')
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
% 12)  we later add a row of unique numbers to tag each trial

for sens = [1 2 3]
    
    switch sens
        case 1
            theseSensors = eegs;
            thisClnData = D_clean_eeg;
            sensorType = 'eegs';
        case 2
            theseSensors = magnetos;
            thisClnData = D_clean_meg;
            sensorType = 'megmags';
        case 3
            theseSensors = planars;
            thisClnData = D_clean_meg;
            sensorType = 'megplanars';
    end
    
    for cond = [1 2 3 4]
        
        switch cond
            case 1
                trialIndices = EcEr_g;
                condName = 'EcEr';
                colour = 'c';
            case 2
                trialIndices = EcHr_g;
                condName = 'EcHr';
                colour = 'b';
            case 3
                trialIndices = HcEr_g;
                condName = 'HcEr';
                colour = 'm';
            case 4
                trialIndices = HcHr_g;
                condName = 'HcHr';
                colour = 'r';
        end
        
        % now we need to loop through the trials
        clear theseRts thisTimeWindow thisSensorDataRaw thisSensorDataClean
        for trialNum = trialIndices
            % find the response time
            theseRts(trialNum) = str2num(megBehaviouralData(trialNum,10));
            % also need to add the variable random dots fixation
            theseRts(trialNum) = theseRts(trialNum) + (1000*fixationDots(trialNum));
            % and add the onset
            theseRts(trialNum) = theseRts(trialNum) + 1000;
            
            thisTimeWindow = timeWindow + theseRts(trialNum);
            % should be sensor by time by trial
            thisSensorDataRaw(trialNum,:) = mean(D_raw([theseSensors(:)],thisTimeWindow(1):thisTimeWindow(2),trialNum));
            thisSensorDataClean(trialNum,:) = mean(thisClnData([theseSensors(:)],thisTimeWindow(1):thisTimeWindow(2),trialNum));
            allSensorRespOnsets(:,trialNum) = thisClnData(:,thisTimeWindow(1),trialNum);
        end; %clear sensorDataRaw sensorDataClean
        if sens == 1
            % plot all sensor topographies for the onset
            topo = mean(allSensorRespOnsets,2);
            plotSensorTopographies(thisClnData,topo,{'MEGMAG','MEGPLANAR','EEG'},1,[inputFolder filesep 'eeg_layout.lay'],cfg.elec.label)
            save_as_png(gcf, [outputFolder filesep 'average_response_onset_topo_' condName '.png'], 150, 1200, 800)
            close all
        end
        % get rid of zeros - wherever this doesn't match the condition but also
        % when it does match the condition, and they didn't answer (properly)
        goodIdx = theseRts>0;
        sensorData(cond).raw = thisSensorDataRaw(goodIdx,:);
        sensorData(cond).clean = thisSensorDataClean(goodIdx,:);
        plot(smoothdata(mean(sensorData(cond).clean),'gaussian',50),colour)
        hold on
    end
    title(['smoothed ' sensorType ' average responses'])
    plot([1000 1000],[ylim],'g'); % mark the response
    save_as_png(gcf, [outputFolder filesep 'smoothed_average_responses_' sensorType '.png'], 150, 800, 800)
    close all
end

getSensors = @(allSensorInds,chanLabels) allSensorInds(find(ismember(D_continuous_raw.chanlabels(allSensorInds),chanLabels)));

% getSensors(eegs,{})
for sens = [1 2 3]
    
    switch sens
        case 1
            theseSensors = getSensors(eegs,{'EEG032' 'EEG033' 'EEG034' 'EEG035' 'EEG036' 'EEG043' 'EEG044' 'EEG045' 'EEG046' 'EEG047'});
            thisClnData = D_clean_eeg;
            sensorType = 'eegs';
        case 2
            theseSensors = getSensors(magnetos,{'MEG0631' 'MEG0711' 'MEG0741' 'MEG1831' 'MEG2011' 'MEG1041' 'MEG0721' 'MEG0731' 'MEG2241' 'MEG2021'});
            thisClnData = D_clean_meg;
            sensorType = 'megmags';
        case 3
            theseSensors = getSensors(planars,{'MEG0632' 'MEG0633' 'MEG0712' 'MEG0713' 'MEG0742' 'MEG0743' 'MEG1832' 'MEG1833' 'MEG2012' 'MEG2013' 'MEG1041' 'MEG1042' 'MEG1043' 'MEG0722' 'MEG0723' 'MEG0732' 'MEG0733' 'MEG2242' 'MEG2243' 'MEG2022' 'MEG2023'});
            thisClnData = D_clean_meg;
            sensorType = 'megplanars';
    end
    
    for cond = [1 2 3 4]
        
        switch cond
            case 1
                trialIndices = EcEr_g;
                condName = 'EcEr';
                colour = 'c';
            case 2
                trialIndices = EcHr_g;
                condName = 'EcHr';
                colour = 'b';
            case 3
                trialIndices = HcEr_g;
                condName = 'HcEr';
                colour = 'm';
            case 4
                trialIndices = HcHr_g;
                condName = 'HcHr';
                colour = 'r';
        end
        
        % now we need to loop through the trials
        clear theseRts thisTimeWindow thisSensorDataRaw thisSensorDataClean
        for trialNum = trialIndices
            % find the response time
            theseRts(trialNum) = str2num(megBehaviouralData(trialNum,10));
            % also need to add the variable random dots fixation
            theseRts(trialNum) = theseRts(trialNum) + (1000*fixationDots(trialNum));
            % and add the onset
            theseRts(trialNum) = theseRts(trialNum) + 1000;
            
            thisTimeWindow = timeWindow + theseRts(trialNum);
            % should be sensor by time by trial
            thisSensorDataRaw(trialNum,:) = mean(D_raw([theseSensors(:)],thisTimeWindow(1):thisTimeWindow(2),trialNum));
            thisSensorDataClean(trialNum,:) = mean(thisClnData([theseSensors(:)],thisTimeWindow(1):thisTimeWindow(2),trialNum));
            
        end; %clear sensorDataRaw sensorDataClean
        % get rid of zeros - wherever this doesn't match the condition but also
        % when it does match the condition, and they didn't answer (properly)
        goodIdx = theseRts>0;
        sensorData(cond).raw = thisSensorDataRaw(goodIdx,:);
        sensorData(cond).clean = thisSensorDataClean(goodIdx,:);
        
        plot(smoothdata(mean(sensorData(cond).clean),'gaussian',50),colour)
        hold on
    end
    title(['smoothed ' sensorType ' average midline parietal responses'])
    plot([1000 1000],[ylim],'g'); % mark the response
    save_as_png(gcf, [outputFolder filesep 'smoothed_average_midline_parietal_responses_' sensorType '.png'], 150, 800, 800)
    close all
end




osl_shutdown



return
end

function figureHandle = plotAgainstArtefact(artefactSensor,sensorToCheck,noisyData,cleanedData,channelTitle)

% artefactSensor: this should be an ECG or EOG sensor - something that is artefactual
% sensorToCheck:  this should be a sensor near your artefact sensor - we will compare it to the artefact before and after ICA
% noisyData:      an spm m/eeg object that the artefact will probably be affecting
% cleanedData:    an spm m/eeg object that the artefact will probably cleaned from
%                 be sure to switch to the appropriate montage, (e.g.
%                 D.montage('switch',1)`) before sending it in! if you have the raw
%                 data as montage 0, then the noisy and clean data won't be visible
% channelTitle:   what to insert before the word 'channel' in the figure titles

figure;
subplot(2,1,1);
plot(noisyData.time(1:10000),noisyData(artefactSensor,1:10000)); % takes first 10000 sample points
title('artefact channel')
xlim([0 10]);
xlabel('Time (s)');

subplot(2,1,2);
plot(noisyData.time(1:10000),noisyData(sensorToCheck,1:10000)); % takes first 10000 sample points
title(['artefact contaminated ' channelTitle ' channel'])
xlim([0 10]);
hold on;
plot(cleanedData.time(1:10000),cleanedData(sensorToCheck,1:10000),'r');
xlim([0 10]);
xlabel('Time (s)');
legend({'pre ICA' 'post ICA'});

figureHandle = gcf;


return
end






function outFiles = addPrefix(files,prefix)

if iscell(files)
    for fileIdx = 1:length(files)
        [pathstr,name,ext] = fileparts(files{fileIdx});
        outFiles{fileIdx} = sprintf(['%s/' prefix '%s%s'],pathstr,name,ext);
    end
else
    [pathstr,name,ext] = fileparts(files);
    outFiles = sprintf(['%s/' prefix '%s%s'],pathstr,name,ext);
end

return
end

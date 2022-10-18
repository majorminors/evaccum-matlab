clear all

%epoch exam at group level

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';
scriptdir = fullfile(rootdir,'tools_analysis'); cd(scriptdir)
datadir = fullfile(rootdir,'data','meg_pilot_4'); addpath(datadir);
behavDataDir = fullfile(datadir,'behavioural');
outputFolder = fullfile(datadir,'epochExam'); if ~exist(outputFolder,'dir'); mkdir(outputFolder); end

oslDir = fullfile('/imaging/woolgar/projects/Dorian/Toolboxes/osl'); % osl wants this if run in shared user mode
addpath(fullfile(oslDir,'osl-core')); % osl_startup can error without this

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
    data(counter).D_raw = D_epoched_deArtefacted.montage('switch',0); % grab the original data
    data(counter).D_clean_meg = D_epoched_deArtefacted.montage('switch',1); % grab the data transformed by the ica
    data(counter).D_clean_eeg = D_epoched_deArtefacted.montage('switch',2); % grab the data transformed by the ica
    
    D_raw = data(counter).D_raw;
    D_clean_meg = data(counter).D_clean_meg;
    D_clean_eeg = data(counter).D_clean_eeg;
    
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
    
    timeWindow = [-1000 500];
    load(behavData,'p')
    fixationDots = p.fixation_dots_duration_vector;
    load(megrtFile,'megBehaviouralData')
    
    disp('collating data for this subject')
    for sens = [1 2 3]
        fprintf('sensor %.0f of 3\n',sens)
        
        switch sens
            case 1
                theseSensors = data(counter).eegs;
                thisClnData = D_clean_eeg;
                sensorType = 'eegs';
            case 2
                theseSensors = data(counter).magnetos;
                thisClnData = D_clean_meg;
                sensorType = 'megmags';
            case 3
                theseSensors = data(counter).planars;
                thisClnData = D_clean_meg;
                sensorType = 'megplanars';
        end
        
        for cond = [1 2 3 4]
            fprintf('condition %.0f of 4\n',cond)
            
            switch cond
                case 1
                    trialIndices = data(counter).EcEr_g.inds;
                    condName = 'EcEr_g';
                    colour = 'c';
                case 2
                    trialIndices = data(counter).EcHr_g.inds;
                    condName = 'EcHr_g';
                    colour = 'b';
                case 3
                    trialIndices = data(counter).HcEr_g.inds;
                    condName = 'HcEr_g';
                    colour = 'm';
                case 4
                    trialIndices = data(counter).HcHr_g.inds;
                    condName = 'HcHr_g';
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

            goodIdx = theseRts>0;
            thisSensorDataRaw = thisSensorDataRaw(goodIdx,:);
            thisSensorDataClean = thisSensorDataClean(goodIdx,:);
            if ~exist('averageCleanOnsets','var')
                averageCleanOnsets = mean(thisSensorDataClean);
            else
                averageCleanOnsets = mean(cat(1,averageCleanOnsets,mean(thisSensorDataClean)));
            end
            if ~exist('averageCleanRaw','var')
                averageRawOnsets = mean(thisSensorDataRaw);
            else
                averageRawOnsets = mean(cat(1,averageRawOnsets,mean(thisSensorDataRaw)));
            end
            data(counter).(condName).allSensorRespOnsets = allSensorRespOnsets;
            % get rid of zeros - wherever this doesn't match the condition but also
            % when it does match the condition, and they didn't answer (properly)
            data(counter).(condName).sensorData.raw = thisSensorDataRaw;
            data(counter).(condName).sensorData.clean = thisSensorDataClean;
            data(counter).(condName).conditionCode = condName;
            
            aves.(sensorType).(condName).averageRawOnsets = averageRawOnsets;
            aves.(sensorType).(condName).averageCleanOnsets = averageCleanOnsets;
        end
        if sens == 1
            
            if ~exist('averageAllSensorRespOnsets','var')
                averageAllSensorRespOnsets = mean(allSensorRespOnsets,2);
            else
                averageAllSensorRespOnsets = mean(cat(2,mean(allSensorRespOnsets,2),averageAllSensorRespOnsets),2);
            end
            aves.(condName).averageAllSensorRespOnsets = averageAllSensorRespOnsets;
            
        end
    end
    
    
    
end

for sens = {'eegs' 'megmags' 'megplanars'}
    sens = sens{:}
    for cond = {'EcEr_g' 'EcHr_g' 'HcEr_g' 'HcHr_g'}
        cond = cond{:};
        switch cond
            case 'EcEr_g'
                colour = 'c';
            case 'EcHr_g'
                colour = 'b';
            case 'HcEr_g'
                colour = 'm';
            case 'HcHr_g'
                colour = 'r';
        end
        
        %         plot(aves.(sens).(cond).averageRawOnsets);
        plot(smoothdata(aves.(sens).(cond).averageCleanOnsets,'gaussian',50),colour);
        hold on
    end
    title(['average ' sens ' onset for all conds'])
    plot([1000 1000],[ylim],'g'); % mark the response
    save_as_png(gcf, [outputFolder filesep 'average_onset_' sens '.png'], 150, 800, 800)
    close all
end



osl_shutdown

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

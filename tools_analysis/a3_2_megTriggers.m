function a3_2_megTriggers(thisSubject,datadir,toolsdir,overwrite)

%addpath /hpc-software/matlab/cbu/
addpath([toolsdir '/fieldtrip-lite-20190410']);

%% directory mapping

megdata = fullfile(datadir,thisSubject.id);
trialEpochFile  = fullfile(megdata, 'Preprocess','run%s_trl.mat'); % file to save trial epoch info
origdatadir = fullfile('/megdata/cbu/evaccum',thisSubject.meg_fld,thisSubject.date_meg);
% take my input behavioural data
behavDataDir = fullfile(datadir,'behavioural');
behavData = fullfile(behavDataDir,[thisSubject.id '_EvAccum.mat']);
% and we'll output this revised version of it to play with
megrtFile = fullfile(behavDataDir,[thisSubject.id '_MEGRTs.mat']);
addpath(genpath(toolsdir))

if ~exist(megrtFile,'file') || overwrite
    

%% load the behavioural data

load(behavData) % matlab behavioural file
runid = thisSubject.runid; % pull information related to the blocks in each run
% runid looks like:
%     runid =
%
%      1     2   - is the run, and which blocks
%      3     4   - is run 2 with those blocks
%      5     6   - etc
%      7     8
%      9    10
% we'll use conditions to code later (though I guess we should really pull
% this from the behavioural file

%% a bit of setup
conditions = {'EcEr', 'EcHr','HcEr','HcHr'}; % 2x2 coherence and rule
megBehaviouralData = []; % init this for later
%init these as part of Ale's preprocessing pipeline
alltrl = [];
tmplabs ={};

% set up to use fieldtrip
ft_defaults;
cfg = [];

%% let's identify what trigger channels we care about for our function below

channels.trialTriggers.labels = {'STI001'};
channels.fixationTriggers.labels = {'STI002'};
channels.photoDiode.labels = {'STI010'};
% we could surely improve this. this might actually be a case for looking
% in ST101 and grabbing the exact values, because I think that actually
% works properly
channels.responseTriggersLeft.labels = {'STI013'};
channels.responseTriggersRight.labels = {'STI016'};

%% and now loop through our runs
for runi = 1:numel(thisSubject.meg_labs)
    
    %% first lets retrieve the MEG data in an intelligible format
    
    % so we can pull things in by running data = ft_preprocessing(cfg) on cfg.dataset
    % (the path to the raw data)
    % look at the data.label vector to see where your channels are coming in,
    % so lets put together the path for the raw data
    %     rawData = fullfile(origdatadir,[thisSubject.meg_runs{runi} '.fif']);
    rawData = fullfile(megdata,[thisSubject.meg_labs{runi} '_raw.fif']);
    cfg.dataset = rawData; % into the cfg structure for ft
    thisData = ft_preprocessing(cfg); % just collect the data
    
    % this is where we send off our channel labels of interest to get our
    % trigger content
    triggerInfo = getTriggersFromChannelLabels(thisData,channels);
    
    %% now let's pull together some RTs from the MEG data
    
    % lets grab our trigger onset times in nice readable variables that
    % will work well with matlabs IDE
    switch thisSubject.id
        case 'S05'
            % photodiode is crazy for s05---lots of extra triggers
            trialOnsets = triggerInfo.trialTriggers.onsetTimes;
            if runi == 2
                skipCorr = 1; % is >90
            end
        case 'S23'
            if runi == 3
                skipCorr = 1; % is ~78...
            end
            trialOnsets = triggerInfo.photoDiode.onsetTimes;
        otherwise
            trialOnsets = triggerInfo.photoDiode.onsetTimes;
    end
    allResponseOnsets = sort([triggerInfo.responseTriggersLeft.onsetTimes;triggerInfo.responseTriggersRight.onsetTimes]);
        
    

    
    % now we go through each trial onset trigger
    for onsetidx = 1:numel(trialOnsets)
        
        % if we're on any trial other than the last trial
        if onsetidx < numel(trialOnsets)
            
            % let's get the trigger timepoints between this trialOnset and the
            % next trialOnset
            trialTimepoints = trialOnsets(onsetidx):(trialOnsets(onsetidx+1)-1);
            % look for the trial timepoints that correspond to any response onset in that timeframe
            theseResponses = trialTimepoints(find(ismember(trialTimepoints,allResponseOnsets)));
            
        else
            
            % or, on the last trial, just the get the last trialOnset
            % timepoint
            trialTimepoints = trialOnsets(onsetidx);
            % and get any responses after that
            theseResponses = allResponseOnsets(allResponseOnsets > trialTimepoints);
            
        end
        
        % let's collate our onset times and calculate an rt
        if ~isempty(theseResponses) % if there is a response (not empty)
            responseTimes(onsetidx) = theseResponses(1); % get the first response timepoint (in case of multiple button presses)
            megRTs(onsetidx) = theseResponses(1)-trialOnsets(onsetidx);
        else % if there is no response
            responseTimes(onsetidx) = 0;
            megRTs(onsetidx) = 0;
        end
        
        
    end
    
    %% now we want to do any postprocessing of our rts
    
    % if we weren't using the photodiode to get accurate presentation times,
    % we'd want to add something (e.g. 34) to our RT,
    % to account for the average projector delay (e.g. @CBU = 34ms) so e.g.:
%     megRTs = megRTs+34;

    % I do want to scrub any rts that are longer than my trials here,
    % because the way I've coded my experiment, trial onsets can straddle a
    % cue, and participants press a button to respond to a cue, so if they
    % missed a response in the trial, but responded to a cue, we would pick
    % that up as an RT
    megRTs(megRTs>1500) = 0;
    responseTimes(megRTs>1500) = 0;
    

    
    %% lets also do some sanity checks here against our behavioural data
    
    % pull the information for the blocks in the run
    vecrunid = runid(runi,:); % select vector of blocks in the run
    vecrunid(isnan(vecrunid)) = []; % get rid of any NaNs if they exist
    block = d.stim_mat_all(:,:,vecrunid); % pull the stimulus matrix info for those blocks
    block = permute(block,[1,3,2]); % swaps 3rd and 2nd dimensions
    block = reshape(block,[],size(d.stim_mat_all,2),1); % reshapes so rows are autocalculated, columns defined by stimulus matrix columns, and forces 3rd dimension into the 2d structure
    
    % first we check we have the right number of triggers
    if size(block,1) ~= numel(trialOnsets)
        warning('the trial onsets and block-expected trials dont match')
        disp('lets do some subject specific editing')
        edited = 0;
        % we need to edit
        % trialOnsets
        % megRTs
        % responseTimes
        
        switch thisSubject.id
            case 'S03'
                % on run 4, there are a bunch of additional photodiodes
                if sum(megRTs(size(block,1)+1:end)) == 0
                    % so delete any megrts that happen after the block
                    megRTs = megRTs(1:size(block,1));
                    trialOnsets = trialOnsets(1:size(block,1));
                    responseTimes = responseTimes(1:size(block,1));
                    edited = 1;
                end
            case 'S10'
                if runi == 1
                    megRTs = megRTs(219:end);
                    trialOnsets = trialOnsets(219:end);
                    responseTimes = responseTimes(219:end);
                    edited = 1;
                end
            case 'S11'
                if runi == 4
                    megRTs = megRTs(1:size(block,1));
                    trialOnsets = trialOnsets(1:size(block,1));
                    responseTimes = responseTimes(1:size(block,1));
                    edited = 1; 
                end
            case 'S15'
                if runi == 2
                    % we missed the first part of the block, so make all
                    % those zeros
                    tmp = size(block,1) - numel(trialOnsets);
                    d.rt(vecrunid,1:tmp) = 0;
                    megRTs = [zeros(1,tmp),megRTs];
                    trialOnsets = [zeros(1,tmp)';trialOnsets];
                    responseTimes = [zeros(1,tmp),responseTimes];
                    edited = 1;
                    skipCorr = 1;
                elseif runi == 5
                    megRTs = megRTs(1:size(block,1));
                    trialOnsets = trialOnsets(1:size(block,1));
                    responseTimes = responseTimes(1:size(block,1));
                    edited = 1;
                end
            case 'S29'
                if runi == 4
                    megRTs = megRTs(1:size(block,1));
                    trialOnsets = trialOnsets(1:size(block,1));
                    responseTimes = responseTimes(1:size(block,1));
                    edited = 1;
                end
            case 'S32'
                if runi == 1
                    megRTs = megRTs(1:size(block,1));
                    trialOnsets = trialOnsets(1:size(block,1));
                    responseTimes = responseTimes(1:size(block,1));
                    edited = 1;
                end
            case 'S36'
                if runi == 4
                    megRTs = megRTs(1:size(block,1));
                    trialOnsets = trialOnsets(1:size(block,1));
                    responseTimes = responseTimes(1:size(block,1));
                    edited = 1;
                end
        end
        
        if ~edited
            error('Inconsistent number of onset triggers vs trials')
        end
    end
    
    if ~exist('skipCorr','var')
        skipCorr = 0;
    end

    behavRTs = d.rt(vecrunid,:)';behavRTs = behavRTs(:); % transpose relevant rows of d.rt---my behavioural rts---and then place each block in this run under one another in the same row
    
    if ~skipCorr
        % I'll also check that my megrts correlate well with my matlab rts
        RTidx = behavRTs~=0 & behavRTs>0; % get an index of where there should be an rt
        corrArray = [behavRTs(RTidx) transpose(megRTs(RTidx))]; % concatenate the two rt vectors
        corrArray = corr(corrArray); % correlate them
        if corrArray(2,1) <0.95 || isnan(corrArray(2,1)) % and error if something seems wrong
            error('large discrepancy between stim RT and MEG RT!')
        end
    end

    % ok, that should do us
    
    %% so now let's combine information we care about from our behavioural data with our MEG-obtained RTs
    % we will save this into a new file, so we don't fuck up our
    % behavioural data
    
    % so to re-compose matrix of behavioural data first we need to grab
    % stuff
    % so block looks like:
    %  1)  cue direction (1-4)
    %  2)  cue direction in degrees
    %  3)  dot motion direction condition (1-8)
    %  4)  dot motion direction in degrees
    %  5)  coherence difficulty (1 = easy, 2 = hard)
    %  6)  matching distance from cue direction (absolute value) - used to calc match and match difficulty
    %  7)  match blue (1) or match orange (2)
    %  8)  matching difficulty (1 = easy, 2 = difficult)
    %  9)  gives you a unique number for each trial condition
    % 10)  gives a number based on 9 for meg triggers
    for i = 1:length(block) % create a row that gives you a number for each condition in your 2x2
        if block(i,5)==1 && block(i,8)==1
            conditionlabels(i,1) = string(conditions{1});
        elseif block(i,5)==1 && block(i,8)==2
            conditionlabels(i,1) = string(conditions{2});
        elseif block(i,5)==2 && block(i,8)==1
            conditionlabels(i,1) = string(conditions{3});
        elseif block(i,5)==2 && block(i,8)==2
            conditionlabels(i,1) = string(conditions{4});
        end
    end
    
    accuracyrow = d.correct(vecrunid,:)'; accuracyrow = accuracyrow(:); % may as well put this in the new matrix - so transpose relevant rows and then place them one after another in a column
    
    %MEG_RT = cat(1,MEG_RT,[block(:,[1:4 12 14 15 ]), tmp_RT']);% 1block# 2coh 3act 4cond 5trigger 6resKey 7acc 8rt % compiles matlab trial information with the new MEG generated rts
    run_tagger = runi.*ones(length(block(:,1)),1);
    megBehaviouralData = cat(1,megBehaviouralData,[block(:,[1 3 5 8 9 10]), conditionlabels, accuracyrow, behavRTs, megRTs',run_tagger]);
    
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
    
    %% do some trial definition epoching
    % also messy, but needed in our pipeline
    
    % set the trial window, including buffer for edge effects during filtering
    % - in the textbook - lowest frequency will determine
    % time before onset for a baseline, and some time before epoch of
    % interest to do time frequency decomposition---padding to avoid edge
    % effects in freqs like alphga and beta
    trialWindow = [-1500 2500];%in samples at 1000Hz sampling % what time window 1500 = 1.5s around the coherence onset
    trl = trialOnsets+trialWindow;% results in [trial-beginning trial-end] !!NOTE we'd want to account for our projector delay here too if we weren't using the photodiode
    trl(:,3) = -1500;% add infor that epoched trials begin 1500ms before coherence onset
    alltrl = cat(1,alltrl,trl); % we want this later too
    tmplabs = cat(1,tmplabs,conditionlabels);
    conditionlabels = cellstr(conditionlabels); % need to convert for the preprocessing
    
    save(sprintf(trialEpochFile,num2str(runi)),'trl','conditionlabels')

    clear conditionlabels megRTs responseTimes % clear these at the end or it'll have trouble on the second run
end

megBehaviouralData(:,end+1) = 1:length(megBehaviouralData(:,1));
% check for consistency
% get the condition labels
allconditionlabels = megBehaviouralData(:,7);
if sum(strcmp(allconditionlabels,tmplabs)) ~= numel(tmplabs); error(['incongruent condition labels:' thisSubject.meg_labs{runi}]);end

save(megrtFile,'megBehaviouralData','alltrl','allconditionlabels')

else
    disp('megrt file found, and not overwriting')
end


return
end

function channels = getTriggersFromChannelLabels(ftData,channels)
% expects structure
% channels.(arbitrary name).labels = {'channel label' 'etc'};

% start by looping through the arbitrary names of our channel groupings
channelGroups = fieldnames(channels);
for groupidx = 1:numel(channelGroups)
    
    thisGroup = channelGroups{groupidx}; % make this more legible
    
    % loop through all the channels in the channel grouping
    % and grab the row number that corresponds to the label
    channels.(thisGroup).rowNum = [];
    for triggeridx = 1:numel(channels.(thisGroup).labels)
        
        thisLabel = channels.(thisGroup).labels{triggeridx}; % make this legible
        
        % convert the channel labels into the row number ft has assigned to this channel
        thisTrigger = find(strcmp(ftData.label,thisLabel));
        % pop that into new trigger structure for safekeeping
        channels.(thisGroup).rowNum = [channels.(thisGroup).rowNum,thisTrigger];
        
    end; clear triggeridx
    
    % get the content of those channels
    theseTriggerChannels = channels.(thisGroup).rowNum; % make this more legible
    channelContent = ftData.trial{1}(theseTriggerChannels,:); % note the {1}---would need to loop through that if this was not continuous data
    
    % so channelContent now has a row per trigger, each column a timepoint
    % channel values will be 0 or 5 (TTL on/off logic
    % https://en.wikipedia.org/wiki/Transistor%E2%80%93transistor_logic)
    % you can see if you plot it:
%         plot(channelContent(1,:));

    % if we have multiple triggers, we can pool them to get binary values
    if size(channelContent,1) > 1
        % in this case, we would treat row-wise stacks of 5 as binary
        % and we would want to convert those into more intelligible decimal values
%         allEventsAsDecimalValues = bi2de([(thisChannel>4)'],'right-msb'); % convert values above 5 as binary into decimal
        % we would then want to look for those in our new decimal-based channel
        % conversion:
%         decimalTriggerValuesOfInterest = [some,numbers,here];
%         triggerOnsetTimes =
%         (find(ismember(round(diff(allEventsAsDecimalValues)),decimalTriggerValuesOfInterest))+1)'; % logic of this explained in the else statement
        % incidentally, this is what STI101 does (although it also throws in button
        % responses as huge values I don't understand)
        % since I'm not pooling here, I'll error and we'd need to adjust this
        % function to accept something like:
        % channels.(group name).decimalCodes = [some,numbers,here];
        % use that to convert as above, then finishing off as we do
        % in the else statement below
        error('we havent coded for pooling trigger channels yet')
    else
        % if we're just using one channel per trigger we can just grab all
        % the onset times by looking for the index
        % of times when there is a difference between values in channelContent
        % (time points) of 5. we add one because this operation would otherwise grab the
        % timepoint BEFORE it changes to five, and we want it WHEN it is five
        triggerOnsetTimes = (find(ismember(round(diff(channelContent)),5))+1)';
        % we can see, if we plot it:
%         plot(channelContent); hold on
%         plot(triggerOnsetTimes,thisChannel(triggerOnsetTimes),'-o'); hold off
        
        % similarly, we can grab the trigger offsets by looking for the
        % indices for when the time points differ by -5
        triggerOffsetTimes = (find(ismember(round(diff(channelContent)),-5))+1)';
        
        
        % let's get these times out now
        channels.(thisGroup).onsetTimes = triggerOnsetTimes;
        channels.(thisGroup).offsetTimes = triggerOffsetTimes;
        
    end; clear triggeridx
    
    
end

return
end
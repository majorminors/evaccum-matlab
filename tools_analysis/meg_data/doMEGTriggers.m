function doMEGTriggers(filenames,behavfile,savebehav,settings,trfiles)


MEG_RT=[];
alltrl = [];
allconditionlabels = {};
tmplabs ={};

load(behavfile) % matlab behavioural file
runid = settings.runid; % pull information related to the blocks in each run

% runid looks like:
%     runid =
%
%      1     2   - is the run, and which blocks
%      3     4   - is run 2 with those blocks
%      5     6   - etc
%      7     8
%      9    10

conditions = {'LcLr', 'LcHr','HcLr','HcHr'}; % 2x2 coherence and rule

for runi = 1:numel(filenames)
    %%
    clear trl conditionlabels
    D = spm_eeg_load(filenames{runi});
    
    
    % double check that trigger responses correspond to behavioral results
    trig = D(D.indchannel('STI101'),:)/1e6;  % 1e6 new scaling introduced by SPM12;
    ttime = D.time;
    
    % pull the information for the blocks in the run
    vecrunid = runid(runi,:); % select vector of blocks in the run
    vecrunid(isnan(vecrunid)) = []; % get rid of any NaNs if they exist
    block = d.stim_mat_all(:,:,vecrunid); % pull the stimulus matrix info for those blocks
    block = permute(block,[1,3,2]); % swaps 3rd and 2nd dimensions
    block = reshape(block,[],size(d.stim_mat_all,2),1); % reshapes so rows are autocalculated, columns defined by stimulus matrix columns, and forces 3rd dimension into the 2d structure
    
    % find onset trials, coherence, trial type, and rt
    % these triggers need to change
    trial_type_triggers = d.stim_mat_all(:,9,1)'; % trial type triggers/coherence onset triggers
    cue_trigger = p.MEGtriggers.cue;
    iti_trigger = p.MEGtriggers.trial;
    keys = [4096 8192 16384 -32768]; % index medium ring little - on right button box
    
    cue_onset = (find(ismember(round(diff(trig)),cue_trigger))+1)'; % cue trigger onset (less useful because the onset might be cut off at the start and thus not always found preceding a sequence of 8 trials)
    cue_offset = (find(ismember(round(diff(trig)),-cue_trigger))+1)'; % cue trigger offset (more useful because it predicts a sequence of 8 trials, and is only likely to be missing if the corresponding sequence is also missing)
    iti_onset = (find(ismember(round(diff(trig)),iti_trigger))+1)'; % trial trigger onset
    trial_type_onset = (find(ismember(round(diff(trig)),trial_type_triggers))+1)'; % trial type trigger onset
    
    
    trialIDtrg  = [0; round(diff(trig)')]; trialID = trialIDtrg(trial_type_onset); % gets the amplitude of the trigger (which should be the thing you specified as the trigger) and then strips anything that is not a trial type trigger
    
    % see if we can't find naughty triggers
    multiTrigs=[];
    for ic = 1:numel(trial_type_onset) % go through each onset trigger
        
        % search before and after the onset trigger for triggers
        searchtwin=[0 round(diff(trig(trial_type_onset(ic)-5:trial_type_onset(ic)+5)))];
        ttrigger = find(searchtwin>1); % get the values in the search window
        if numel(ttrigger) > 1 % if there's more than one value
            multiTrigs = [multiTrigs;[ttrigger,ic]]; % catch them, and the 'trial' they represent
        end
        
    end
            
    if settings.deleteMultiTrigs(runi)
        trial_type_onset(multiTrigs(:,3)') = [];
        trialID(multiTrigs(:,3)') = [];
    end
    
    if settings.checkTrigs(runi) || (size(block,1) ~= numel(trial_type_onset))
        
        % plot the triggers
        figure;
        % plot the value of the triggerline at every timepoint
        plot(trig);
        hold on;
        % plot the triggers we discovered by searching for any valid trial type trigger values
        plot(trial_type_onset,trig(trial_type_onset),'-o');
        % plot what SHOULD be happening
        whatShouldHappen = reshape(squeeze(d.stim_mat_all(:,9,:)),[],1); % could also do block(:,9)
        plot(trial_type_onset,whatShouldHappen(1:numel(trial_type_onset)),'-g');
        % plot the trial IDs it thinks it found
        % notice that at the erroneous triggers, this thinks it has found a 2
        % and a 26, although that is not reflected in the actual trigger value
        % which is at 250 both times. Can't think of a reason for this unless
        % it's doing that thing where it took too long to get up to 250 and so
        % it thinks there's a trigger in the ascent of the ITI trigger value?
        plot(trial_type_onset,trialID(1:numel(trial_type_onset)),'-m'); % plot what trialID thinks it found
        plot(trial_type_onset([multiTrigs(:,3)']),trig(trial_type_onset([multiTrigs(:,3)'])),'b*'); % plot our naughty triggers
        % multiple triggers looks right - corresponds to erroneous triggers
        % notice though, that the value of these multi triggers (Var1) are both
        % identical and do not in either case appear to relate to what trialID (Var2) thinks they are:
        viewNaughties = table(multiTrigs(:,1:2),trialID(multiTrigs(:,3))); disp(viewNaughties);

        
        howManySequences = floor(numel(trial_type_onset)/8); % how many full sequences: num trials divided by 8 (8 trials between cues)
        howManyBlocks = floor(howManySequences/8); % how many full blocks: num full sequences divided by 8 (8 sequences of 8 trials [64 trials] in a block)
    
    end
    
    if settings.reduceTriggers(runi) == 1
        
        trial_type_onset = trial_type_onset(1:size(block,1));
        
    elseif setting.reduceTriggers(runi) == 2
        
        trial_type_onset = trial_type_onset(size(block,1):end);
        
    end
    
    
    % check now that we have the right number of triggers
    if size(block,1) ~= numel(trial_type_onset)
        error('Inconsistent number of onset triggers vs trials')
    end
    
    
    clear tmp_RT tKeyp responset
    
    % now we go through each trial onset trigger
    for ions = 1:numel(trial_type_onset)
        if ions < numel(trial_type_onset) % save trigger stream between this and the next trial onset trigger
            tmp    =trig(trial_type_onset(ions):(trial_type_onset(ions+1)-1));
        else
            tmp    =trig(trial_type_onset(ions):end);
        end
        
        % look for responses in that timeframe
        tmpres = (find(ismember(diff(tmp),keys),1,'first')+1)';%first to avoid multiple button presses % time of button press within trial period (iti to iti)
        if ~isempty(tmpres) % if response
            responset(ions)  = trial_type_onset(ions)+tmpres-1; % get time of button press within block
            if length(trial_type_onset) < ions % if this is after the number of trials (onset triggers) that it should be
                responset(ions) = 0; % zero it all out
                tmp_RT(ions) = 0;
            else % otherwise get the time between the coherence onset and the response onset
                tmp_RT(ions)       = ttime(responset(ions))-(ttime(trial_type_onset(ions))+.034);% 34ms projector's delay % this is reaction time (i.e. time of button press from coherence onset)
            end
        else % if there's no response, zero everything out
            responset(ions) = 0;
            tmp_RT(ions)      = 0;
        end
        
    end
    
    
    rts = d.rt(vecrunid,:)';rts = rts(:); % transpose relevant rows of d.rt and then place each block under one another in the same row
    
    idx = rts~=0 & rts>0; % where there is a rt
    
    tmp = [rts(idx) transpose(tmp_RT(idx))]; % this compares rts here with rts there
    tmp = corr(tmp); % make sure they correlate
    if tmp(2,1) <0.95 || isnan(tmp(2,1))   % takes the correlation from the output of corr
        error('large discrepancy between stim RT and MEG RT!')
    end
    
    % re-compose matrix of behavioural data
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
    MEG_RT = cat(1,MEG_RT,[block(:,[1 3 5 8 9 10]), conditionlabels, accuracyrow, rts, tmp_RT',run_tagger]);
    
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
    
    %trial definition for epoching
    
    % set the trial window, including buffer for edge effects during filtering
    % - in the textbook - lowest frequency will determine
    twind = [-1500 2500];%in samples at 1000Hz sampling % what time window 1500 = 1.5s around the coherence onset
    trl      = (trial_type_onset+34)+twind;%[trial-beginning trial-end] accounts for 34ms projector's delay.
    trl(:,3) = -1500;% epoched trials begin 1500ms before coherence onset
    
    %add to session trial def
    alltrl = cat(1,alltrl,trl);
    
    tmplabs = cat(1,tmplabs,conditionlabels);
    conditionlabels = cellstr(conditionlabels); % need to convert for the preprocessing
    save(trfiles{runi},'trl','conditionlabels')
end

%%

% add trial numbers to MEG_RT
MEG_RT(:,end+1) = 1:length(MEG_RT(:,1));
% now use those to create something that we can tag the meeg data with

% get the condition labels
allconditionlabels = MEG_RT(:,7);

% check for consistency
if sum(strcmp(allconditionlabels,tmplabs)) ~= numel(tmplabs); error(['incongruent condition labels:' filenames{runi}]);end

save(savebehav,'MEG_RT','settings','alltrl','allconditionlabels')

return
end

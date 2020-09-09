function fun_MEGtrgRT(filenames,behavfile,savebehav,settings,trfiles)


MEG_RT=[];
alltrl = [];
allconditionlabels = {};
tmplabs ={};

load(behavfile) % matlab behavioural file

% runid = %%%% we want to pull this from the organisation %%%%

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
    triggers = d.stim_mat_all(:,10,1)'; % trial type triggers
    coh_trigger = 1; % value of coherence onset trigger
    keys = [4096 8192 16384 -32768]; % index medium ring little - on right button box
    trigger_onset = (find(round(diff(trig)) == 1)+1)'; % onset of all triggers - SUPERFLOUS
    resp_onset_neg = (find(ismember(round(diff(trig)),-keys))+1)';
    coh_onset = (find(ismember(round(diff(trig)),coh_trigger))+1)'; % coherence trigger onset
    trial_type_onset = (find(ismember(round(diff(trig)),triggers))+1)'; % trial type trigger onset
    trial_type_onset_neg = (find(ismember(round(diff(trig)),-triggers))+1)'; % trial type trigger onset
    
    cohtrigsID  = [0; round(diff(trig)')]; cohtrigsID = cohtrigsID(coh_onset); % gets the amplitude of the trigger (which should be the thing you specified as the trigger) and then strips anything that is not a coherence onset trigger
    trialIDtrg  = [0; round(diff(trig)')]; trialID = trialIDtrg(trial_type_onset); % same for trial type
    
    if numel(trialID) ~= numel(block(:,10)) || sum(trialID ~= block(:,10))
        
        if numel(coh_onset) == numel(block(:,10)) % if coh  onset is reliable
            tmp_tronset = [];
            tmp_trialID =[];
            for ic = 1:numel(coh_onset)
                
                searchtwin=[0 round(diff(trig(coh_onset(ic)-400:coh_onset(ic))))];
                ttrigger = find(searchtwin>1 & searchtwin <= max(triggers));
                ttrigger = ttrigger(1);
                
                tmp_trialID(ic) = block(ic,10);
                
                
                tmp_tronset(ic) = coh_onset(ic)-400 + ttrigger-1;
                
            end
            trial_type_onset   = tmp_tronset';
            trialIDtrg  =  tmp_trialID';
        end

    end
    
    % plot the triggers
    figure; plot(trig); hold on; plot(trial_type_onset,trig(trial_type_onset),'-o');

    %check
    if numel(trigger_onset) ~= size(block,1)
        error('Inconsistent number of onset triggers vs trials')
    end

    % check
    if numel(trialID) ~= numel(block(:,10)) || sum(trialID ~= block(:,10))
        if sum(trialIDtrg ~= block(:,10)) % compare with the d.stim_mat_all column with trigger amplitudes
            error('Discrepancy between triggers coded conditions and trials')
        end
    end
    
    % find responses and empty trials
    clear tmp_RT tKeyp responset
    
    for ions = 1:numel(trial_type_onset)
        if ions < numel(trial_type_onset)
            tmp    =trig(trial_type_onset(ions):(trial_type_onset(ions+1)-1));
        else
            tmp    =trig(trial_type_onset(ions):end);
        end
        
        tmpres = (find(ismember(diff(tmp),keys),1,'first')+1)';%first to avoid multiple button presses % time of button press within trial period (iti to iti)
        if ~isempty(tmpres)
            responset(ions)  = trial_type_onset(ions)+tmpres-1;%#ok % time of button press within block
            if length(coh_onset) < ions
                responset(ions) = 0;
                tmp_RT(ions) = 0;
            else
                tmp_RT(ions)       = ttime(responset(ions))-(ttime(coh_onset(ions))+.034);%#ok 34ms projector's delay % this is reaction time (i.e. time of button press from coherence onset)
            end
        else
            responset(ions) = 0;%#ok
            tmp_RT(ions)      = 0;%#ok
        end
        
    end
    
    %check
    % rts = block(:,13); % will need to get rts from d.rts (this is a column in this script)
    
    rts = d.rt(vecrunid,:)';rts = rts(:); % transpose relevant rows of d.rt and then place each block under one another in the same row
    
    idx = rts~=0 & rts>0; % where there is a rt
    
    tmp = [rts(idx) transpose(tmp_RT(idx))]; % this compares rts here with rts there
    tmp = corr(tmp);
    
    
    if tmp(2,1) <0.99 || isnan(tmp(2,1))   % takes the correlation from the output of corr
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
    MEG_RT = cat(1,MEG_RT,[block(:,[1 3 5 8 10]), conditionlabels, accuracyrow, rts, tmp_RT',run_tagger]);

    % re-composed matrix looks like:
    %  1)  cue direction (1-4)
    %  2)  dot motion direction condition (1-8)
    %  3)  coherence difficulty (1 = easy, 2 = hard)
    %  4)  matching difficulty (1 = easy, 2 = difficult)
    %  5)  gives a number based on 9 for meg triggers
    %  6)  conditions (HcHr, HcLr, LcHr, LcLr)
    %  7)  accuracy (0, 1, or -1 for missed trials)
    %  8)  rts from behavioural data
    %  9)  rts from MEG
    % 10)  something to tag which run this sequence of data belongs to
    % 11)  we later add a row of unique numbers to tag each trial
    
    %trial definition for epoching
    
    % set the trial window, including buffer for edge effects during filtering
    % - in the textbook - lowest frequency will determine
    twind = [-1500 2500];%in samples at 1000Hz sampling % what time window 1500 = 1.5s around the coherence onset
    trl      = (coh_onset+34)+twind;%[trial-beginning trial-end] accounts for 34ms projector's delay.
    trl(:,3) = -1500;% epoched trials begin 1500ms before coherence onset
    
    %add to session trial def
    alltrl = cat(1,alltrl,trl);
    
    tmplabs = cat(1,tmplabs,conditionlabels);
    conditionlabels = cellstr(conditionlabels); % need to convert for the preprocessing
    save(trfiles{runi},'trl','conditionlabels')
end

% don't think we need this
% idbad = MEG_RT(:,7)<=0 | MEG_RT(:,7)>2.1; % defines bad trials based on reaction time will need to check long RT and column number
% MEG_RT(idbad,8) = 0;
% MEG_RT(idbad,[6 7]) = -1;

%%

% add trial numbers to MEG_RT
MEG_RT(:,end+1) = 1:length(MEG_RT(:,1));
% now use those to create something that we can tag the meeg data with

%allconditionlabels = {conditions{MEG_RT(:,8)}}';%#ok
allconditionlabels = MEG_RT(:,6);

%check for consistency
if sum(strcmp(allconditionlabels,tmplabs)) ~= numel(tmplabs); error(['incongruent condition labels:' filenames{runi}]);end

save(savebehav,'MEG_RT','settings','alltrl','allconditionlabels')

function fun_MEGtrgRT(filenames,behavfile,savebehav,settings,trfiles)%#ok


MEG_RT=[];alltrl = []; allconditionlabels = {};tmplabs ={};%#ok
load(behavfile) % matlab behavioural file

runid = repmat((1:numel(filenames))',[1,2]).*2-repmat([1 0],[numel(filenames),1]); % tells the block ids within one run

% runid looks like:
%     runid =
% 
%      1     2   - is the run, and which blocks
%      3     4   - is run 2 with those blocks
%      5     6   - etc
%      7     8
%      9    10

conditions = {'pLaL', 'pHaL','pLaH','pHaH'}; % 2x2 coherence and rule


for runi = 1:numel(filenames)
%%
clear trl conditionlabels
D = spm_eeg_load(filenames{runi});


% double check that trigger responses correspond to behavioral results
trig = D(D.indchannel('STI101'),:)/1e6;  % 1e6 new scaling introduced by SPM12;
ttime = D.time;


    

block = ResponseArray(ismember(ResponseArray(:,1),runid(runi,:)),:); % (d.stim_mat_all - this selects the elements from the behavioural file from the block/s associated with the run with runid + the reaction time data - so runi on the 3rd element

%find onset trials, coherence, trial type, and rt
triggers = [4 8 16 32]; % trial type triggers
coh_trigger = []; % value of coherence onset trigger
keys     = [4096 8192 16384 -32768]; % index medium ring little - on right button box
trigger_onset  = (find(round(diff(trig)) == 1)+1)'; % onset of all triggers
coh_onset = (find(ismember(round(diff(trig)),coh_trigger))+1)'; 
trial_type_onset = (find(ismember(round(diff(trig)),triggers))+1)'; 
trigsID  = [0; round(diff(trig)')]; trigsID = trigsID(coh_onset); % gets the amplitude of the trigger (which should be the thing you specified as the trigger) and then strips anything that is not a coherence onset trigger
trialID  = [0; round(diff(trig)')]; trialID = trialID(trial_type_onset); % same for trial type
%check
if numel(trigger_onset)~= size(block,1)
    error('Inconsistent number of onset triggers vs trials')
end
if sum(trialID ~= block(:,12)) % compare with the d.stim_mat_all column with trigger amplitudes
     error('Discrepancy between triggers coded conditions and trials')
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
        tmp_RT(ions)       = ttime(responset(ions))-(ttime(coh_onset(ions))+.034);%#ok 34ms projector's delay % this is reaction time (i.e. time of button press from coherence onset)
    else
        responset(ions) = 0;%#ok
        tmp_RT(ions)      = 0;%#ok
    end
        
end

%check
rts = block(:,13); % will need to get rts from d.rts (this is a column in this script)

idx = rts~=0 & rts~=Inf & rts>0; % inf is no reply - you will want -1

tmp = [rts(idx) transpose(tmp_RT(idx))]; % this compares rts here with rts there
tmp = corr(tmp);

if tmp(2,1) <0.99 || isnan(tmp(2,1))   % takes the correlation from the output of corr  
     error('large discrepancy between stim RT and MEG RT!')     
end    


%re-compose matrix of behavioural data
MEG_RT = cat(1,MEG_RT,[block(:,[1:4 12 14 15 ]), tmp_RT']);% 1block# 2coh 3act 4cond 5trigger 6resKey 7acc 8rt % compiles matlab trial information with the new MEG generated rts

%trial definition for epoching
if sum(block(:,12) ~= trigsID);error('discrepancy between MEGRT and MEGbehav condition labels!');end

% set the trial window, including buffer for edge effects during filtering
% - in the textbook - lowest frequency will determine
twind = [-1500 2500];%in samples at 1000Hz sampling % what time window 1500 = 1.5s around the coherence onset
trl      = (coh_onset+34)+twind;%[trial-beginning trial-end] accounts for 34ms projector's delay.
trl(:,3) = -1500;% epoched trials begin 1500ms before coherence onset

%add to session trial def
alltrl = cat(1,alltrl,trl);

conditionlabels = {conditions{block(:,4)}}';%#ok % from the d.stim_mat.all - give the condition names for each trial in the block (h/l coherence and rule)
tmplabs = cat(1,tmplabs,conditionlabels);

save(trfiles{runi},'trl','conditionlabels')
end

idbad = ResponseArray(:,13)<=0 | MEG_RT(:,end)>2.1 | isinf(MEG_RT(:,6)); % defines bad trials based on reaction time
MEG_RT(idbad,8) = 0;
MEG_RT(idbad,[6 7]) = -1;

%%

allconditionlabels = {conditions{MEG_RT(:,4)}}';%#ok

%check for consistency
if sum(strcmp(allconditionlabels,tmplabs)) ~= numel(tmplabs); error(['incongruent condition labels:' filenames{runi}]);end


save(savebehav,'MEG_RT','settings','alltrl','allconditionlabels')
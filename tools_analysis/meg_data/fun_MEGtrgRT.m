function fun_MEGtrgRT(filenames,behavfile,savebehav,settings,trfiles)%#ok


MEG_RT=[];alltrl = []; allconditionlabels = {};tmplabs ={};%#ok
load(behavfile) % matlab behavioural file

%runid = repmat((1:numel(filenames))',[1,2]).*2-repmat([1 0],[numel(filenames),1]); % tells the block ids within one run

% runid looks like:
%     runid =
% 
%      1     2   - is the run, and which blocks
%      3     4   - is run 2 with those blocks
%      5     6   - etc
%      7     8
%      9    10

% S01 - 6 runs, with 6 blocks per run
runid = [1:6;7:12;13:18;19:24;25:30;31:36];

% S02 - 3 runs, run 1 no data, run 2 is two blocks, run 3 is 4 blocks
% runid = [[];1:2;3:6];
% runid = NaN(3,4); runid(2,1:2) = 1:2; runid(3,:) = 3:6;

% S03 - 5 runs, 6 blocks per run ex last run is 1 block, and first 16 trials of run 4 belong to the end of the last block in run 3
% runid = [1:6;7:12;13:18;19:24;25]; % so block 19 has 16 trials from end of block 18
% runid = [1:6;7:12;13:18;19:24;25,NaN,NaN,NaN,NaN,NaN];

conditions = {'LcLr', 'LcHr','HcLr','HcHr'}; % 2x2 coherence and rule

for runi = 1:numel(filenames)
%%
clear trl conditionlabels
D = spm_eeg_load(filenames{runi});


% double check that trigger responses correspond to behavioral results
trig = D(D.indchannel('STI101'),:)/1e6;  % 1e6 new scaling introduced by SPM12;
ttime = D.time;


    

% ale's code - block = ResponseArray(ismember(ResponseArray(:,1),runid(runi,:)),:); % (d.stim_mat_all - this selects the elements from the behavioural file from the block/s associated with the run with runid + the reaction time data - so runi on the 3rd element
% block = d.stim_mat_all(:,:,runid(runi,:));
% so can do the above and work with 3 dimensions, but requires us to change block refs to suit, eg :
% if numel(trigger_onset)~= size(block,1)*size(block,3)
% if sum(trialID ~= block(:,10,:))
% or can reshape d.stim_mat_all as below, but concerned I'll lose block related information? maybe can add a row...

block = permute(d.stim_mat_all,[1,3,2]); % swaps 3rd and 2nd dimensions
block = reshape(block,[],size(d.stim_mat_all,2),1); % reshapes so rows are autocalculated, columns defined by testmat columns, and forces 3rd dimension into the 2d structure


% find onset trials, coherence, trial type, and rt
% for S01 (no separate trial type and coh onset) need to swap refs to other triggers for trigger_onset, and figure out what responset should be, also triggers are in column 9, not 10
triggers = d.stim_mat_all(:,10)'; % trial type triggers
coh_trigger = 1; % value of coherence onset trigger
keys     = [4096 8192 16384 -32768]; % index medium ring little - on right button box
trigger_onset  = (find(round(diff(trig)) == 1)+1)'; % onset of all triggers
coh_onset = (find(ismember(round(diff(trig)),coh_trigger))+1)'; % coherence trigger onset
trial_type_onset = (find(ismember(round(diff(trig)),triggers))+1)'; % trial type trigger onset
cohtrigsID  = [0; round(diff(trig)')]; cohtrigsID = cohtrigsID(coh_onset); % gets the amplitude of the trigger (which should be the thing you specified as the trigger) and then strips anything that is not a coherence onset trigger
trialID  = [0; round(diff(trig)')]; trialID = trialID(trial_type_onset); % same for trial type

%check
if numel(trigger_onset)~= size(block,1)
    error('Inconsistent number of onset triggers vs trials')
end
if sum(trialID ~= block(:,10)) % compare with the d.stim_mat_all column with trigger amplitudes
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
% rts = block(:,13); % will need to get rts from d.rts (this is a column in this script)
rts = reshape(d.rt,[],1); % transpose d.rt and place each block under one another in the same row

idx = rts~=0 & rts>0; % where there is a rt

tmp = [rts(idx) transpose(tmp_RT(idx))]; % this compares rts here with rts there
tmp = corr(tmp);

if tmp(2,1) <0.99 || isnan(tmp(2,1))   % takes the correlation from the output of corr  
     error('large discrepancy between stim RT and MEG RT!')     
end    


%re-compose matrix of behavioural data
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
condnumbs = []; for i = 1:length(block) % create a row that gives you a number for each condition in your 2x2
    if block(i,5)==1 && block(i,8)==1
        condnumbs(i,:) = 1;
    elseif block(i,5)==1 && block(i,8)==2
        condnumbs(i,:) = 2;
    elseif block(i,5)==2 && block(i,8)==1
        condnumbs(i,:) = 3;
    elseif block(i,5)==2 && block(i,8)==2
        condnumbs(i,:) = 4;
    end
end
accuracyrow = reshape(d.correct,[],1); % may as well put this in the new matrix
%MEG_RT = cat(1,MEG_RT,[block(:,[1:4 12 14 15 ]), tmp_RT']);% 1block# 2coh 3act 4cond 5trigger 6resKey 7acc 8rt % compiles matlab trial information with the new MEG generated rts
MEG_RT = cat(1,MEG_RT,[block(:,[1 3 5 8 10]), accuracyrow, rts, condnumbs, tmp_RT']);

%trial definition for epoching
if sum(block(:,1) ~= cohtrigsID);error('discrepancy between MEGRT and MEGbehav condition labels!');end

% set the trial window, including buffer for edge effects during filtering
% - in the textbook - lowest frequency will determine
twind = [-1500 2500];%in samples at 1000Hz sampling % what time window 1500 = 1.5s around the coherence onset
trl      = (coh_onset+34)+twind;%[trial-beginning trial-end] accounts for 34ms projector's delay.
trl(:,3) = -1500;% epoched trials begin 1500ms before coherence onset

%add to session trial def
alltrl = cat(1,alltrl,trl);

conditionlabels = {conditions{MEG_RT(:,8)}}';%#ok % will need to ensure correct - not sure how MEG_RT works... from the d.stim_mat.all - give the condition names for each trial in the block (h/l coherence and rule)
tmplabs = cat(1,tmplabs,conditionlabels);

save(trfiles{runi},'trl','conditionlabels')
end

% don't think we need this
% idbad = MEG_RT(:,7)<=0 | MEG_RT(:,7)>2.1; % defines bad trials based on reaction time will need to check long RT and column number
% MEG_RT(idbad,8) = 0;
% MEG_RT(idbad,[6 7]) = -1;

%%

allconditionlabels = {conditions{MEG_RT(:,8)}}';%#ok

%check for consistency
if sum(strcmp(allconditionlabels,tmplabs)) ~= numel(tmplabs); error(['incongruent condition labels:' filenames{runi}]);end


save(savebehav,'MEG_RT','settings','alltrl','allconditionlabels')
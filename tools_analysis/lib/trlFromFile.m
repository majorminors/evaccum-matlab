function trl = trlFromFile(cfg)
% a custom function to make a fieldtrip definetrial trl structure from an
% spm/osl trl file with condition codes
% assumes you have saved the file with
%    trl = [trialStart trialEnd trialOffset; etc] for each trial
%    conditionlabels = {'condition1' 'condition2'... for each trial
% so you can specify:
%    add ms pre trial window (cfg.pre)
%    add ms post trial window (cfg.pre)
%    lock the 'onset' time to response (cfg.lockTo = 'response') (uses the response timing
%        from a behavioural file and then adds any post timing)
%    lock the 'onset' time to fixation (cfg.lockTo = 'fixation') (similarly
%        uses the fixation timing from a behavioural file), CURRENTLY
%        OVERWRITES POST TIMING
%    doesn't do anything for locking to onset because this assumes you
%        already did that in the trl file

%% check fields
if ~isfield(cfg,'pre')
    cfg.pre = 0;
end
if ~isfield(cfg,'post')
    cfg.post = 0;
end
if ~isfield(cfg,'lockTo')
    cfg.lockTo = 'onset';
end

%% load our behavioural data
load(cfg.trlfile);
tmp = load(cfg.behavfile);
behavTable = tmp.megBehaviouralData;
clear tmp

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
% 12)  we later add the fixation dots timings for each trial
% 13)  we later add a row of unique numbers to tag each trial

%% filter by run and good trials
if isfield(cfg,'run')
    runIdx = find(behavTable(:,11)==cfg.run);
    goodTrials = str2double(behavTable(runIdx,10))>0;
else
    goodTrials = str2double(behavTable(:,10))>0;
end

%% alter trial 'start'
% my trl file is already locked to onset, so I don't need to define this
% here
if strcmp(cfg.lockTo,'response')
    trl(:,1) = trl(:,1)+str2double(behavTable(runIdx,10));
    cfg.post = -(trl(:,2) - (trl(:,1)+cfg.post));
elseif strcmp(cfg.lockTo,'fixation')
    trl(:,1) = trl(:,1)-(str2double(behavTable(runIdx,12))*1000);
    cfg.post = -(str2double(behavTable(runIdx,12))*1000);
end

%% now we need to convert condition labels to codes
for idx = 1:numel(conditionlabels)
    if strcmp(conditionlabels{idx},'EcEr')
        conditioncodes(idx) = 1;
    elseif strcmp(conditionlabels{idx},'EcHr')
        conditioncodes(idx) = 2;
    elseif strcmp(conditionlabels{idx},'HcEr')
        conditioncodes(idx) = 3;
    elseif strcmp(conditionlabels{idx},'HcHr')
        conditioncodes(idx) = 4;
    else
        error('have you edited this properly to code your condition labels?')
    end
end

%% add any pre and post timings
trl(:,1) = trl(:,1)+cfg.pre;
trl(:,2) = trl(:,2)+cfg.post;
trl(:,3) = cfg.pre;

%% recompose trl variable
trl = [trl conditioncodes' goodTrials];
trl = fix(trl); % since matlab does it's absolute best to confuse stuff with e-notation, thus fucking with its own functions which don't, of course, understand e-notation

end
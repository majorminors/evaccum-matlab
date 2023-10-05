function trl = trlFromFile(cfg)
% a custom function to make a fieldtrip definetrial trl structure from an
% spm/osl trl file with condition codes
% assumes you have saved the file with
%    trl = [trialStart trialEnd trialOffset; etc] for each trial
% you can optionally specify:
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
% filter by run
if isfield(cfg,'run')
    behavTable = tmp.megBehaviouralData(find(tmp.megBehaviouralData(:,11)==cfg.run),:);
else
    behavTable = tmp.megBehaviouralData;
end
clear tmp

%% now we need to add some stuff in from the behavioural results
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
% 12)  the fixation dots timings for each trial
% 13)  a row of unique numbers to tag each trial

% get good trials (not missed, and not zero rt)
goodTrials = str2double(behavTable(:,8))>-1 & str2double(behavTable(:,10))>0;

% get condition codes
for idx = 1:size(trl,1)
    if strcmp(behavTable(idx,7),'EcEr')
        conditioncodes(idx) = 1;
    elseif strcmp(behavTable(idx,7),'EcHr')
        conditioncodes(idx) = 2;
    elseif strcmp(behavTable(idx,7),'HcEr')
        conditioncodes(idx) = 3;
    elseif strcmp(behavTable(idx,7),'HcHr')
        conditioncodes(idx) = 4;
    else
        error('have you edited this properly to code your condition labels?')
    end
end

%% alter trial 'start'
% my trl file is already locked to onset, so I don't need to define this
% here
if strcmp(cfg.lockTo,'response')
    trl(:,1) = trl(:,1)+str2double(behavTable(:,10));
    cfg.post = -(trl(:,2) - (trl(:,1)+cfg.post));
elseif strcmp(cfg.lockTo,'fixation')
    trl(:,1) = trl(:,1)-(str2double(behavTable(:,12))*1000);
    cfg.post = -(str2double(behavTable(:,12))*1000);
end


%% add any pre and post timings
trl(:,1) = trl(:,1)+cfg.pre;
trl(:,2) = trl(:,2)+cfg.post;
trl(:,3) = cfg.pre;

%% recompose trl variable
trl = [trl conditioncodes' goodTrials];
trl = fix(trl); % since matlab does it's absolute best to confuse stuff with e-notation, thus fucking with its own functions which don't, of course, understand e-notation

end
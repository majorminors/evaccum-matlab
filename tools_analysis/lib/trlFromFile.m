function trl = trlFromFile(cfg)
% a custom function to make a fieldtrip definetrial trl structure from an
% spm/osl trl file with condition codes
% assumes you have saved the file with
%    trl = [trialStart trialEnd trialOffset; etc] for each trial
%    conditionlabels = {'condition1' 'condition2'... for each trial

if ~isfield(cfg,'pre')
    cfg.pre = 0;
end
if ~isfield(cfg,'post')
    cfg.post = 0;
end

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

% let's grab the run idx

if isfield(cfg,'run')
    runIdx = find(behavTable(:,11)==cfg.run);
    goodTrials = str2double(behavTable(runIdx,10))>0;
else
    goodTrials = str2double(behavTable(:,10))>0;
end

if ischar(cfg.pre)
    if strcmp(cfg.pre,'fixation')
        cfg.pre = -(str2double(behavTable(runIdx,12))*1000); % start the 'trial' when the fixation dots started by subtracting the fixation times from the onset
    elseif strcmp(cfg.pre,'response')
        cfg.pre = (trl(:,1)-(trl(:,1)+str2double(behavTable(runIdx,10))))-500; % start the trial at trial onset+RT -500ms (so we can look for a plateau etc)
    end
end
if ischar(cfg.post)
    if strcmp(cfg.post,'response')
        cfg.post = -(trl(:,2)-(trl(:,1)+str2double(behavTable(runIdx,10))))+200; % finish the trial when the response happened by getting the negative of trial end - (trial start+rt)
%         cfg.pre = -(trl(:,1)-((trl(:,2)+cfg.post)-1000)); % get whatever 1000 ms prior to response is relative to the trl onset
    end
end

% now we need to convert condition labels to codes
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

trl(:,1) = trl(:,1)+cfg.pre;
trl(:,2) = trl(:,2)+cfg.post;
trl(:,3) = cfg.pre;

trl = [trl conditioncodes' goodTrials];
trl = fix(trl); % since matlab does it's absolute best to confuse stuff with e-notation, thus fucking with its own functions which don't, of course, understand e-notation

end
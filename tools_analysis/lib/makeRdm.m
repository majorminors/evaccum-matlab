function outRdm = makeRdm(model,theseTrials)
% this loads rdms generated with generate_rdms file

assert(isvector(theseTrials), 'trial ids must be a vector');
if  size(theseTrials, 1) > 1; theseTrials = theseTrials'; end % make it a row

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';
toolsdir = fullfile(rootdir,'tools_analysis');
modeldir = fullfile(toolsdir,'rdms');

%% grab model
switch model
    case 'stim'
        % Stimulus (motion) representation fit through time
        templateRdm = csvread(fullfile(modeldir,'rdm_stim.csv'));
    case 'decbdry'
        % Categorisation (decision boundary) representation fit through time
        templateRdm = csvread(fullfile(modeldir,'rdm_decbdry.csv'));
    case 'dec_detail_null'
        % Categorisation (decision boundary) motion moving on same or diff side of decision boundary with no prediction for perpendicular boundary
        templateRdm = csvread(fullfile(modeldir,'rdm_dec_detail_null.csv'));
    case 'dec_detail_pred'
        % Categorisation (decision boundary) motion moving on same or diff side of decision boundary with predicition for perpendicular boundary
        templateRdm = csvread(fullfile(modeldir,'rdm_dec_detail_pred.csv'));
    case 'dec_simple'
        % Categorisation (decision boundary) motion moving on same or diff side of decision boundary but cues are all different to each other
        templateRdm = csvread(fullfile(modeldir,'rdm_dec_simple.csv'));
    case 'resp'
        % Categorisation (?) (button press) representation fit through time
        templateRdm = csvread(fullfile(modeldir,'rdm_resp.csv'));
    otherwise
        error('no model found')
end

%% load trial ids
load([modeldir filesep 'trialIds.mat'])
rows = trialId';
cols = trialId;
clear trialId;

outRdm = zeros(length(theseTrials)); % preallocate
rowCount = 0; colCount = 0; % start counters
for thisRow = theseTrials % go row by row
    rowCount = rowCount+1; % iterate row count
    for thisCol = theseTrials % go column by column
        colCount = colCount+1; % iterate column counter
        if colCount > numel(theseTrials); colCount = 1; end % if we've exceeded the number of trials, we're on a new row so reset this

        outRdm(rowCount,colCount) = templateRdm(find(rows==thisRow),find(cols==thisCol)); % grab the corresponding value
        
    end
end

return
end

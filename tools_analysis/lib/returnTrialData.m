function [trialAverage, acrossSensors] = returnTrialData(dataObject,channelLabels,timePoints,trialIndices)

if ~exist('timePoints','var') % get all timepoints
    timePoints = size(dataObject,2);
end
if ~exist('trialIndices','var') % get all epochs/trials
    trialIndices = size(dataObject,3);
end
if ~isempty(channelLabels)
    if ~iscell(channelLabels)
        labelFunction = @contains;
    else
        labelFunction = @ismember;
    end
    sensors = find(labelFunction(chanlabels(dataObject),channelLabels));
else
    sensors = 1:size(dataObject,1);
end

% so spm data objects are sensor by time by trial

acrossSensors = mean(dataObject(sensors,timePoints,trialIndices),3); % will mean across trials
trialAverage = mean(acrossSensors);

end


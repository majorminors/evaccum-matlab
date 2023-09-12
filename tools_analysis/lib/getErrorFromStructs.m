
function errorMat = getErrorFromStructs(data)
% for getting the std err (errorMat.stderr) or std deviation (errorMat.stddev)
% from multiple fieldtrip data structures
% assumes they're all the same channels + labels etc

% grab some initial data, so we don't have to refer to the original
% structure for that
errorMat.time = data{1}.time;
errorMat.label = data{1}.label;


numStructures = numel(data); % get the number of structures in data
matrixSize = size(data{1}.avg); % get the size of the avg matrix in the first structure

% init the matrix to store average values from all structures
accumulatedMatrix = zeros(matrixSize);

% accumulate the values from all structures
for i = 1:numStructures
    avgMatrix = data{i}.avg; % Get the avg matrix for the current structure
    accumulatedMatrix = accumulatedMatrix + avgMatrix; % Accumulate the values
end

% calc the average across all structures
avgAcrossStructures = accumulatedMatrix / numStructures;

% calc the standard deviation across all structures
errorMat.stddev = zeros(matrixSize);
for i = 1:numStructures
    avgMatrix = data{i}.avg; % get the avg matrix for the current structure
    errorMat.stddev = errorMat.stddev + (avgMatrix - avgAcrossStructures).^2; % accumulate the squared differences
end
errorMat.stddev = sqrt(errorMat.stddev / numStructures); % take the square root to calculate the standard deviation

% calc the standard error across all structures
errorMat.stderr = zeros(matrixSize);
for i = 1:numStructures
    avgMatrix = data{i}.avg; % get the avg matrix for the current structure
    errorMat.stderr = errorMat.stderr + (avgMatrix - avgAcrossStructures).^2; % accumulate the squared differences
end
errorMat.stderr = sqrt(errorMat.stderr / (numStructures-1) / numStructures); % calc the standard error

return
end
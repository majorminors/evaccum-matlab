
function errorMat = getErrorFromStructs(data,type)
% for getting the std err or std deviation from multiple fieldtrip data structures

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

switch type
    case 'deviation'
        % calc the standard deviation across all structures
        errorMat = zeros(matrixSize);
        for i = 1:numStructures
            avgMatrix = data{i}.avg; % get the avg matrix for the current structure
            errorMat = errorMat + (avgMatrix - avgAcrossStructures).^2; % accumulate the squared differences
        end
        errorMat = sqrt(errorMat / numStructures); % take the square root to calculate the standard deviation
    case 'error'   
        % calc the standard error across all structures
        errorMat = zeros(matrixSize);
        for i = 1:numStructures
            avgMatrix = data{i}.avg; % get the avg matrix for the current structure
            errorMat = errorMat + (avgMatrix - avgAcrossStructures).^2; % accumulate the squared differences
        end
        errorMat = sqrt(errorMat / (numStructures-1) / numStructures); % calc the standard error
end

return
end
function outTable = makeTableWithNans(varargin)

titles = varargin{1};

varargin(1) = [];

if numel(titles) ~= numel(varargin); error('titles and variables not aligned'); end

numVars = numel(varargin);

% first find out how long the longest row will need to be:
for i = 1:numVars
    
    % check for outliers and NaN them
    if ~ischar(varargin{i}) && ~iscellstr(varargin{i})
%         varargin{i}(isoutlier(varargin{i})) = NaN;
    end
    
    % get the upper limit of elements from the first participant (i.e. how long
    % our longest row will have to be)
    if i == 1
        upperLimit = numel(varargin{i});
        continue
    end
    
    % now check whether we need to increase that upper limit with the
    % subsequent participants
    
    thisSize = numel(varargin{i});
    
    if thisSize > upperLimit
        upperLimit = thisSize;
    end
    
end

% now we know how long our longest row will need to be:
for i = 1:numVars
    
    % make sure everything is in row form
    if isrow(varargin{i})
        varargin{i} = varargin{i}';
    end
    
    % add nans to fill out the row if it's too short
    if numel(varargin{i}) < upperLimit
        nansToAdd = upperLimit-numel(varargin{i});
        varargin{i} = [varargin{i};nan(nansToAdd,1)];
    end

end

% format the table
outTable = table(varargin{:});
outTable.Properties.VariableNames = titles;

return
end
function plotErrorFromStructs(dataStructs,colours,channels,errorMat)

% loop through them

for i = 1:numel(errorMat)
    
    
    % now grab the error from the channels we care about
    selected_data = errorMat{i}.stderr(find(ismember(errorMat{i}.label,channels)),:);
    
    % get the mean of this, if you're doing multiple channels
    if numel(channels) > 1
        selected_data = mean(selected_data);
    end
    selected_time = dataStructs{i}.time;
    
    % grab the data to centre the error on
    centre_data = dataStructs{i}.avg(find(ismember(dataStructs{i}.label,channels)),:);
    if size(centre_data,1) > 1 % average it if there is more than one channel
        centre_data = mean(centre_data);
    end
    
    % define the lower and upper bounds of the area
    lowerBound = centre_data - selected_data/2;
    upperBound = centre_data + selected_data/2;
    
    % create the area plot using the fill function
    fill([selected_time, fliplr(selected_time)],...
        [lowerBound, fliplr(upperBound)], colours(i,:),...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
end

return
end

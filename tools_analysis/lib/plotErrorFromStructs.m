function plotErrorFromStructs(figHandle,dataStructs,colours,channels,errorMat)

axesData = get(figHandle.CurrentAxes, 'Children');

for i = 1:numel(axesData)
    
    % since the line axes go in reverse plot order, get a reverse index
    ir = numel(axesData) - i + 1;
    
    % get the data to centre on from the y-axis
    centre_data = get(axesData(i), 'YData');
    
    % now grab the error from the channels we care about
    selected_data = errorMat{ir}(find(ismember(dataStructs{ir}.label,channels)),:);
    
    % get the mean of this, if you're doing multiple channels
    if numel(channels) > 1
        selected_data = mean(selected_data);
    end
    selected_time = dataStructs{ir}.time;
    
    % just for record, we can re-do the mean from the data instead of
    % pulling it from the plot
    % centre_data = ecerOnsAve.avg(find(ismember(ecerOnsAve.label,CPP)),:);
    % centre_data = mean(centre_data);
    
    % define the lower and upper bounds of the area
    lowerBound = centre_data - selected_data/2;
    upperBound = centre_data + selected_data/2;
    
    % create the area plot using the fill function
    fill([selected_time, fliplr(selected_time)],...
        [lowerBound, fliplr(upperBound)], colours(ir,:),...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
end

return
end

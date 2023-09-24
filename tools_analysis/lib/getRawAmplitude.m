function rawAmplitude = getRawAmplitude(dataStructure, channels)

    selectedChannelData = dataStructure.avg(ismember(dataStructure.label, channels), :);
    if size(selectedChannelData, 1) > 1
        rawAmplitude = mean(selectedChannelData, 1);
    else
        rawAmplitude = selectedChannelData;
    end
    
end
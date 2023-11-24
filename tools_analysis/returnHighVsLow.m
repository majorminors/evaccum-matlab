function [high low] = returnHighVsLow(data,channels)

for i = 1:numel(data)
    tmp(i,:) = mean(data{i}.avg(find(ismember(data{i}.label,channels)),:));    
end

% Calculating row-wise mean
row_means = mean(tmp, 2); % row-wise mean (2 specifies operating along the 2nd dimension)

% For highest 5 mean rows
[sorted_means, sorted_indices] = sort(row_means, 'descend');
high = data(sorted_indices(1:5));

% For lowest 5 mean rows
[sorted_means, sorted_indices] = sort(row_means, 'ascend');
low = data(sorted_indices(1:5));

for i = 1:numel(low)
    if mean(high{i}.avg(find(ismember(high{i}.label,channels)),:)) == mean(low{i}.avg(find(ismember(low{i}.label,channels)),:))
        error('same!')
    end
end

return
end
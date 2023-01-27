function [mn,sem] = returnStats(matrix)

mn = mean(matrix,1);
sem = std(matrix, [], 1)./ sqrt(size(matrix,1));                                % Calculate Standard Error Of The Mean
% ci = bsxfun(@plus, mean(matrix,1), bsxfun(@times, [-1  1]*1.96, sem));   % 95% Confidence Intervals


end
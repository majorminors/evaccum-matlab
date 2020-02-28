function [data2Fit_Free,data2Fit_Spec]=data_stats_indv(dataRaw)
% Get behavioural statistics from indvidual subject
% Used for the Simple RT task in CamCan
% Output:
%   data_to_fit - RT quantiles and other stuff
% Input:
%   dataRaw - the raw informatio from excel result file; Column = ['State' 'Response' 'RT (ms)']

% Data organisation:
% Column(1): button pressed (1 or 2, or 0 for invalid)
% Column(2): accuracy (1 = correct, 0 = incorrect, -1 = invalid)
% Column(3): reaction time (ms)
% Column(4): trial type (1-64)
% Column(5): condition (LcLr,LcHr,HcLr,HcHr)

data=[];
% mappings={{'Index',4,1},...
%             {'Middle',5,2},...
%             {'Ring',6,3},...
%             {'Little',7,4},...
%             {'Free',[4 5 6 7],5}};   %{'state','available response','task cond'}

mappings = {{'LcLr'},...
            {'LcHr'},...
            {'HcLr'},...
            {'HcHr'};

for i=1:size(dataRaw,1)
    for j=1:size(mappings,1)
       if strcmp(dataRaw(i,5),mappings{j}{1})
           if ismember(dataRaw{i,2},mappings{j}{2}) && dataRaw{i,3}>0
               data(end+1,:)=[mappings{j}{3} dataRaw{i,2} dataRaw{i,3}/1000];  % [cond, response, rt in sec]
           else
               disp(['invalid response or rt on trial ' num2str(i-1)]);
           end
       end
    end
end

SpecData=data(find(data(:,1)<5),:); %don't need this
FreeData=data(find(data(:,1)==5),:);

data2Fit_Spec=data_stats(SpecData(:,2:3)); % or this
data2Fit_Free=data_stats_FT4(FreeData(:,[2 3 1]));

end

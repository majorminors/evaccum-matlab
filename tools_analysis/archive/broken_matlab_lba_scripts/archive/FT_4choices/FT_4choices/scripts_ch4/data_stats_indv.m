function [data2Fit_Free,data2Fit_Spec]=data_stats_indv(dataRaw)
% Get behavioural statistics from indvidual subject
% Used for the Simple RT task in CamCan
% Output:
%   data_to_fit - RT quantiles and other stuff
% Input:
%   dataRaw - the raw informatio from excel result file; Column = ['State' 'Response' 'RT (ms)']

% Get data
% ---Response id---
% Index = 4;
% Middle = 5;
% Ring = 6;
% Little = 7;
% ---Task condition---
% Index = 1;
% Middle = 2;
% Ring = 3;
% Little = 4;
% Free =5;
data=[];
mappings={{'Index',4,1},...
            {'Middle',5,2},...
            {'Ring',6,3},...
            {'Little',7,4},...
            {'Free',[4 5 6 7],5}};   %{'state','available response','task cond'}
        
for i=2:size(dataRaw,1)  % first row in xls file is the titles
    for j=1:size(mappings,2)
       if strcmp(dataRaw(i,1),mappings{j}{1})
           if ismember(dataRaw{i,2},mappings{j}{2}) && dataRaw{i,3}>0
               data(end+1,:)=[mappings{j}{3} dataRaw{i,2} dataRaw{i,3}/1000];  % [cond, response, rt in sec]
           else
               disp(['invalid response or rt on trial ' num2str(i-1)]);
           end
       end
    end
end

SpecData=data(find(data(:,1)<5),:);
FreeData=data(find(data(:,1)==5),:);

data2Fit_Spec=data_stats(SpecData(:,2:3));
data2Fit_Free=data_stats_FT4(FreeData(:,[2 3 1]));

end

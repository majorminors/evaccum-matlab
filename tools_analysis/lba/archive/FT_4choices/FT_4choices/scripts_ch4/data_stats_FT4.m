function obs_stats=data_stats_FT4(inpData)
%   Input:
%       inpData: column 1 -response
%       inpData: column 2 -response time
%       inpData: column 3 -condition
%   Output:
%       allObs: observations across responses
%       indvObs: observations for individual response
%       priorProb: probability for individual response and number of trials for individual response


qobs_all=data_stats_simple(inpData(:,2));

resSet=unique(inpData(:,1));

% historical code for 3-choice decisions
ava_cond=[0 0 0 0];  %[8 7 6 5]
num_ava=zeros(1,length(ava_cond));
for i=1:size(inpData,1)
    num_ava=num_ava+1;
    num_ava(ava_cond==inpData(i,3))=num_ava(ava_cond==inpData(i,3))-1;
end
    
    

for i=1:length(resSet)  % loop for statistics of each response
    prob_each{i}=[sum(inpData(:,1)==resSet(i))/num_ava(i) sum(inpData(:,1)==resSet(i))];          % probability for each response
    qobs_each{i}=data_stats_simple(inpData(inpData(:,1)==resSet(i),2));
    rt_median{i}=[median(inpData(inpData(:,1)==resSet(i),2)) sum(inpData(:,1)==resSet(i))];
    prob_noavail(i)=length(find(inpData(:,3)==ava_cond(i)))/size(inpData,1);
end

obs_stats=struct('allObs',qobs_all,'indvObs',{qobs_each},'priorProb',{prob_each},'indvObs_med',{rt_median},'cond_ratio',{prob_noavail});

end


function qobs=data_stats_simple(rt_data)
% Compute the quantile statistics of the experiment data
% Input:
%   rt_data - vector of rt
%   defect_p - defective probability
% output:
%   data_stats:
%        column 1 - quantils
%        column 2 - quantils boundaries
%        column 3 - the probability mass in each bin
%        column 4 - the observed frequencies in each bin

q=[.1; .3; .5; .7; .9;];

if max(size(rt_data))<max(size(q))
    disp('The rt_data set is too small');
    qobs=nan;
    return
end

qobs=quantile(rt_data,q);
qobs(:,2)=q;
qobs(end+1,1:2)=[inf 1]; %Final quantile (corresponding to an infinitely long latency)

%Third column for the probability mass in each bin
qobs(:,3)=qobs(:,2);
qobs(2:end,3)=diff(qobs(:,2));
%Fourth column for the observed frequencies in each bin
qobs(:,4)=qobs(:,3)*length(rt_data);
end

% 
% figure;
% plot(data_to_fit_free_rep{1}.allObs(1:end-1,1),'r')
% hold on
% plot(data_to_fit_free_norep{1}.allObs(1:end-1,1))
% plot(data_to_fit_spec_norep{1}.allObs(1:end-1,1),'g')
% plot(data_to_fit_spec_rep{1}.allObs(1:end-1,1),'k')

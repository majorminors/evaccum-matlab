function [signal,t]=EA_LBA(t01,param,rt_i,n_units,finger,avoid,fsample,shape)
%shape 1; spike
%shape 2; spline
%shape 3; just ramping part
%[pval,pidx]=min(bestval{subj});%best fit 
%parms = bestpar{subj}(pidx,:);%corresponding parameter

%[num_param,param_Low,param_High]=getModelParam_RDK(design,4,parms);


%rt_i = RT{subj}.data{cond}(trial_num,2);finger = RT{subj}.data{cond}(trial_num,3);avoid = RT{subj}.data{cond}(trial_num,4);
if ~exist('shape','var');shape = 1;end
%expected activity
X = {param.Ame, param.Astd, param.B, param.T0,param.C0,rt_i};%mu,sigma,theta,t0,c0,rt,n_units,fsample,b01      
%X = {param_Low.Ame, param_Low.Astd, param_Low.B, param_Low.T0,param_Low.C0,rt_i};%mu,sigma,theta,t0,c0,rt,n_units,fsample,b01      

%NOTE: only 1 accumulator for perceptual condition
[~,~,~,~,signal,t] = expected_activityLBA(X{1},X{2},X{3},X{4},X{5},X{6},n_units,finger,avoid,fsample,t01,shape);


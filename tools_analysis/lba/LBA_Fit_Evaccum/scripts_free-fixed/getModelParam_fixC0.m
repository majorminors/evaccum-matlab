function [num_param,param_Fnrep,param_Snrep]=getModelParam_fixC0_old(mod_feature,param)
% generate model parameter for difference conditions
% total free parameter for a particular design
% if param is null, will return the number of free parameters that satistfy
% the model design
% bssic parameter:
% B, [\mu_1:\mu_4], t0
% feature:
% ------------------ 1 - Std differ across fingers [3 0] not used
% 1 - Std differ across fingers [4 1]
% ------------------ 2 - C0 differ across fingers  [4 1]           not used
% 3 - B differ across conditions [1 0]
% 4 - Drift rate differ across conditions [1 0]
% 5 - t0 differ across conditions [1 0]
% 6 - drift rates are the same for 4 fingers in spec conditions [0 0] not
% used
% 7 - drift rate ratio applied to all fingers rather than only repetition
% not used
% finger [0 0]

N=4; % accumulators
basic_param=6;  %[B \mu_1 : \mu_4, t0]
feature_space=[4 0 1 1 1 0 0; 1 0 0 0 0 0 0]; % [feature ticked; feature unticked]
defaultC0 = 0.5; % upper mean starting point fixed at 0.5

% get index for feature from param vector
feature_len=zeros(1,size(feature_space,2));
for i=1:length(feature_len)
    if ismember(i,mod_feature)
        feature_len(i)=feature_space(1,i);
    else
        feature_len(i)=feature_space(2,i);
    end
end


num_param=sum(feature_len)+basic_param;

if nargin < 2 %if param is not specified then just give number of free params
    return;
end

if num_param~=length(param)
    error('parameter vector length error!');
end

feature_indx=cell(1,size(feature_space,2));
for i=1:length(feature_len)
    if feature_len(i)~=0;
        feature_indx{i}=basic_param+sum(feature_len(1:(i-1)))+(1:feature_len(i)); % C0, B, \mu ratio, t0, \mu in spec trials
    end
end
% 
% if ~isempty(feature_indx{1})
%     Astd=[defaultAstd param(feature_indx{1})];
% else
%     Astd=defaultAstd.*ones(1,N);
% end

param_Fnrep=struct('N',N,'B',{ones(1,N)*param(1)},'C0',defaultC0.*ones(1,N),'Ame',{param(2:5)},'Astd',{param(feature_indx{1}).*ones(1,N)},'T0',param(6));

param_Frep=param_Fnrep;     % free selection, repetition available
param_Snrep=param_Fnrep;    % specific, no repetation



% construct model parameters
% change B
if ~isempty(feature_indx{3})
    param_Snrep.B=param_Snrep.B*param(feature_indx{3}(1));
end
% param_Frep.alpha=1;

% change mu
if ~isempty(feature_indx{4})
    param_Snrep.Ame=param_Snrep.Ame*param(feature_indx{4}(1));
end
% param_Frep.beta=1;

% change t0
if ~isempty(feature_indx{5})
    param_Snrep.T0=param_Snrep.T0*param(feature_indx{5}(1));
end
% param_Frep.theta=1;



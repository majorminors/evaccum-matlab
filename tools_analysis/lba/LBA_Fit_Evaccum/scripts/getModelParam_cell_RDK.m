function [num_param,paramLL,paramHL,paramLH,paramHH]=getModelParam_cell_RDK(mod_feature,N,param)
% generate model parameter for difference conditions
% total free parameter for a particular design
% if param is null, will return the number of free parameters that satistfy
% the model design
% basic parameter:
% B, [\mu_1:\mu_4], t0s
% feature:
% 
% 1 - Std differ across fingers [4 1]
% 2 - C0 differ across fingers  [4 1]          
% 3 - B differ across conditions [1 0]
% 4 - Drift rate differ across conditions [1 0]
% 5 - t0 differ across conditions [1 0]
% 6 - drift rates are the same for 4 fingers in spec conditions [0 0] not
% used
% 7 - drift rate ratio applied to all fingers rather than only repetition
% not used
% finger [0 0]


basic_param=N+2;  %[B \mu_1 : \mu_4, t0] - mu = response options (accumulators reqd)


%note 1 in first row: 2 were used to create different starting params for
%rep and no-rep conditions.

varC0 =ismember(2,mod_feature);

% so the below (at 3) is saying we need 3 additional parameters above the
% one we will generate
if varC0
    feature_space=[N N 3 3 3 0 0; 1 0 0 0 0 0 0]; % [feature ticked; feature unticked]
else
    feature_space=[N 0 3 3 3 0 0; 1 0 0 0 0 0 0];
    defaultC0 = 0.1;
end
% get indexfor feature from param vector: mod_feature selects what feature
% to use e.g[1 4] = std differ across fingers and Drift rate differs across
% conditions

feature_len=zeros(1,size(feature_space,2));% flag used params
for i=1:length(feature_len)
    if ismember(i,mod_feature)
        feature_len(i)=feature_space(1,i);%use it
    else
        feature_len(i)=feature_space(2,i);%don't use it
    end
end



num_param=sum(feature_len)+basic_param;%add to basic params

if nargin < 3
    return;
end

if num_param~=length(param)
    error('parameter vector length error!');
end

feature_indx=cell(1,size(feature_space,2));
for i=1:length(feature_len)
    if feature_len(i)~=0
        feature_indx{i}=basic_param+sum(feature_len(1:(i-1)))+(1:feature_len(i)); % C0, B, \mu ratio, t0, \mu in spec trials
    end
end
% 
% if ~isempty(feature_indx{1})
%     Astd=[defaultAstd param(feature_indx{1})];
% else
%     Astd=defaultAstd.*ones(1,N);
% end

% this part gives you the starting parameters for the fitting, so we set
% them for condition low/high (LH) (arbitrarily) and then assign those same
% values to the other three conditions so all have the same starting
% parameters
if varC0
    paramLH=struct('N',N,'B',{ones(1,N)*param(1)},'C0',param(feature_indx{2}).*ones(1,N),'Ame',{param(2:(basic_param-1))},'Astd',{param(feature_indx{1}).*ones(1,N)},'T0',param(basic_param));
else
    
    paramLH=struct('N',N,'B',{ones(1,N)*param(1)},'C0',defaultC0.*ones(1,N),'Ame',{param(2:(basic_param-1))},'Astd',{param(feature_indx{1}).*ones(1,N)},'T0',param(basic_param));
end


paramLL=paramLH;
paramHL=paramLH;
paramHH=paramLH;


% construct model parameters
% change B
if ~isempty(feature_indx{3})
    paramLL.B=paramLL.B*param(feature_indx{3}(1));%note: Low is a proportion of High
    paramHL.B=paramHL.B*param(feature_indx{3}(2));%note: Low is a proportion of High
    paramHH.B=paramHH.B*param(feature_indx{3}(3));%note: Low is a proportion of High

end
% param_Frep.alpha=1;

% change mu
if ~isempty(feature_indx{4})
    paramLL.Ame=paramLL.Ame*param(feature_indx{4}(1));
    paramHL.Ame=paramHL.Ame*param(feature_indx{4}(2));
    paramHH.Ame=paramHH.Ame*param(feature_indx{4}(3));
    
end
% param_Frep.beta=1;

% change t0
if ~isempty(feature_indx{5})
    paramLL.T0=paramLL.T0*param(feature_indx{5}(1));
    paramHL.T0=paramHL.T0*param(feature_indx{5}(2));
    paramHH.T0=paramHH.T0*param(feature_indx{5}(3));
    
end
% param_Frep.theta=1;



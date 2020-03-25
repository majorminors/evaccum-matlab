function [f_error, BIC, pmod_LL,pmod_LH,priorMod_HA,qobs] = fiterror_cell_EASYCOH(param, mod_feature, data2fit)
% fit the model to freechoice and specified selection
% param [boundary, X0 range, drift 1, drift 2, drift 3, drift 4, drift std,
% T0, theta,alpha] 

if length(data2fit)==2
    goalstat_LL=data2fit{1};
    goalstat_LH=data2fit{2};   
end

minRT = goalstat_LL.minRT;

Ntrials= totalTrials(data2fit);

% compute the statistics for model
%N=1;
% B_num=0;    % 1. fixed boundary for all choice
BIC=0;
f_error=0;
% model_param=struct('N',N,'B',{param(1:B_num).*ones(1,N)},'C0',param(B_num+1),'Ame',{[param(B_num+2:(N+B_num)+1) ]},'Astd',param(N+B_num+2),'T0',param(N+B_num+3));

[num_param,parLL,~,parLH,~]=getModelParam_cell_RDK(mod_feature,2,param);

qobs=cell(1,4);

if sum(param<0)>0
    f_error=1000000;  % negative parameter
    BIC=1000000;
    return;
elseif sum(parLL.B<=parLL.C0) || sum(parLH.B<=parLH.C0)
    f_error=1000000;  % Threshold is smaller than starting point
    BIC=1000000;
    return;
    
else    
%%     
    % Choice: Low PU
    pmod_LH = nan; priorMod_LH = nan;%#ok
    [pmod_LH,priorMod_LH,qobs{3}]=mod_stats_sim('LBA_spec_FT3_c_even_template',parLH,1,goalstat_LH.allObs(:,1),2,0); %,goalstat_free_norep.cond_ratio);
    
    f_error=f_error+2*sum(goalstat_LH.allObs(:,4).*log(goalstat_LH.allObs(:,3)./pmod_LH'));
    BIC=BIC-2*sum(goalstat_LH.allObs(:,4).*log(pmod_LH'));
    for i=1:parLH.N
        f_error=f_error+2*sum(goalstat_LH.priorProb{i}(2).*log(goalstat_LH.priorProb{i}(1)./priorMod_LH(i)));
        BIC=BIC-2*sum(goalstat_LH.priorProb{i}(2).*log(priorMod_LH(i)));
    end
    
    
    

    %%
    % Specified: Low PU 
    pmod_LL =nan;priorMod_LL=nan;%#ok
    [pmod_LL,foo1,qobs{1}]=mod_stats_sim('LBA_spec_FT3_c_even_template',parLL,1,goalstat_LL.allObs(:,1),2,0); %#ok ,goalstat_spec_norep.cond_ratio);
    
    f_error=f_error+2*sum(goalstat_LL.allObs(:,4).*log(goalstat_LL.allObs(:,3)./pmod_LL'));
    BIC=BIC-2*sum(goalstat_LL.allObs(:,4).*log(pmod_LL'));    
        
%%    
    %error is zero unless special conditions hold, as above;
    %gaolstat_H.allObs refers to the prob of the RT quanitles for the
    %observed data;  and pmod_H refers to the probability of the
    %quanitles from th simulation.  with reference to Ludwig et al
    %formula (1) for G^2 statistic.  ,4 refers to number of trials per
    %condition.   ,3 refers to distance between the quanitle boundaries. 
    
    
    %not yet real BIC, but start to build BIC: see line 64
    %and again for the p distribution of choices, again wrt Ludwig et al
    %g^2 formula
    
    BIC=BIC+length(param)*log(Ntrials); % penalty for num. of model parameters
end
end


function Ntrials=totalTrials(data2fit)
Ntrials=0;
% Nfreetrials=0;
for i=1:length(data2fit)
    Ntrials=Ntrials+sum(data2fit{i}.allObs(:,4));    
%     if i>2
%         Nfreetrials=Nfreetrials+sum(data2fit{i}.allObs(:,4));  
%     end
end
end



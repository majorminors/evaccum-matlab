function [error, BIC, pmod_FT, priorMod_FT, pmod_FT_rep,priorMod_FT_rep,pmod_spec,pmod_spec_norep,qobs] = fiterror_cs_overall_FT4_template (param, mod_feature, data2fit)
% fit the model to freechoice and specified selection
% param [boundary, X0 range, drift 1, drift 2, drift 3, drift 4, drift std,
% T0, theta,alpha] 

if length(data2fit)==2
    goalstat_free=data2fit{1};
    goalstat_spec=data2fit{2};
else
    error('datafile in wrong format');
end

Ntrials=totalTrials(data2fit);

% compute the statistics for model
N=4;
% B_num=0;    % 1. fixed boundary for all choice
BIC=0;
error=0;
% model_param=struct('N',N,'B',{param(1:B_num).*ones(1,N)},'C0',param(B_num+1),'Ame',{[param(B_num+2:(N+B_num)+1) ]},'Astd',param(N+B_num+2),'T0',param(N+B_num+3));
[num_par,par_Fnrep,par_Snrep]=getModelParam_fixC0(mod_feature,param);

pmod_FT=nan;
priorMod_FT=nan;
% pmod_FT_rep=nan;
% priorMod_FT_rep=nan;
pmod_spec=nan;
% pmod_spec_norep=nan;
qobs=cell(1,2);

if sum(param<0)>0
    error=1000000;  % negative parameter
    BIC=1000000;
    return;
elseif sum(par_Fnrep.B<=par_Fnrep.C0) || sum(par_Snrep.B<=par_Snrep.C0)
    error=1000000;  % Threshold is smaller than starting point
    BIC=1000000;
    return;
else
    % Free-4 selection, repetition not available
    [pmod_FT,priorMod_FT,qobs{1}]=mod_stats_sim('LBA_spec_FT3_c_even_template',par_Fnrep,1,goalstat_free.allObs(:,1),4,0); %,goalstat_free_norep.cond_ratio);
    %for RT
    %error is zero unless special conditions hold, as above;
    %gaolstat_free.allObs referes to the prob of the RT quanitles for the
    %observed data;  and pmod_FT refers to th eprobabiolity of the
    %quanitles from th esimulation.  with reference to Ludwig et al
    %formula (1) for G^2 statistic.  ,4 refers to number of trials per
    %condition.   ,3 refers to distnace between the quanitle boundaries. 
    error=error+2*sum(goalstat_free.allObs(:,4).*log(goalstat_free.allObs(:,3)./pmod_FT'));
    BIC=BIC-2*sum(goalstat_free.allObs(:,4).*log(pmod_FT')); %not yet real BIC, but start to build BIC: see line 64
    %and again for the p distribution of choices, again wrt Ludwig et al
    %g^2 formula
    for i=1:N
        error=error+2*sum(goalstat_free.priorProb{i}(2).*log(goalstat_free.priorProb{i}(1)./priorMod_FT(i)));
        BIC=BIC-2*sum(goalstat_free.priorProb{i}(2).*log(priorMod_FT(i)));
    end

    %specified trials
    if exist('goalstat_spec','var')
        [pmod_spec,foo2,qobs{2}]=mod_stats_sim('LBA_spec_c_even_template',par_Snrep,1,goalstat_spec.allObs(:,1),4,0); %,goalstat_spec_norep.cond_ratio);
        error=error+2*sum(goalstat_spec.allObs(:,4).*log(goalstat_spec.allObs(:,3)./pmod_spec'));
        BIC=BIC-2*sum(goalstat_spec.allObs(:,4).*log(pmod_spec'));
    end
    
    BIC=BIC+length(param)*log(Ntrials); % panalty for num. of model parameters, now true BIC
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



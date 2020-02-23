function [pmod,probmod,qobs,rtquantile,cumRT]=mod_stats_sim(model,model_param,posRTonly,quantBound,numAct,getcumRT,cond_ratio)
%   numAct - number of possible action (4: 4 selection, 3: 3 selection)
%   quantBound -  number of qunatiles up to Inf

probmod=nan(1,length(model_param.N));
% quantBound=quantBound*1000;
% UpperRT=3000;

if exist('cond_ratio','var')
    [rt_iter,act_iter,cond]= feval (model, model_param.N,model_param.B,model_param.C0,model_param.Ame,model_param.Astd,model_param.T0,cond_ratio);
else
    [rt_iter,act_iter]= feval (model, model_param.N,model_param.B,model_param.C0,model_param.Ame,model_param.Astd,model_param.T0,ones(1,numAct)./numAct);
end
if posRTonly==1     % only use positive RT
    if sum(rt_iter<=0 | act_iter==0)~=0
        indx=find(rt_iter<=0 | act_iter==0 );
%         disp([num2str(length(indx)) ' negative RT trials!']);
        act_iter(indx)=[];
        rt_iter(indx)=[];
%         cond(indx)=[];
    end
end

q=[.1; .3; .5; .7; .9;];
rtquantile=zeros(length(q),1);
rtquantile=quantile(rt_iter,q);
rtquantile(rtquantile==0)=eps;


% rt_iter=rt_iter.*1000;
nmodel=histc(rt_iter,[-Inf;quantBound]); %count the number of simulated latencies in between the quantile boundaries
pmod=nmodel(1:end-1)/length(rt_iter); %convert counts into probabilities, ignoring the last bin as it will have a count of 0 (corresponding to inf)
pmod(pmod==0)=eps; %avoid numerical problems with predicted probabilities of 0

for i=1:model_param.N
    if exist('cond','var')
        if strcmp(model,'LBA_spec_c_even_template')
            probmod(i)=sum(act_iter==i)/length(find(cond==(i-1)));
        elseif strcmp(model,'LBA_spec_FT3_c_even_template')
            probmod(i)=sum(act_iter==i)/length(find(cond~=(i-1)));
        else
            probmod(i)=sum(act_iter==i)/(length(cond)-length(find(cond==(i-1))));
        end
    else
        probmod(i)=sum(act_iter==i)/(length(act_iter)*numAct/model_param.N);
    end
end
probmod(probmod==0)=eps;

if exist('getcumRT','var')     % get cumulative RT distribtion from simulation results
    
    if getcumRT==1
        for i=1:length(pmod)
            cumRT(i)=sum(pmod(1:i));
        end
%         nmodel_step=histc(rt_iter,0:UpperRT);
%         pmod_step=nmodel_step(1:end-1)/length(rt_iter);
%         for i=1:UpperRT-1
%             cumRT(i)=sum(pmod_step(1:i));
%         end

    end
end

qobs=data_stats_simple(rt_iter);
end

function [F_all,F]=LBA_rt_pdf(params,rt)
N=params.N; % Num of choice
B=params.B; % boundary/threshold
C0=params.C0; % X0 range [0 C_X0]
Ame=params.Ame;  % mean drift rate for the N choices
Astd=params.Astd; % SAME stand diviation of the drift rate for the N choice
T0=params.T0;     % non-decision time
F_all=0;
t=rt-T0;
for i=1:N
    F(i)=1+((B(i)-C0-t*Ame(i))/C0)*normcdf((B(i)-C0-t*Ame(i))/(t*Astd))...
        -((B(i)-t*Ame(i))/C0)*normcdf((B(i)-t*Ame(i))/(t*Astd))...
        +(t*Astd/C0)*normpdf((B(i)-C0-t*Ame(i))/(t*Astd))...
        -(t*Astd/C0)*normpdf((B(i)-t*Ame(i))/(t*Astd));
    F_all=F_all+F(i)*1/N;
end
end

function qobs=data_stats_simple(rt_data,defect_p)
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

if exist('defect_p','var')
%     qobs(:,2)=qobs(:,2)*defect_p ; %normalise to the defective probability
    
    
    %if the distribution does not add to 1, create an additional bin for
    %the remaining alternative responses (e.g. errors)
%     if(defect_p<1)
%         qobs(end+1,1:2)=[NaN 1];
%     end
end
%Third column for the probability mass in each bin
qobs(:,3)=qobs(:,2);
qobs(2:end,3)=diff(qobs(:,2));
%Fourth column for the observed frequencies in each bin
qobs(:,4)=qobs(:,3)*length(rt_data);
end

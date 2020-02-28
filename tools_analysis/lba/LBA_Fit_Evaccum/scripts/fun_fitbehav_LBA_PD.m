function fun_fitbehav_LBA_PD(settings,data)
% Fit LBA model to RDK-tapping task
% Use for  3-choice task: CONDITIONS segregated
% JZ 15/10/2014 - 4-choice tapping task
% AT 17/02/2016 adapted to 3-choice RDK task
% DM 02/2020 adapted to 2-choice RDK task

rseed = settings.rseed;
rng(rseed,'twister') % for reproducibility

%clear all
%%
savename = settings.savename;
Model_Feature = settings.modfeat; % model variant for this job

% optional but doesn't save time (thought it would be faster than normal
% gradient descent, but isn't)
% if ~isfield(settings,'bayesOptim')
%     bayesian = 0;
% else
%     bayesian = settings.bayesOptim;
%     if bayesian
%         addpath(genpath('/imaging/at07/Matlab/Toolboxes/bads'));
%     end
% end


% experimental designs
%contrasts(:,:,1) = [1 3; 2 4];%main effect P
%contrasts(:,:,2) = [1 2; 3 4];%main effect A

%contrasts = 1:4;

num_subjs = sum(~cellfun('isempty',{data.id})); % get the number of not empty arrays from the first field of the 'data' structure

%initialise cells
data2fit={};Ava_opt = {}; CohLevs = {};%#ok

for idxsubj = 1:length(num_subjs)
    %%
    clear rt conds resp acc % clear vars from behavioural data to avoid reading mixed subject data
    fprintf('loading %s \n',data.id(idxsubj))
    
    error('youre here, trying to get condLabs to pull the condition labels from subjdata');
    % pull the data
    subjdata = data(idxsubj).data; % put the data in a readable variable
    dataValid = []; % now we'll strip invalid responses out
    condLabs = [];
    for i = 1:size(subjdata,1)
        if subjdata{i,5} >= 0 % if there's a valid response
            dataValid(end+1,:) = [subjdata{i,2} subjdata{i,3} subjdata{i,4} subjdata{i,5}]; % add the following rows to this new variable in order: condition code, response, rt, accuracy
            condLabs{end+1,:} = subjdata{i,1}; % add the valid condition labels
        end
    end; clear i;
    
    % converting again to readable variables
    conds = dataValid(:,1);
    resp = dataValid(:,2);
    rt = dataValid(:,3);
    acc = dataValid(:,4);    
    
    % do descriptive stats    
    minRT(idxsubj) = min(rt(rt>=0.1)); % not used
    

    %%
    %Pool conditions according to design matrix & calculate stats
    for level = 1:4 % four conditions
        
        trialn = conds == level;
        dataRaw{idxsubj,level}  = [resp(trialn) rt(trialn)];
        
        data2fit{idxsubj,level} = data_stats(dataRaw{idxsubj,level});%#ok

        % add some info to data2fit (the condition labels as a string, and
        % the minRT - in the case that the min rt is smaller than the
        % non-decision time - let's leave this for now, just in case)
        data2fit{idxsubj,level}.cond = condLabs{trialn};%#ok
        data2fit{idxsubj,level}.minRT = round(minRT(idxsubj),2);
        %parange(end,1) = round(minRT(idxsubj),2);%constrain the lower bound of T0 to the shortest RT
    end
    %%
end



% fit the basic model
% startpar=[1,1,1,1,1, .1 .3]; % initial parameter [boundary,four drift rate, non-decisiontime, drift rate std]
% upper staring-point is fixed
%randiter = 500;%100; % randome search 100 iters before optimization
%nosession = 100;%100;%20; % 20 optimization iterations

randiter  = settings.randiter;
nosession = settings.nosession;


% Model fitting

% model features that are available
% 1- std differs across fingers
% 2- C0  differs across fingers
% 3- B differs across conditions
% 4- Drift rate differs across conditions
% 5- t0 differs across conditions

%design_space={};
% design_space={[1,3],[1,4],[1,3,4],[1,3,4,5],[1,2],[1,2,3],[1,2,4],[1,2,3,4],[1,2,3,4,5]};
% Model_Feature=[1  4];
%
% numParam   =getModelParam_cell_RDK(Model_Feature,4);
%
% parange=[zeros(1,numParam);zeros(1,numParam)+10]';


%%


numParam   =getModelParam_cell_RDK(Model_Feature,2); % you will change this 4 to two - reflects response options
parange=[zeros(1,numParam);zeros(1,numParam)+10]';


for idxsubj = 1:length(num_subjs)
    
    [bestpar{idxsubj,1},bestval{idxsubj,1},BIC{idxsubj,1}]=fitparams_refine_template_RDK('fiterror_cell_RDK',Model_Feature,{data2fit{idxsubj,1:4}},randiter,nosession,[],parange,bayesian);%#ok
    
end

save(savename,'bestpar','bestval','BIC','rseed','settings');%,'bestpar_PA','bestval_PA','BIC_PA')

end

%
%
% for idxsubj = 1:length(sname)    % Best mdoel
%     [fooL(idxsubj),best_indxL(idxsubj)]=min(bestval_Low{idxsubj,1});%#ok
%     disp(['best parameter vector : ' num2str(bestpar_Low{idxsubj,1}(best_indxL(idxsubj),:))]);
%     disp(['best error funcion value: ' num2str(bestval_Low{idxsubj,1}(best_indxL(idxsubj)))]);
%
%     % Re-run the best model to get simulation results
%
%     [num_parL(idxsubj,1),par_FreeL(idxsubj,1),par_SpecL(idxsubj,1)]=getModelParam_fixC0(Model_Feature,bestpar_Low{idxsubj,1}(best_indxL(idxsubj),:));%#ok
%     % [pmod_simpRT,qobs{1},cumRT,propInvalidTrial]=mod_stats_sim_simpRT('LBA_simpRT_template',par_simpRT,1,data2fit.allObs(:,1),1,1); %,goalstat_spec_norep.cond_ratio);
%     [pmod_FT_L{idxsubj,1},priorMod_FT_L{idxsubj,1},qobs_Free_L{idxsubj, 1},foo3L{idxsubj,1},cumRT_Free_L{idxsubj,1}]=mod_stats_sim('LBA_spec_FT3_c_even_template',par_FreeL(idxsubj,1),1,data2fit{idxsubj, 1}.allObs(:,1),3,1);%#ok %,goalstat_free_norep.cond_ratio);
%     [pmod_spec_norep_L{idxsubj,1},foo2L{idxsubj,1},qobs_Spec_L{idxsubj,1},foo4L{idxsubj,1},cumRT_Spec_L{idxsubj,1}]=mod_stats_sim('LBA_spec_c_even_template',par_SpecL(idxsubj,1),1,data2fit{idxsubj, 3}.allObs(:,1),4,1);%#ok %,goalstat_spec_norep.cond_ratio);
%
%
%     [fooH(idxsubj),best_indxH(idxsubj)]=min(bestval_High{idxsubj,1});%#ok
%     disp(['best parameter vector : ' num2str(bestpar_High{idxsubj,1}(best_indxH(idxsubj),:))]);
%     disp(['best error funcion value: ' num2str(bestval_High{idxsubj,1}(best_indxH(idxsubj)))]);
%
%     % Re-run the best model to get simulation results
%
%     [num_parH(idxsubj,1),par_FreeH(idxsubj,1),par_SpecH(idxsubj,1)]=getModelParam_fixC0(Model_Feature,bestpar_High{idxsubj,1}(best_indxH(idxsubj),:));%#ok
%     % [pmod_simpRT,qobs{1},cumRT,propInvalidTrial]=mod_stats_sim_simpRT('LBA_simpRT_template',par_simpRT,1,data2fit.allObs(:,1),1,1); %,goalstat_spec_norep.cond_ratio);
%     [pmod_FT_H{idxsubj,1},priorMod_FT_H{idxsubj,1},qobs_Free_H{idxsubj, 1},foo3H{idxsubj,1},cumRT_Free_H{idxsubj,1}]=mod_stats_sim('LBA_spec_FT3_c_even_template',par_FreeH(idxsubj,1),1,data2fit{idxsubj, 2}.allObs(:,1),3,1);%#ok %,goalstat_free_norep.cond_ratio);
%     [pmod_spec_norep_H{idxsubj,1},foo2H{idxsubj,1},qobs_Spec_H{idxsubj,1},foo4H{idxsubj,1},cumRT_Spec_H{idxsubj,1}]=mod_stats_sim('LBA_spec_c_even_template',par_SpecH(idxsubj,1),1,data2fit{idxsubj, 4}.allObs(:,1),4,1);%#ok %,goalstat_spec_norep.cond_ratio);
%
%
% end
%
% for s = 1:numel(sname)
%
%     driftSpec.L(s,:,1)=par_SpecL(s,1).Ame;%mean
%     driftSpec.L(s,:,2)=par_SpecL(s,1).Astd;%std
%
%     boundSpec.L(s,:,1)=par_SpecL(s,1).B;%mean
%     TSpec.L(s,:,1)=par_SpecL(s,1).T0;%mean
%
%     driftSpec.H(s,:,1)=par_SpecH(s,1).Ame;%mean
%     driftSpec.H(s,:,2)=par_SpecH(s,1).Astd;%std
%
%     boundSpec.H(s,:,1)=par_SpecH(s,1).B;%mean
%     TSpec.H(s,:,1)=par_SpecH(s,1).T0;%mean
%
%     driftFree.L(s,:,1)=par_FreeL(s,1).Ame;%mean
%     driftFree.L(s,:,2)=par_FreeL(s,1).Astd;%std
%
%     boundFree.L(s,:,1)=par_FreeL(s,1).B;%mean
%     TFree.L(s,:,1)=par_FreeL(s,1).T0;%mean
%
%     driftFree.H(s,:,1)=par_FreeH(s,1).Ame;%mean
%     driftFree.H(s,:,2)=par_FreeH(s,1).Astd;%std
%
%     boundFree.H(s,:,1)=par_FreeH(s,1).B;%mean
%     TFree.H(s,:,1)=par_FreeH(s,1).T0;%mean
% end
%
%
% save('RDK_fits_130516_3_4N_A.mat')
%
% %% Plotting...
% figure;
% for idxsubj = 1: length(sname)
%     subplot(4,4,idxsubj)
%     plot(data2fit{idxsubj,1}.allObs(1:end-1,1),data2fit{idxsubj,1}.allObs(1:end-1,2),'k','Marker','o'); hold on;
%     plot(data2fit{idxsubj,1}.allObs(1:end-1,1),cumRT_Free_L{idxsubj,1}(1:end-1),'r','Marker','*');
%     legend({'data','model'},'Location','SouthEast');
%     % title({['Model params (B, A, Astd,To): [' num2str(bestpar(best_indx,:),2) ']'] });
%     title('Chosen condition');
%     ylim([0 1]);
%     xlim([0.3 1.2]);
%     % xlim([0.2,0.6]);
%     ylabel('Cumulative probability');
%     xlabel('Response Time (s)');
% end
%
% figure;
% for idxsubj = 1: length(sname)
%     subplot(4,4,idxsubj)
%     plot(data2fit{idxsubj,3}.allObs(1:end-1,1),data2fit{idxsubj,3}.allObs(1:end-1,2),'k','Marker','o'); hold on;
%     plot(data2fit{idxsubj,3}.allObs(1:end-1,1),cumRT_Spec_L{idxsubj,1}(1:end-1),'r','Marker','*');
%     legend({'data','model'},'Location','SouthEast');
%     % title({['Model params (B, A, Astd,To): [' num2str(bestpar(best_indx,:),2) ']'] });
%     title('Specified condition');
%     ylim([0 1]);
%     xlim([0.3 1.2]);
%     % xlim([0.2,0.6]);
%     ylabel('Cumulative probability');
%     xlabel('Response Time (s)');
% end
%
%
% figure;
% for idxsubj = 1: length(sname)
%     subplot(4,4,idxsubj)
%     plot(data2fit{idxsubj,2}.allObs(1:end-1,1),data2fit{idxsubj,2}.allObs(1:end-1,2),'k','Marker','o'); hold on;
%     plot(data2fit{idxsubj,2}.allObs(1:end-1,1),cumRT_Free_H{idxsubj,1}(1:end-1),'r','Marker','*');
%     legend({'data','model'},'Location','SouthEast');
%     % title({['Model params (B, A, Astd,To): [' num2str(bestpar(best_indx,:),2) ']'] });
%     title('Chosen condition');
%     ylim([0 1]);
%     %xlim([0.3 1.2]);
%     % xlim([0.2,0.6]);
%     ylabel('Cumulative probability');
%     xlabel('Response Time (s)');
%     set(gcf,'Color','White')
%
% end
%
% figure;
% for idxsubj = 1: length(sname)
%     subplot(4,4,idxsubj)
%     plot(data2fit{idxsubj,4}.allObs(1:end-1,1),data2fit{idxsubj,4}.allObs(1:end-1,2),'k','Marker','o'); hold on;
%     plot(data2fit{idxsubj,4}.allObs(1:end-1,1),cumRT_Spec_H{idxsubj,1}(1:end-1),'r','Marker','*');
%     legend({'data','model'},'Location','SouthEast');
%     % title({['Model params (B, A, Astd,To): [' num2str(bestpar(best_indx,:),2) ']'] });
%     title('Specified condition');
%     ylim([0 1]);
%     xlim([0.3 1.2]);
%     % xlim([0.2,0.6]);
%     ylabel('Cumulative probability');
%     xlabel('Response Time (s)');
%     set(gcf,'Color','White')
% end
%
%
%
% figure; %hold on;
% for idxsubj = 1: length(sname)
%     temp=[data2fit{idxsubj,1}.priorProb{1:4}];
%     choiceProb=[temp(1:2:8);priorMod_FT_L{idxsubj,1}];
%     subplot(4,4,idxsubj)
%     bar(choiceProb');
%      ylim([0 0.6]);
%     set(gca,'XTick',1:4);
%     set(gca,'XTickLabel',{'Index' 'Middle' 'Ring' 'Little'});
%     title('Choice probability in Free condition');
%     set(gcf,'Color','White')
% end
% legend({'data','model'});
%
% figure; %hold on;
% for idxsubj = 1: length(sname)
%     temp=[data2fit{idxsubj,2}.priorProb{1:4}];
%     choiceProb=[temp(1:2:8);priorMod_FT_H{idxsubj,1}];
%     subplot(4,4,idxsubj)
%     bar(choiceProb');
%      ylim([0 0.6]);
%     set(gca,'XTick',1:4);
%     set(gca,'XTickLabel',{'Index' 'Middle' 'Ring' 'Little'});
%     title('Choice probability in Free condition');
%     set(gcf,'Color','White')
% end
% legend({'data','model'});
%
%
% %Plot estimated drifts
% MatDDMPar = [];
% for idxsubj = 1: length(sname)
%
%     MattDDMPar(idxsubj,1,1) = mean(par_FreeL(idxsubj,1).Ame);%#ok
%     MattDDMPar(idxsubj,2,1) = mean(par_FreeH(idxsubj,1).Ame);%#ok
%     MattDDMPar(idxsubj,3,1) = mean(par_SpecL(idxsubj,1).Ame);%#ok
%     MattDDMPar(idxsubj,4,1) = mean(par_SpecH(idxsubj,1).Ame);%#ok
%
%     MattDDMPar(idxsubj,1,2) = mean(par_FreeL(idxsubj,1).Astd);%#ok
%     MattDDMPar(idxsubj,2,2) = mean(par_FreeH(idxsubj,1).Astd);%#ok
%     MattDDMPar(idxsubj,3,2) = mean(par_SpecL(idxsubj,1).Astd);%#ok
%     MattDDMPar(idxsubj,4,2) = mean(par_SpecH(idxsubj,1).Astd);%#ok
%
% end
%
% doh
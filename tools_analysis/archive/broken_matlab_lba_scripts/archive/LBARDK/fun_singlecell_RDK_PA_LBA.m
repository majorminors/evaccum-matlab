function fun_singlecell_RDK_PA_LBA(settings)
% Fit LBA model to RDK-tapping task
% Use for  3-choice task: CONDITIONS segregated
% JZ 15/10/2014 - 4-choice tapping task
% AT 17/02/2016 adapted to 3-choice RDK task

rng(17,'twister') % for reproducibility

%clear all
%%
droot = '/imaging/at07/Matlab/Projects/CBU2015/RDKUnc/MEGData/';
datafld = '/imaging/at07/Matlab/Projects/CBU2015/RDKUnc/Behav/ModelData/';
stemf = '%s_LBAbehav.mat';
savename = settings.savename;

addpath(genpath('/imaging/at07/Matlab/Projects/CBU2015/RDKUnc/Process_code/'));

Model_Feature = settings.modfeat;

condLabs = {'pLaL','pHaL','pLaH','pHaH'};




% experimental designs
%contrasts(:,:,1) = [1 3; 2 4];%main effect P
%contrasts(:,:,2) = [1 2; 3 4];%main effect A

%contrasts = 1:4;


%% Get sample IDs
% --------------------
sname = at_getsubjs(droot);

%initialise cells
data2fit={};Ava_opt = {}; dataRaw = {}; CohLevs = {};%#ok

for idxsubj = 1:length(sname)   
    
    clear ResponseArray LocStim muRT rt conds badT resp loc
    
    load(sprintf([datafld stemf],sname{idxsubj}))
    fprintf(['loading ' stemf '\n'],sname{idxsubj})
    



% do descriptive stats

rt    = muRT(1).megRT(:);
loc   = muRT(1).LocStim;
if size(rt,1)~=size(loc,1);
    fprintf('%s: correcting number of loc \n',sname{idxsubj});
    loc = loc(1:size(rt,1),:);
end
conds = muRT(1).cond;
badT  = muRT(1).badtrials;
badT  = unique([badT;find(rt==0)]);
resp  = muRT(1).finger; %note fingers inverted in MEG: 4 = index


%remove bad trials
if ~isempty(badT)
    disp([sname{idxsubj}, ' invalid trials (',num2str(length(badT)),'/',num2str(length(rt)), '):',num2str(badT')]);
    
    resp(badT) = [];rt(badT) = []; conds(badT)=[];loc(badT,:)=[];
    
end

%check responses
check = [];for i = 1:size(loc,1);check(i) =loc(i,resp(i));end%#ok
if sum(check)~=size(loc,1); error('stim location and response do not correspond!');end


%Pool conditions according to design matrix & calculate stats
  for levels = 1:4
       
        
        trialn = strcmp(conds,condLabs(levels));%identify condition
        dataRaw{idxsubj,levels}  = [resp(trialn) rt(trialn)]; %#ok
        
        if  sum(ismember(levels,[1 2]));%single stim
            data2fit{idxsubj,levels} = data_stats(dataRaw{idxsubj,levels});%#ok
        else
            tmp = loc(trialn,:);tmp(:,5) = nan;for i = 1:size(tmp,1);tmp(i,5)=find(tmp(i,:)==0);end
            dataRaw{idxsubj,levels}(:,3) = tmp(:,5);%#ok
            data2fit{idxsubj,levels}     = data_stats_FT4RDK(dataRaw{idxsubj,levels});%#ok
             
        end
        data2fit{idxsubj,levels}.cond= condLabs(levels);%#ok
   end

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


numParam            =getModelParam_PA_cell_RDK(Model_Feature,4);
parange             =[zeros(1,numParam);zeros(1,numParam)+10]';
  [~,~,~,~,~,ratios]  =getModelParam_PA_cell_RDK(Model_Feature,4,parange(:,1)');
  parange(ratios ,:) = repmat([0.01 0.99],[numel(ratios),1]);%P ratio needs to be constrainted between 0 and 1;
%  parange(ratios+1 ,:) = repmat([0.01 0.99],[numel(ratios),1]);%P ratio needs to be constrainted between 0 and 1;

for idxsubj = 1:length(sname)
          
    [bestpar{idxsubj,1},bestval{idxsubj,1},BIC{idxsubj,1}]=fitparams_refine_template_RDK('fiterror_cell_PA_RDK',Model_Feature,{data2fit{idxsubj,1:4}},randiter,nosession,[],parange);%#ok
   
end  

save(savename,'bestpar','bestval','BIC','settings');%,'bestpar_PA','bestval_PA','BIC_PA')

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
%      [~,~,~,~,cumRTLL]=mod_stats_sim('LBA_spec_c_even_template',paramLL,1,data2fit{1,1}.allObs(:,1),1,1);
%      [~,~,~,~,cumRTHL]=mod_stats_sim('LBA_spec_c_even_template',paramHL,1,data2fit{1,2}.allObs(:,1),1,1);
%      
%      [~,~,~,~,cumRTLH]=mod_stats_sim('LBA_spec_FT3_c_even_template',paramLH,1,data2fit{1,3}.allObs(:,1),3,1);
%      [~,~,~,~,cumRTHH]=mod_stats_sim('LBA_spec_FT3_c_even_template',paramHH,1,data2fit{1,4}.allObs(:,1),3,1);     
% 
%      figure;
%      subplot(2,2,1)
%      plot(data2fit{1,1}.allObs(1:end-1,1),data2fit{1,1}.allObs(1:end-1,2),'ko');hold on;plot(data2fit{1,1}.allObs(1:end-1,1),cumRTLL(1:end-1),'r*');hold on;
%         subplot(2,2,2)
%      plot(data2fit{1,2}.allObs(1:end-1,1),data2fit{1,2}.allObs(1:end-1,2),'ko');hold on;plot(data2fit{1,2}.allObs(1:end-1,1),cumRTHL(1:end-1),'r*');hold on;
%         subplot(2,2,3)
%      plot(data2fit{1,3}.allObs(1:end-1,1),data2fit{1,3}.allObs(1:end-1,2),'ko');hold on;plot(data2fit{1,3}.allObs(1:end-1,1),cumRTLH(1:end-1),'r*');hold on;
%         subplot(2,2,4)
%      plot(data2fit{1,4}.allObs(1:end-1,1),data2fit{1,4}.allObs(1:end-1,2),'ko');hold on;plot(data2fit{1,4}.allObs(1:end-1,1),cumRTHH(1:end-1),'r*');hold on;
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
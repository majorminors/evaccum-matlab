%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);
p = struct(); % keep some of our parameters tidy
d = struct(); % set up a structure for the data info
t = struct(); % set up a structure for temp data

% set up variables
rootdir = '\\cbsu\data\Group\Woolgar-Lab\projects\Dorian\EvAccum'; %'C:\Users\doria\Nextcloud\desiderata\desiderata\04 Research\05 Evidence Accumulation\01 EvAccum Code'; %

datadir = fullfile(rootdir,'data','behav_pilot_1-hr');
num_subjects = 14;

modeldir = fullfile(datadir,'lba_fit','results'); % expects to find your modelling results here
toolsdir = fullfile(rootdir, 'tools_analysis'); % where are all your scripts/tools?

% directory mapping
addpath(genpath(toolsdir)); % add tools folder to path
%%
% 1 - Std differ across fingers 
% 2 - C0 differ across fingers            
% 3 - B differ across conditions 
% 4 - Drift rate differ across conditions 
% 5 - t0 differ across conditions 

design_space={[1,3],[1,4],[1,3,4],[1,3,4,5],[1,2],[1,2,3],[1,2,4],[1,2,3,4],[1,2,3,4,5],[1,5],[1,3,5],[1,4,5],[1,2,5],[1,2,3,5],[1,2,4,5]};
mod_num = 2;
Model_Feature = design_space{mod_num};

[dname,IDnum] = getnames(modeldir,11);

flabs = fullfile(rootdir,'Model_%s.mat');
%%
BIC_all =[];
for m = [1:4 6:15]%1:length(design_space)%
    try
    flabs =  sprintf('Model_%s.mat',num2str(m)); 
    flabs = fullfile(modeldir,flabs);    
    load(flabs);
    catch
        continue
    end

    bestfit = []; bestfitIdx =[];iict=0;iipd=0;
    for i = 1:num_subjects
        
        [bestfit(i,1) bestfitIdx(i,1)]= min(bestval{i});
        
        BIC_all(m,i) = -0.5*BIC{i}(bestfitIdx(i,1)); %BIC needs to be negative
        
    end

end
%%
% Random Bayesian effect selection as implemented by VBA toolbox
%BIC_all(:,all(~BIC_all,1))=[];%make sure there all values different from zero

%Under H= (no difference between groups) teh group specific dataset can be
%pooled to perform a standard RFX-BMS yielding a single evidence p(y | H=)
%   - h: test decision about a difference between the group (rejection of
%        the null hypothesis of equality)
%   - pH_eq: the posterior probability that the two groups have the same model
%        frequencies



%% plot results
BIC_all(5,:)=[]; % delete removed models
[posterior_all, out_all] = VBA_groupBMC(BIC_all) ;
%out_pd;
%out_ctr;

K = size(BIC_all,1);
n = size(BIC_all,2);

% plot BICs by subject
figure; b = bar(BIC_all'*-0.5,'FaceColor',[0 0.4470 0.7410]);
% ylim([1300 1410]);
b(1).FaceColor = [.25 0 .25];
b(2).FaceColor = [.5 0 .5];
% b(5).FaceColor = [.9 0 .9];
b(9).FaceColor = [.5 0 .5];
% b(10).FaceColor = [.5 0 .5];
export_fig(fullfile(modeldir,'subjectBICs.jpeg'),'-transparent');

% plot BICs by model
figure; b = bar(BIC_all*-0.5,'FaceColor',[0 0.4470 0.7410]);
export_fig(fullfile(modeldir,'modelBICs.jpeg'),'-transparent');

%plot BIC (mean)
for i = 1:K
    meanBIC(i,1) = mean(BIC_all(i,:));
    semBIC(i,1) = nansem(BIC_all(i,:));
    meanBIC(i,1) = meanBIC(i,1)*-0.5;
end

figure; b = bar(meanBIC);
ylim([min(meanBIC)-10 max(meanBIC)+10]);
b.FaceColor = 'flat';
b.CData(1,:) = [.25 0 .25];
b.CData(2,:) = [.5 0 .5];
% b.CData(5,:) = [.9 0 .9];
b.CData(9,:) = [.5 0 .5];
% b.CData(10,:) = [.5 0 .5];
hold on
er = errorbar(1:K,meanBIC,semBIC,semBIC);
er.Color = [0 0 0];                            
er.LineStyle = 'none';
export_fig(fullfile(modeldir,'meanBICs.jpeg'),'-transparent');

fig75 = figure(75);
pos2 = [205   429   912   377];
   %set(fig75,'Render','OpenGL','Units','pixels','Color',[.95 .95 .95], 'Position', [100 100 650 525],'PaperPosition',[100 100 500 375])
   set(fig75,'Render','OpenGL','Units','pixels','Color',[.95 .95 .95], 'Position',pos2 ,'PaperPosition',[pos2(1) pos2(2) pos2(3)-150 pos2(3)-150])
   ex_Prob = out_all.ep;
   ex_Prob(ex_Prob<0.05) = 0.05;% for plotting purposes
   hs = bar(ex_Prob);
   ylim([0,1])
   set(hs,'FaceColor',[0.3843    0.4784    0.6157],'LineStyle','none')
   %set(hs,'FaceColor',[0.7    0.7    0.7],'LineStyle','none')
% 
    [~,winner] = max(out_all.ep);
    hold on;
    bar(winner,out_all.ep(winner),'r','LineStyle','none')
%    
   set(gca,'Box','off','TickDir','out','TickLength',[.03 .03],'XMinorTick','off','YMinorTick','on',...
    'YGrid','on','XColor',[.2 .2 .2],'LineWidth',1);
    set(gca,'FontName','Helvetica','Fontsize',19);
    hXLabel = xlabel('\fontsize{20} Model #');
    hYLabel = ylabel('\fontsize{20} Exceedance probability');
    xlim([0,16])
%   
     set(gcf,'Color','white')
%     
%     %export_fig(['ProbModComp','.png'],'-png','-transparent','-painters')
export_fig(fullfile(modeldir,'model_comparison.jpeg'),'-transparent');
% 
% plot frequencies
fig75 = figure(75);
clf
pos2 = [205   429   912   377];
   %set(fig75,'Render','OpenGL','Units','pixels','Color',[.95 .95 .95], 'Position', [100 100 650 525],'PaperPosition',[100 100 500 375])
   set(fig75,'Render','OpenGL','Units','pixels','Color',[.95 .95 .95], 'Position',pos2 ,'PaperPosition',[pos2(1) pos2(2) pos2(3)-150 pos2(3)-150])
   [haf,hf,hp] = plotUncertainTimeSeries(out_all.Ef,diag(out_all.Vf),[],gca);
   xlim([0,16])

   plot(gca,[0,16],[1,1]/K,'--k','LineWidth',2.5)

xlabel('models')
set(gca,'xtick',1:K,'xlim',[0.5,K+0.5],'ylim',[0 1],'ygrid','on')
if ~isempty(out_all.options.modelNames)
    set(gca,'xticklabel',out_all.options.modelNames)
end
title('estimated model frequencies')
   
set(hf,'Marker','none','LineWidth',1.5,'CapSize',0,'Color','k')   
  
set(hp,'FaceColor',[0.3843    0.4784    0.6157],'LineStyle','none')
   %set(hs,'FaceColor',[0.7    0.7    0.7],'LineStyle','none')
% 
%set(hp(2),'FaceColor','g')
    [~,winner] = max(out_all.ep);
    hold on;
    bar(winner,out_all.Ef(winner),'r','LineStyle','none')
%    
   set(gca,'Box','off','TickDir','out','TickLength',[.03 .03],'XMinorTick','off','YMinorTick','on',...
    'YGrid','on','XColor',[.2 .2 .2],'LineWidth',1);
    set(gca,'FontName','Helvetica','Fontsize',19);
   
%   
     set(gcf,'Color','white')
     h = get(gca,'Children');
     set(gca,'Children',[h(3) h(2) h(1) h(4)]);
%     %export_fig(['ProbModComp','.png'],'-png','-transparent','-painters')
export_fig(fullfile(modeldir,'model_freq.jpeg'),'-transparent');
% % plot CTR attributions
% fig76 = figure(76);
% pos2 = [ 183   257   876   638];
%    %set(fig75,'Render','OpenGL','Units','pixels','Color',[.95 .95 .95], 'Position', [100 100 650 525],'PaperPosition',[100 100 500 375])
%    set(fig76,'Render','OpenGL','Units','pixels','Color',[.95 .95 .95], 'Position',pos2 ,'PaperPosition',[-2.2767  1.4323 13.0521 8.1354])
%    
%    hi = imagesc(posterior_ctr.r');
%    axis('square')
%    xlabel('models')
% ylabel('subjects')
% title('model attributions')
% set(gca,'xlim',[0.5,K+0.5],'xtick',[1:K],'ylim',[0.5,n+0.5],'ytick',[1:n],'ydir','reverse','clim',[0 1])
% if ~isempty(out_ctr.options.modelNames)
%     set(gca,'xticklabel',out_ctr.options.modelNames)
% end
% set(gca,'visible','on')
% colormap(gca,flipud(bone))   
% blues = cbrewer('seq','Blues',100);   
% colormap(gca,(blues))   
% hc = colorbar('peer',gca,'location','NorthOutside');
% 
%    
%  set(gca,'FontName','Helvetica','Fontsize',19);
%     
% %   
%  set(gcf,'Color','white')
% %     
% %     %export_fig(['ProbModComp','.png'],'-png','-transparent','-painters')
% export_fig(fullfile(modeldir,'ctr_attributions.jpeg'),'-transparent');
% % 
% 
% % plot PD attributions
% fig76 = figure(76);
% clf
% pos2 = [ 183   257   876   638];
%    %set(fig75,'Render','OpenGL','Units','pixels','Color',[.95 .95 .95], 'Position', [100 100 650 525],'PaperPosition',[100 100 500 375])
%    set(fig76,'Render','OpenGL','Units','pixels','Color',[.95 .95 .95], 'Position',pos2 ,'PaperPosition',[-2.2767  1.4323 13.0521 8.1354])
%    
%    hi = imagesc(posterior_pd.r');
%    axis('square')
%    xlabel('models')
% ylabel('subjects')
% title('model attributions')
% set(gca,'xlim',[0.5,K+0.5],'xtick',[1:K],'ylim',[0.5,n+0.5],'ytick',[1:n],'ydir','reverse','clim',[0 1])
% if ~isempty(out_pd.options.modelNames)
%     set(gca,'xticklabel',out_pd.options.modelNames)
% end
% set(gca,'visible','on')
% %colormap(gca,flipud(bone))   
% blues = cbrewer('seq','Blues',100);   
% colormap(gca,(blues))   
% hc = colorbar('peer',gca,'location','NorthOutside');
% 
%    
%  set(gca,'FontName','Helvetica','Fontsize',19);
%     
% %   
%  set(gcf,'Color','white')
% %     
% %     %export_fig(['ProbModComp','.png'],'-png','-transparent','-painters')
% export_fig( fullfile(modeldir,'pd_attributions.jpeg'),'-transparent');
% % 
close all

%% sem function
function semval = nansem(vector_data)
% Recall that s.e.m. = std(x)/sqrt(length(x));
nonan_std = nanstd(vector_data);
nonan_len = length(vector_data(~isnan(vector_data)));
% Plug in values
semval = nonan_std / sqrt(nonan_len);
end

%% Demo for fitting LBA to one subject
% Adapted from Alessandro>Holly's 4-choice task
% DM Last Edit: Feb 2020


%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);
t = struct(); % for temp vars

% set up variables
rootdir = 'C:\Users\doria\Nextcloud\desiderata\desiderata\04 Research\05 Evidence Accumulation\01 EvAccum Code';%'\\cbsu\data\Group\Woolgar-Lab\projects\Dorian\EvAccum'; % root directory - used to inform directory mappings
datadir = fullfile(rootdir,'data\behav_pilot_2');
t.testing = 1; % if you want to use the testing data, then switch to 1 and add the data folder to the path, else to 0
t.subject = 1;

if ~t.testing
    % directory mappings
    addpath(genpath(fullfile(rootdir, 'tools_analysis'))); % add tools folder to path (includes LBA scripts)
    lbadatadir = fullfile(datadir,'lba_fit'); % find directory with prepped data
    t.fileinfo = dir(fullfile(lbadatadir,'prepped_data.mat'));
    t.datapath = fullfile(lbadatadir,t.fileinfo.name);
    
    % get the data
    t.alldata = load(t.datapath);
    t.data = t.alldata.d.subject(t.subject).data;
end

%% testing dataset
%% --------------------
if t.testing
   t.alldata = load('lba_test_data.mat');
   t.data = t.alldata.d.subject(t.subject).data;
end

%% -----------------
dataRaw = t.data; % here's the data
dataValid = []; % strip invalid responses out
for i = 1:size(dataRaw,1)
    if dataRaw{i,5} >= 0 % if there's a valid response
        dataValid(end+1,:) = [dataRaw{i,2} dataRaw{i,3} dataRaw{i,4} dataRaw{i,5}]; % add the following rows to this new variable in order: condition code, response, rt, accuracy
    end
end
data2fit={};
% get quantiles from RT data
data2fit = data_stats_cond(dataValid(:,[2 3 1])); % requires: response, rt, condition code

% fit the basic model
% startpar=[1,1,1,1,1, .1 .3]; % initial parameter [boundary,four drift rate, non-decisiontime, drift rate std]
% upper staring-point is fixed 
randiter = 100; % random search 100 iters before optimization
nosession = 1; % 20 optimization iterations

% Model fitting
Model_Feature=[1 4];
numParam=getModelParam_fixC0(Model_Feature);
parange=[zeros(1,numParam);zeros(1,numParam)+10]';
% parange(6,2)=0.4; % T0 non-decition time
[bestpar,bestval,BIC]=fitparams_refine_template('fiterror_cs_overall_FT4_template',Model_Feature,data2fit,randiter,nosession,[],parange);

% Best mdoel
[foo,best_indx]=min(bestval);
disp(['best parameter vector : ' num2str(bestpar(best_indx,:))]);
disp(['best error funcion value: ' num2str(bestval(best_indx))]);

% Re-run the best model to get simulation results

[num_par,par_Free,par_Spec]=getModelParam_fixC0(Model_Feature,bestpar(best_indx,:));
% [pmod_simpRT,qobs{1},cumRT,propInvalidTrial]=mod_stats_sim_simpRT('LBA_simpRT_template',par_simpRT,1,data2fit.allObs(:,1),1,1); %,goalstat_spec_norep.cond_ratio);
[pmod_FT,priorMod_FT,qobs_Free{1},foo3,cumRT_Free]=mod_stats_sim('LBA_spec_FT3_c_even_template',par_Free,1,data2fit.allObs(:,1),2,1); %,goalstat_free_norep.cond_ratio);


% Plotting...
figure; hold on;
plot(data2fit.allObs(1:end-1,1),data2fit.allObs(1:end-1,2),'k','Marker','o');
plot(data2fit.allObs(1:end-1,1),cumRT_Free(1:end-1),'r','Marker','*');
legend({'data','model'},'Location','SouthEast');
% title({['Model params (B, A, Astd,To): [' num2str(bestpar(best_indx,:),2) ']'] });
title('Chosen condition');
ylim([0 1]);
% xlim([0.2,0.6]);
ylabel('Cumulative probability');
xlabel('Response Time (s)');

figure; hold on;
temp=[data2fit.priorProb{1:end}];
choiceProb=[temp;priorMod_FT];
bar(choiceProb');
legend({'data','model'});
set(gca,'XTick',[1:2]);
set(gca,'XTickLabel',{'Button 1' 'Button 2'});
title('Choice probability');

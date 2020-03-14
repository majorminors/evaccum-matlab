%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);
p = struct(); % keep some of our parameters tidy
d = struct(); % set up a structure for the data info
t = struct(); % set up a structure for temp data

% set up variables
rootdir = 'C:\Users\doria\Nextcloud\desiderata\desiderata\04 Research\05 Evidence Accumulation\01 EvAccum Code'; %'\\cbsu\data\Group\Woolgar-Lab\projects\Dorian\EvAccum'; % % root directory - used to inform directory mappings

datadir = fullfile(rootdir,'data','behav_pilot_2');

modeldir = fullfile(datadir,'lba_fit','results'); % expects to find your modelling results here
toolsdir = fullfile(rootdir, 'tools_analysis'); % where are all your scripts/tools?

p.datafilepattern =  'Model_*.mat';
%p.savefilename = '';

% directory mapping
addpath(genpath(toolsdir)); % add tools folder to path (don't think we need this, but in case)
%save_file = fullfile(modeldir, p.savefilename);

%% run the data2fit analysis again

lbadatadir = fullfile(datadir,'lba_fit'); % expects to find your data here and will save results in a sub-folder here
p.data_name = 'prepped_data.mat'; % data file name

t.fileinfo = dir(fullfile(lbadatadir,p.data_name));
t.datapath = fullfile(lbadatadir,t.fileinfo.name);

% get the data
t.alldata = load(t.datapath);
t.data = t.alldata.d.subjects;

num_subjs = sum(~cellfun('isempty',{t.data.id})); % get the number of not empty arrays from the first field of the 'data' structure
%initialise cells
d.data2fit={};

for idxsubj = 1:num_subjs
    fprintf('loading subject %1.0f \n',t.data(idxsubj).id);
    
    % pull the data
    subjdata = t.data(idxsubj).data; % put the data in a readable variable
    dataValid = []; % now we'll strip invalid responses out
    validLabs = {};
    for i = 1:size(subjdata,1)
        if subjdata{i,5} >= 0 % if there's a valid response
            dataValid(end+1,:) = [subjdata{i,2} subjdata{i,3} subjdata{i,4} subjdata{i,5}]; % add the following rows to this new variable in order: condition code, response, rt, accuracy
            temp = subjdata{i,1}; % extract valid condition label (can't do in one step with multi-level structures and non-scalar indexing)
            validLabs{end+1} = temp; % add the valid condition label
        end
    end; clear i temp subjdata;
    
    % converting again to readable variables
    conds = dataValid(:,1);
    resp = dataValid(:,2);
    rt = dataValid(:,3);
    acc = dataValid(:,4);
    
    % do descriptive stats
    minRT(idxsubj) = min(rt(rt>=0.1)); % not used
    
    
    %%
    %Pool conditions according to design matrix & calculate stats
    unique_conds = unique(conds);
    for level = 1:length(unique_conds) % for all conditions
        
        trialn = conds == level;
        dataRaw{idxsubj,level}  = [resp(trialn) rt(trialn)];
        
        d.data2fit{idxsubj,level} = data_stats(dataRaw{idxsubj,level});
        
        % add some info to d.data2fit (the condition labels as a string, and
        % the minRT - in the case that the min rt is smaller than the
        % non-decision time - let's leave this for now, just in case)
        d.data2fit{idxsubj,level}.cond = validLabs{trialn};
        d.data2fit{idxsubj,level}.minRT = round(minRT(idxsubj),2);
        %parange(end,1) = round(minRT(idxsubj),2);%constrain the lower bound of T0 to the shortest RT
    end; clear level trialn
    
    clear rt conds resp acc dataRaw dataValid validLabs minRT % clear vars from behavioural data to avoid reading mixed subject data
end; clear idxsubj;

%% get model data

d.fileinfo = dir(fullfile(modeldir, p.datafilepattern)); % find all the datafiles and get their info
for i = 1%:length(d.fileinfo) % loop through each
    t.path = fullfile(modeldir, d.fileinfo(i).name);% get the full path to the file
    fprintf(1, 'working with %s\n', t.path); % print that so you can check
    
    t.model_results = load(t.path); % load in the data
    
    %% plotting
    
    for idxsubj = 1:num_subjs
        [~,t.id] = min(t.model_results.bestval{idxsubj});
        t.bestparams = t.model_results.bestpar{idxsubj}(t.id,:);
        
        for level = 1:length(unique_conds)
            t.data2fit{level} = d.data2fit{idxsubj,level};
        end
        
        [~,parLL,parHL,parLH,parHH]=getModelParam_cell_RDK(t.model_results.settings.modfeat,2,t.bestparams);
        t.params(1) = parLL;
        t.params(2) = parLH;
        t.params(3) = parHL;
        t.params(4) = parHH;
        clear parLL parLH parHL parHH
        
        t.probmod = []; t.cumRT = [];
        for level = 1:length(unique_conds) % for all conditions
            [~,t.probmod(:,level),~,~,t.cumRT(:,level)]=mod_stats_sim('LBA_spec_FT3_c_even_template',t.params(level),1,t.data2fit{:,level}.allObs(:,1),2,1);
        end; clear level;
        
        % Plotting...
        if length(unique_conds) ~= length(t.data2fit); error('you dont appear to have as many conditions as sets of t.data2fit'); end % sanity check
        
        % hold onto T0 to plot later
        non_dec_time(idxsubj) = t.params.T0;
        
%         % plot parameters
%         %   t.params(1) = LL; t.params(2) = LH; t.params(3) = HL; t.params(4) = HH;
%         %   param names = B C0 Ame Astd T0
%         figure; hold on;
%         temp=[mean(t.params(1).B),mean(t.params(2).B),mean(t.params(3).B),mean(t.params(4).B)];
%         bar(temp');
%         ylim([min(temp)-1 max(temp)+1]);
%         %legend({'Decision Boundary'});
%         set(gca,'XTick',[1:4]);
%         set(gca,'XTickLabel',{'LL' 'LH' 'HL' 'HH'});
%         title('Decision Boundary');
        
        for level = 1:length(unique_conds) % for all conditions

%             % plot quintiles
%             figure; hold on;
%             plot(t.data2fit{level}.allObs(1:end-1,1),t.data2fit{level}.allObs(1:end-1,2),'k','Marker','o');
%             plot(t.data2fit{level}.allObs(1:end-1,1),t.cumRT(1:end-1,level),'r','Marker','*');
%             legend({'data','model'},'Location','SouthEast');
%             % title({['Model params (B, A, Astd,To): [' num2str(t.model_results.bestpar(t.id,:),2) ']'] });
%             title('Chosen condition');
%             ylim([0 1]);
%             % xlim([0.2,0.6]);
%             ylabel('Cumulative probability');
%             xlabel('Response Time (s)');
%             
%             % plot choices
%             figure; hold on;
%             temp=[t.data2fit{level}.priorProb{1,1}(1),t.data2fit{level}.priorProb{1,2}(1)];
%             t.choiceProb=[temp;t.probmod(:,level)'];
%             bar(t.choiceProb');
%             legend({'data','model'});
%             set(gca,'XTick',[1:2]);
%             set(gca,'XTickLabel',{'Button 1' 'Button 2'});
%             title('Choice probability');
            
        end; clear level;
    end; clear idxsubj;
    
    % plot T0
    figure; hold on;
    bar(non_dec_time);
    ylim([min(non_dec_time)-0.5 max(non_dec_time)+0.5]);
    %legend({'Decision Boundary'});
    set(gca,'XTick',[1:length(non_dec_time)]);
    %set(gca,'XTickLabel',{'LL' 'LH' 'HL' 'HH'});
    title('Non-Decision Time');
    
end; clear i;


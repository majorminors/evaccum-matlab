%% Get some descriptives
% Dorian Minors
% Created: MAR20
% Last Edit: APR20
%
% 
% you specify your datadir and it will save output data into
% datadir\descriptives (can change in 'outdir') as whatever you specify as p.savefilename
% produces:
%
% d.subjects struct saved as .mat file
% d.allsubjs matrix saved as .csv file
% notes file describing output in .txt
% 
% d.fileinfo = information about matlab files read in
% d.subjects = rowise subject data
%   fields:
%       1) subject id
%       2) all data organised into a cell matrix thus: condition* | condition number* | button press | reaction time (in ms) | accuracy | trial type** | 
%       3) data stripped of condition labels and invalid responses organised into a normal matrix thus: condition number* | button press | reaction time (in ms) | accuracy | trial type** | 
%       4) mean according to condition* (1-4)
%       5) sem according to condition* (1-4)
%       6) percent correct according to condition* (1-4)
%       7) all descriptives collated as [mean;sem;accuracy] for easy access
%
%       *conditions are:
%               LcLr = 1 = low coherence, easy matching (low rule)
%               LcHr = 2 = low coherence, hard matching (hard rule)
%               HcLr = 3 = high coherence, easy matching (low rule)
%               HcHr = 4 = high coherence, hard matching (hard rule)
%       **trial type (1-64) is each unique trial condition - 2
%               coherence levels x 2 matching difficulties x 8 coherence
%               directions x 2 button presses
%
% d.allsubjs = cell matrix of rowise subject data, but with all subjects collated in one matrix, and
%              subject id appended as first column
% d.allsubjs = the same, but in a normal matrix and stripped of condition labels and invalid
%              responses
% d.mean = overall mean according to condition* (1-4)
% d.sem = overall sem according to condition* (1-4)
% d.accuracy = overall percent correct according to condition* (1-4)
% d.alldescriptives = all overall descriptives collated as [mean;sem;accuracy] for easy access

% saves all into d into whatever you call p.savefilename, and saves a figure for the mean, and a figure for the
% overall accuracy

%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);
p = struct(); % keep some of our parameters tidy - these are things you might want to edit
d = struct(); % set up a structure for the data info - these are outputs we want to save
t = struct(); % set up a structure for temp data - these are things that are liable to change throughout the script

% set up variables
rootdir = '\\cbsu\data\Group\Woolgar-Lab\projects\Dorian\EvAccum'; % % root directory - used to inform directory mappings
datadir = fullfile(rootdir,'data','behav_pilot_1');
p.strip_invalid = 1; % will strip invalid responses from output if 1 from allsubjs (can't get individual subjs working)
p.datafilepattern = '*_EvAccum.mat';
p.vars = {'p','d','block'}; % which variables do you need
p.savefilename = 'descriptives';
p.notesfilename = [p.savefilename,'_notes.txt'];
p.notes = 'see top of descriptives.m for organisation of data';
p.conditions = {'LcLr','LcHr','HcLr','HcHr'}; % 2x2 coherence and rule
p.conditioncodes = {1,2,3,4};

% directory mapping
addpath(genpath(fullfile(rootdir, 'tools_analysis'))); % add tools folder to path (don't think we need this, but in case)
behavdatadir = fullfile(datadir,'behavioural'); % add matlab behavioural data
outdir = fullfile(datadir,'descriptives'); % find or make directory to output lba fit results
if ~exist(outdir,'dir')
    mkdir(outdir);
end
save_file = fullfile(outdir,p.savefilename);

%% get that data

d.fileinfo = dir(fullfile(behavdatadir, p.datafilepattern)); % find all the datafiles and get their info
d.allsubjs = {}; d.allsubjsstripped = []; % create these so we can work with it
for i = 1:length(d.fileinfo) % loop through each
  t.path = fullfile(behavdatadir, d.fileinfo(i).name); % get the full path to the file
  fprintf(1, 'working with %s\n', t.path); % print that so you can check
  
  t.alldata = load(t.path,p.vars{:}); % load in the data
  
  t.id = t.alldata.d.participant_id;
  t.blocks = t.alldata.block; % how many blocks?
  
  t.data = {}; t.strippeddata = []; % create these so we can work with it
  for block = 1:t.blocks % go through each block of data
      t.rts = t.alldata.d.rt(block,:)';
      t.accuracy = t.alldata.d.correct(block,:)';
      % this is the wrong button - correct button t.button = t.alldata.d.stim_mat_all(:,7,block); % pull the correct button
      % t.trialtype = t.alldata.d.stim_mat_all(:,9,block);
      
      for icond = 1:length(t.alldata.d.stim_mat_all(:,5,block))
          % get the button pressed
          if strcmp(t.alldata.p.resp_keys{1},t.alldata.d.resp_key_name{block,icond})
              t.button(icond,1) = 1;
          elseif strcmp(t.alldata.p.resp_keys{2},t.alldata.d.resp_key_name{block,icond})
              t.button(icond,1) = 2;
          else
              t.button(icond,1) = 0;
          end
          
          % create a row that gives you a number for each condition in your 2x2
          if t.alldata.d.stim_mat_all(icond,5)==1 && t.alldata.d.stim_mat_all(icond,8)==1
              t.condition(icond,1) = string(p.conditions{1});
              t.conditioncode(icond,1) = p.conditioncodes{1};
          elseif t.alldata.d.stim_mat_all(icond,5)==1 && t.alldata.d.stim_mat_all(icond,8)==2
              t.condition(icond,1) = string(p.conditions{2});
              t.conditioncode(icond,1) = p.conditioncodes{2};
          elseif t.alldata.d.stim_mat_all(icond,5)==2 && t.alldata.d.stim_mat_all(icond,8)==1
              t.condition(icond,1) = string(p.conditions{3});
              t.conditioncode(icond,1) = p.conditioncodes{3};
          elseif t.alldata.d.stim_mat_all(icond,5)==2 && t.alldata.d.stim_mat_all(icond,8)==2
              t.condition(icond,1) = string(p.conditions{4});
              t.conditioncode(icond,1) = p.conditioncodes{4};
          end
      end
      clear icond

      % consolidate all that data
      t.consolidata(:,1) = num2cell(t.condition);
      t.consolidata(:,2) = num2cell(t.conditioncode);
      t.consolidata(:,3) = num2cell(t.button);
      t.consolidata(:,4) = num2cell(t.rts);
      t.consolidata(:,5) = num2cell(t.accuracy);
      % t.consolidata(:,6) = num2cell(t.trialtype);
      
      % stack it up
      t.data = [t.data;t.consolidata];
      t.strippeddata = [t.strippeddata;t.consolidata(:,2:end)];
  end
  clear block
  
  % strip invalid
  t.strippeddata = cell2mat(t.strippeddata);
  t.invalid = t.strippeddata(:,4) == -1;
  t.strippeddata(t.invalid,:) = [];
  
  % put it all in a persistent variable
  d.subjects(i).id = t.id;
  d.subjects(i).data = t.data;
  d.subjects(i).strippeddata = t.strippeddata;
  
  t.idcolumn = ones(length(d.subjects(i).data),1)*t.id;
  t.idcolstripped = ones(length(d.subjects(i).strippeddata),1)*t.id;
  
  t.thissubj = [num2cell(t.idcolumn) d.subjects(i).data];
  t.thissubjstripped = [t.idcolstripped d.subjects(i).strippeddata];
  
  d.allsubjs = [d.allsubjs;t.thissubj];
  d.allsubjsstripped = [d.allsubjsstripped;t.thissubjstripped];
end; clear i

%% get some descriptives
for subj = 1:length(d.fileinfo)
    for condidx = 1:4
        t.condition = find(d.subjects(subj).strippeddata(:,1) == condidx);
        d.subjects(subj).mean(condidx) = nanmean(d.subjects(subj).strippeddata(t.condition,3));
        d.subjects(subj).sem(condidx) = nansem(d.subjects(subj).strippeddata(t.condition,3));
        t.accuracy = d.subjects(subj).strippeddata(t.condition,4);
        d.subjects(subj).accuracy(condidx) = sum(t.accuracy)/length(t.accuracy)*100;
    end; clear condidx t
    d.subjects(subj).alldescriptives = [d.subjects(subj).mean;d.subjects(subj).sem;d.subjects(subj).accuracy];
end; clear subj

for condidx = 1:4
    t.condition = find(d.allsubjsstripped(:,2) == condidx);
    d.mean(condidx) = nanmean(d.allsubjsstripped(t.condition,4));
    d.sem(condidx) = nansem(d.allsubjsstripped(t.condition,4));
    t.accuracy = d.allsubjsstripped(t.condition,5);
    d.accuracy(condidx) = sum(t.accuracy)/length(t.accuracy)*100;
end; clear condidx t
d.alldescriptives = [d.mean;d.sem;d.accuracy];

%% figures
figure; b = bar(d.mean);
ylim([min(d.mean)-0.01 max(d.mean)+0.01]);
b.FaceColor = 'flat';
hold on
er = errorbar(1:size(d.mean,2),d.mean,d.sem);
er.Color = [0 0 0];                            
er.LineStyle = 'none';
export_fig(fullfile([save_file 'mean_fig.jpeg']),'-transparent');

figure; b = bar(d.accuracy);
ylim([min(d.accuracy)-1 max(d.accuracy)+1]);
b.FaceColor = 'flat';
export_fig(fullfile([save_file 'accuracy_fig.jpeg']),'-transparent');

%% saving data and notes
fprintf('saving output and notes from %s\n', mfilename);
txtfile = fopen(fullfile(outdir,p.notesfilename),'w');
fprintf(txtfile,p.notes);
fclose(txtfile); clear txtfile;
save(save_file,'-struct','d'); % save all data to a .mat file
% writecell(d.allsubjs,[save_file '_alldata.csv'])
% writematrix(d.allsubjsstripped,[save_file '_strippeddata.csv']);
% writematrix(d.alldescriptives,[save_file '_descriptivesonly.csv']);

%% sem function
function semval = nansem(vector_data)
% Recall that s.e.m. = std(x)/sqrt(length(x));
nonan_std = nanstd(vector_data);
nonan_len = length(vector_data(~isnan(vector_data)));
% Plug in values
semval = nonan_std / sqrt(nonan_len);
end
  
%% prep data for input to LBA fitting
% Dorian Minors
% Created: FEB20
% Last Edit: FEB20
%
% 
% you specify your datadir and it will save prepped data into datadir\lba_fit
% as whatever you specify as p.savefilename
%
% produces:
% 
% d.fileinfo = information about matlab files read in
% d.subjects = rowise subject data, with subject id in first field and
%             lba-relevent data in second field organised thus:
%             condition | condition number | button press | reaction time (in ms) | accuracy | trial type | 
%
%             trial type (1-64) is each unique trial condition - 2
%               coherence levels x 2 matching difficulties x 8 coherence
%               directions x 2 button presses
%             conditions are:
%               LcLr = 1 = low coherence, easy matching (low rule)
%               LcHr = 2 = low coherence, hard matching (hard rule)
%               HcLr = 3 = high coherence, easy matching (low rule)
%               HcHr = 4 = high coherence, hard matching (hard rule)

%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);
p = struct(); % keep some of our parameters tidy
d = struct(); % set up a structure for the data info
t = struct(); % set up a structure for temp data

% set up variables
rootdir = '\\cbsu\data\Group\Woolgar-Lab\projects\Dorian\EvAccum'; %'C:\Users\doria\Nextcloud\desiderata\desiderata\04 Research\05 Evidence Accumulation\01 EvAccum Code'; % root directory - used to inform directory mappings
datadir = fullfile(rootdir,'data\behav_pilot_2');
p.datafilepattern = '*_EvAccum.mat';
p.savefilename = 'prepped_data';
p.notesfilename = [p.savefilename,'_notes.txt'];
p.notes = 'd.fileinfo = information about matlab files read in \r\nd.subjects = rowise subject data, with subject id in first field and lba-relevent data in second field organised thus: \r\ncondition | condition number | button pressed | reaction time (ms) | accuracy | trial type (1-64) | condition (LcLr,LcHr,HcLr,HcHr)';
t.conditions = {'LcLr','LcHr','HcLr','HcHr'}; % 2x2 coherence and rule
t.conditioncodes = {1,2,3,4};

% directory mapping
addpath(genpath(fullfile(rootdir, 'tools_analysis'))); % add tools folder to path (don't think we need this, but in case)
behavdatadir = fullfile(datadir,'behavioural'); % add matlab behavioural data
lbadatadir = fullfile(datadir,'lba_fit'); % find or make directory to output lba fit results
if ~exist(lbadatadir,'dir')
    mkdir(lbadatadir);
end
save_file = fullfile(lbadatadir, p.savefilename);

%% get that data

d.fileinfo = dir(fullfile(behavdatadir, p.datafilepattern)); % find all the datafiles and get their info
for i = 1:length(d.fileinfo) % loop through each
  t.path = fullfile(behavdatadir, d.fileinfo(i).name); % get the full path to the file
  fprintf(1, 'working with %s\n', t.path); % print that so you can check
  
  t.alldata = load(t.path); % load in the data
  
  t.id = t.alldata.d.participant_id;
  t.blocks = t.alldata.block; % how many blocks?
  
  t.data = {}; % create this so we can work with it
  for block = 1:t.blocks % go through each block of data
      t.rts = t.alldata.d.rt(block,:)';
      t.accuracy = t.alldata.d.correct(block,:)';
      % this is the wrong button - correct button t.button = t.alldata.d.stim_mat_all(:,7,block); % pull the correct button
      t.trialtype = t.alldata.d.stim_mat_all(:,9,block);
      
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
              t.condition(icond,1) = string(t.conditions{1});
              t.conditioncode(icond,1) = t.conditioncodes{1};
          elseif t.alldata.d.stim_mat_all(icond,5)==1 && t.alldata.d.stim_mat_all(icond,8)==2
              t.condition(icond,1) = string(t.conditions{2});
              t.conditioncode(icond,1) = t.conditioncodes{2};
          elseif t.alldata.d.stim_mat_all(icond,5)==2 && t.alldata.d.stim_mat_all(icond,8)==1
              t.condition(icond,1) = string(t.conditions{3});
              t.conditioncode(icond,1) = t.conditioncodes{3};
          elseif t.alldata.d.stim_mat_all(icond,5)==2 && t.alldata.d.stim_mat_all(icond,8)==2
              t.condition(icond,1) = string(t.conditions{4});
              t.conditioncode(icond,1) = t.conditioncodes{4};
          end
      end
      clear icond

      % consolidate all that data
      t.consolidata(:,1) = num2cell(t.condition);
      t.consolidata(:,2) = num2cell(t.conditioncode);
      t.consolidata(:,3) = num2cell(t.button);
      t.consolidata(:,4) = num2cell(t.rts);
      t.consolidata(:,5) = num2cell(t.accuracy);
      t.consolidata(:,6) = num2cell(t.trialtype);
      
      % stack it up
      t.data = [t.data;t.consolidata];  
  end
  clear block
  
  d.subjects(i).id = t.id;
  d.subjects(i).data = t.data;

end
clear i

%% saving data and notes
fprintf('saving output and notes from %s\n', mfilename);
txtfile = fopen(fullfile(lbadatadir,p.notesfilename),'w');
fprintf(txtfile,p.notes);
fclose(txtfile); clear txtfile;
save(save_file,'d'); % save all data to a .mat file

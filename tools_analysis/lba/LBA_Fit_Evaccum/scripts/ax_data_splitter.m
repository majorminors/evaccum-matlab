%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);
p = struct(); % for parameters
t = struct(); % for temp vars

%% set up variables

% required
rootdir = 'C:\Users\doria\Nextcloud\desiderata\desiderata\04 Research\05 Evidence Accumulation\01 EvAccum Code';%'\\cbsu\data\Group\Woolgar-Lab\projects\Dorian\EvAccum'; % root directory - used to inform directory mappings
datadir = fullfile(rootdir,'data','behav_pilot_1','lba_fit'); % expects to find your data here and will save results in a sub-folder here
p.data_name = 'prepped_data.mat'; % data file name

%% get the data
t.fileinfo = dir(fullfile(datadir,p.data_name));
t.datapath = fullfile(datadir,t.fileinfo.name);
t.alldata = load(t.datapath);

t.subjects = t.alldata.d.subjects; % here is your data

%% split the data

for subj = 1:size(t.subjects,2)
    t.data = t.subjects(subj).data;
    
    
    
    t.easy_data = []; t.hard_data = []; t.easy_labs = {}; t.hard_labs = {};
    for i = 1:size(t.data,1)
        if t.data{i,2} == 1 || t.data{i,2} == 3; % if there's a valid response
            t.easy_data(end+1,:) = [t.data{i,2} t.data{i,3} t.data{i,4} t.data{i,5}];% t.data{i,6}]; % add the following rows to this new variable in order: condition code, response, rt, accuracy
            temp = t.data{i,1}; % extract valid condition label (can't do in one step with multi-level structures and non-scalar indexing)
            t.easy_labs{end+1} = temp; % add the valid condition label
        elseif t.data{i,2} == 2 || t.data{i,2} == 4;
            t.hard_data(end+1,:) = [t.data{i,2} t.data{i,3} t.data{i,4} t.data{i,5}];% t.data{i,6}]; % add the following rows to this new variable in order: condition code, response, rt, accuracy
            temp = t.data{i,1}; % extract valid condition label (can't do in one step with multi-level structures and non-scalar indexing)
            t.hard_labs{end+1} = temp; % add the valid condition label
        end; clear temp;
    end; clear i;
    
    %% these are the ones we've messed with
    t.subjects(subj).easy_data = [t.easy_labs' num2cell(t.easy_data)];
    t.subjects(subj).hard_data = [t.hard_labs' num2cell(t.hard_data)];
    
    %% outputs
    d.subjects(subj).id = t.subjects(subj).id;
    d.subjects(subj).data = t.subjects(subj).easy_data;

end; clear subj;

p.save_name = 'easy_data.mat';
p.save_file = fullfile(datadir, p.save_name);


save(p.save_file,'d'); % save all data to a .mat file

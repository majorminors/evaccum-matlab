%% stitches two matlab datafiles together

%% set up

close all;
clearvars;
clc;

% enter subject
subject = 'S02';

% directory mapping
stitcherrootdir = '\\cbsu\data\Group\Woolgar-Lab\projects\Dorian\EvAccum\data\meg_pilot_1\behavioural\S02'; % named thus because we have a 'rootdir' in the matlab files we're loading in

% enter filename
data_file_1 = strcat(subject,'_1_EvAccum');
data_file_2 = strcat(subject,'_2_EvAccum');

% how many blocks in each file
numblocks1 = 6;
numblocks2 = 4;

% what do you want to keep it as
stitcher_save_file = strcat(subject, '_EvAccum_stitched'); % named thus because we have a 'save_file' in the matlab files we're loading in

% filenames
ds1 = fullfile(stitcherrootdir,[data_file_1 '.mat']);
ds2 = fullfile(stitcherrootdir,[data_file_2 '.mat']);

% load first dataset
load(ds1);

% pull relevant data from first dataset
for i = 1:numblocks1
    new_stim_mat(:,:,i) = d.stim_mat_all(:,:,i);
    new_rt(i,:) = d.rt(i,:);
    new_correct(i,:) = d.correct(i,:);
end

clearvars -except stitcherrootdir numblocks1 numblocks2 ds2 stitcher_save_file new_stim_mat new_rt new_correct

% load second dataset
load(ds2);

% pull relevant data from second dataset
for i = 1:numblocks2
    new_stim_mat(:,:,i+numblocks1) = d.stim_mat_all(:,:,i);
    new_rt(i+numblocks1,:) = d.rt(i,:);
    new_correct(i+numblocks1,:) = d.correct(i,:);
end

clearvars -except stitcherrootdir stitcher_save_file new_stim_mat new_rt new_correct

d.stim_mat_all = new_stim_mat;
d.rt = new_rt;
d.correct = new_correct;

save(fullfile(stitcherrootdir, stitcher_save_file));
%% stitches multiple matlab datafiles together

%% set up

close all;
clearvars;
clc;

% enter filename
data_file_1 = 'S20_EvAccum';
data_file_2 = 'S22_EvAccum';

% what do you want to keep it as
stitcher_save_file = 'S20_EvAccum_stitched'; % named thus because we have a 'save_file' in the matlab files we're loading in

% directory mapping
stitcherrootdir = 'C:\Users\doria\Desktop\to stitch'; % named thus because we have a 'rootdir' in the matlab files we're loading in

% filenames
ds1 = fullfile(stitcherrootdir,[data_file_1 '.mat']);
ds2 = fullfile(stitcherrootdir,[data_file_2 '.mat']);

% load first dataset
load(ds1);

% pull relevant data from first dataset
for i = 1:5
    new_stim_mat(:,:,i) = d.stim_mat_all(:,:,i);
    new_rt(i,:) = d.rt(i,:);
    new_correct(i,:) = d.correct(i,:);
end

clearvars -except ds1 ds2 stitcherrootdir stitcher_save_file new_stim_mat new_rt new_correct

% load second dataset
load(ds2);

% pull relevant data from second dataset
for i = 1:5
    new_stim_mat(:,:,i+5) = d.stim_mat_all(:,:,i);
    new_rt(i+5,:) = d.rt(i,:);
    new_correct(i+5,:) = d.correct(i,:);
end

clearvars -except ds1 ds2 stitcherrootdir stitcher_save_file new_stim_mat new_rt new_correct

% load first dataset
load(ds1);

% pull relevant data from first dataset
for i = 6:10
    new_stim_mat(:,:,i+5) = d.stim_mat_all(:,:,i);
    new_rt(i+5,:) = d.rt(i,:);
    new_correct(i+5,:) = d.correct(i,:);
end

clearvars -except ds1 ds2 stitcherrootdir stitcher_save_file new_stim_mat new_rt new_correct

% load second dataset
load(ds2);

% pull relevant data from second dataset
for i = 6:10
    new_stim_mat(:,:,i+10) = d.stim_mat_all(:,:,i);
    new_rt(i+10,:) = d.rt(i,:);
    new_correct(i+10,:) = d.correct(i,:);
end

clearvars -except ds1 ds2 stitcherrootdir stitcher_save_file new_stim_mat new_rt new_correct

d.stim_mat_all = new_stim_mat;
d.rt = new_rt;
d.correct = new_correct;

save(fullfile(stitcherrootdir, stitcher_save_file));
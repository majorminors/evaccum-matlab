%% set up

close all;
clearvars;
clc;

% set up variables
data_file = 'S20_EvAccum_stitched'; % enter filename
vars = {'p','d'}; % which variables do you need
rootdir = 'C:\Users\doria\Google Drive\04 Research\05 Evidence Accumulation\01 EvAccum Code';
num_blocks = 20;

% directory mapping
datadir = fullfile(rootdir, 'data');
backupdir = fullfile(datadir, 'backup');
if ~exist(backupdir,'dir')
    mkdir(backupdir);
end
analysisdir = fullfile(datadir, 'analysis');
if ~exist(analysisdir,'dir')
    mkdir(analysisdir);
end
subdir = fullfile(analysisdir, data_file);
if ~exist(fullfile(analysisdir, data_file),'dir')
    mkdir(fullfile(analysisdir, data_file));
end




% load file
load(fullfile(datadir,[data_file '.mat']),vars{:});
if exist(fullfile(backupdir, [data_file '_backup.mat']),'file') % check if the backup already exists and throw a warning if it does
    warning('the following backup file already exists - overwrite? (y/n)\n %s.mat', data_file);
    while 1 % loop forever until y or n
        ListenChar(2);
        [secs,keyCode] = KbWait; % wait for response
        key_name = KbName(keyCode); % find out name of key that was pressed
        if strcmp(key_name, 'y')
            fprintf('instructed to overwrite: overwriting and continuing with %s\n', mfilename)
            ListenChar(0);
            clear secs keyCode key_name
            break % break the loop and continue
        elseif strcmp(key_name, 'n')
            ListenChar(0);
            clear secs keyCode key_name
            error('instructed not to overwrite: aborting %s\n', mfilename); % error out
        end
    end % end response loop
end % end check backup file exist
save(fullfile(backupdir, [data_file '_backup'])); % backup subject data before we mess with it (will try to overwrite every time you run the script on the same data)

% reminder of stimulus condition columns in p.stim_mat
%  1)  cue direction (1-4)
%  2)  cue direction in degrees
%  3)  dot motion direction condition (1-8)
%  4)  dot motion direction in degrees
%  5)  coherence difficulty (1 = easy, 2 = hard)
%  6)  matching distance from cue direction (absolute value) - used to calc match and match difficulty
%  7)  match (1) or non-match (2)
%  8)  matching difficulty (1 = easy, 2 = difficult)

%% average across trial conditions

% make them into a single matrix (see below for layout)
stim_mat = d.stim_mat_all;
[data_summary_all,sem_summary_all] = get_rts(d,stim_mat,num_blocks);
clear stim_mat;
for i = 1:num_blocks
    temp_block = d.stim_mat_all(:,7,i)==1;
    stim_mat(:,:,i) = d.stim_mat_all(temp_block,:,i);
end
[data_summary_match,sem_summary_match] = get_rts(d,stim_mat,num_blocks);
clear stim_mat temp_block i;
for i = 1:num_blocks
    temp_block = d.stim_mat_all(:,7,i)==2;
    stim_mat(:,:,i) = d.stim_mat_all(temp_block,:,i);
end
[data_summary_nonmatch,sem_summary_nonmatch] = get_rts(d,stim_mat,num_blocks);
clear stim_mat temp_block i;

data_summary = [data_summary_all,data_summary_match,data_summary_nonmatch];

sem_summary = [nanmean(sem_summary_all.easycohcorr),nanmean(sem_summary_all.hardcohcorr),nanmean(sem_summary_all.easyrulecorr),nanmean(sem_summary_all.hardrulecorr);...
    nanmean(sem_summary_all.easycoh_easyrulecorr),nanmean(sem_summary_all.easycoh_hardrulecorr),nanmean(sem_summary_all.hardcoh_easyrulecorr),nanmean(sem_summary_all.hardcoh_hardrulecorr)];

fig1 = figure(1);
lab = categorical({'easy coh', 'hard coh', 'easy rule', 'hard rule'});
lab = reordercats(lab,{'easy coh', 'hard coh', 'easy rule', 'hard rule'});
bar = bar(lab,data_summary_all(3,:));
bar.FaceColor = 'flat';
bar.CData(1,:) = [.5 0 .5];
bar.CData(2,:) = [.5 0 .10];
bar.CData(3,:) = [0 0 255];
bar.CData(4,:) = [.10 0 .5];
hold on
er = errorbar(lab,data_summary_all(3,:),sem_summary(1,:));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off

fig2 = figure(2);
%lab = categorical({'easy coh', 'hard coh', 'easy rule', 'hard rule'});
%lab = reordercats(lab,{'easy coh', 'hard coh', 'easy rule', 'hard rule'});
bar = bar(lab,data_summary_all(6,:));
bar.FaceColor = 'flat';
bar.CData(1,:) = [.5 0 .5];
bar.CData(2,:) = [.5 0 .10];
bar.CData(3,:) = [0 0 255];
bar.CData(4,:) = [.10 0 .5];
hold on
er = errorbar(lab,data_summary_all(6,:),sem_summary(2,:));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off

% you'll end up with a table like this (without the labels):
%                 | easy coh |  hard coh | easy rule | hard rule |
% average rt      |          |           |           |           |
% percent correct |          |           |           |           |
% average correct |          |           |           |           |
%                      easy coherence        hard coherence
%                 |easy rule | hard rule | easy rule | hard rule |
% average rt      |          |           |           |           |
% percent correct |          |           |           |           |
% average correct |          |           |           |           |

% repeated 3x horizontally: | all results | results for blue | results for orange |

% save all that
dlmwrite(fullfile(subdir, [data_file '_datasummary.csv']), data_summary);
save(fullfile(subdir, [data_file '_datasummary']), 'data_summary');

%% averaging function
function [data_summary,sem_summary] = get_rts(d,stim_mat,num_blocks)

for i = 1:num_blocks % loop around the blocks
    
    % get all the rts for coherence and rule conditions separately
    rt_easy_coh(i,:) = d.rt(i,find(stim_mat(:,5,i)==1));
    rt_hard_coh(i,:) = d.rt(i,find(stim_mat(:,5,i)==2));
    rt_easy_rule(i,:) = d.rt(i,find(stim_mat(:,8,i)==1));
    rt_hard_rule(i,:) = d.rt(i,find(stim_mat(:,8,i)==2));
    
    % rts for the combinations
    rt_easy_coh_easy_rule(i,:) = d.rt(i,find(stim_mat(:,5,i)==1 & stim_mat(:,8,i)==1));
    rt_easy_coh_hard_rule(i,:) = d.rt(i,find(stim_mat(:,5,i)==1 & stim_mat(:,8,i)==2));
    rt_hard_coh_easy_rule(i,:) = d.rt(i,find(stim_mat(:,5,i)==2 & stim_mat(:,8,i)==1));
    rt_hard_coh_hard_rule(i,:) = d.rt(i,find(stim_mat(:,5,i)==2 & stim_mat(:,8,i)==2));

    
    % average all those
    av_rt_easy_coh(i,:) = nanmean(rt_easy_coh(i,:));
    sem_rt_easy_coh(i,:) = nansem(rt_easy_coh(i,:));
    av_rt_hard_coh(i,:) = nanmean(rt_hard_coh(i,:));
    sem_rt_hard_coh(i,:) = nansem(rt_hard_coh(i,:));
    av_rt_easy_rule(i,:) = nanmean(rt_easy_rule(i,:));
    sem_rt_easy_rule(i,:) = nansem(rt_easy_rule(i,:));
    av_rt_hard_rule(i,:) = nanmean(rt_hard_rule(i,:));
    sem_rt_hard_rule(i,:) = nansem(rt_hard_rule(i,:));
    av_rt_easy_coh_easy_rule(i,:) = nanmean(rt_easy_coh_easy_rule(i,:));
    sem_rt_easy_coh_easy_rule(i,:) = nansem(rt_easy_coh_easy_rule(i,:));
    av_rt_easy_coh_hard_rule(i,:) = nanmean(rt_easy_coh_hard_rule(i,:));
    sem_rt_easy_coh_hard_rule(i,:) = nansem(rt_easy_coh_hard_rule(i,:));
    av_rt_hard_coh_easy_rule(i,:) = nanmean(rt_hard_coh_easy_rule(i,:));
    sem_rt_hard_coh_easy_rule(i,:) = nansem(rt_hard_coh_easy_rule(i,:));
    av_rt_hard_coh_hard_rule(i,:) = nanmean(rt_hard_coh_hard_rule(i,:));
    sem_rt_hard_coh_hard_rule(i,:) = nansem(rt_hard_coh_hard_rule(i,:));
    
    % find the percent correct and rts of correct trials for conditions seperately    
    temp = d.correct(i,find(stim_mat(:,5,i)==1)); % d.correct for easy coherence trials
    rt_easy_coh_correcttrials = rt_easy_coh(find(temp==1)); % rt for correct trials
    correct = sum(temp == 1);
    incorrect = sum(temp == 0);
    total = correct+incorrect;
    pc_easy_coh(i,:) = correct/total*100;
    clear temp correct incorrect total;
    
    temp = d.correct(1,find(stim_mat(:,5,i)==2)); % d.correct for hard coherence trials
    rt_hard_coh_correcttrials = rt_hard_coh(find(temp==1)); % rt for correct trials
    correct = sum(temp == 1);
    incorrect = sum(temp == 0);
    total = correct+incorrect;
    pc_hard_coh(i,:) = correct/total*100;
    clear temp correct incorrect total;

    temp = d.correct(1,find(stim_mat(:,8,i)==1)); % d.correct for easy rule trials
    rt_easy_rule_correcttrials = rt_easy_rule(find(temp==1)); % rt for correct trials
    correct = sum(temp == 1);
    incorrect = sum(temp == 0);
    total = correct+incorrect;
    pc_easy_rule(i,:) = correct/total*100;
    clear temp correct incorrect total;
    
    temp = d.correct(1,find(stim_mat(:,8,i)==2)); % d.correct for hard rule trials
    rt_hard_rule_correcttrials = rt_hard_rule(find(temp==1)); % rt for correct trials
    correct = sum(temp == 1);
    incorrect = sum(temp == 0);
    total = correct+incorrect;
    pc_hard_rule(i,:) = correct/total*100;
    clear temp correct incorrect total;
    
    % find the percent correct and rts of correct trials for combinations  
    temp = d.correct(i,find(stim_mat(:,5,i)==1 & stim_mat(:,8,i)==1)); % d.correct for easy coherence and easy rule trials
    rt_easy_coh_easy_rule_correcttrials = rt_easy_coh_easy_rule(find(temp==1)); % rt for correct trials
    correct = sum(temp == 1);
    incorrect = sum(temp == 0);
    total = correct+incorrect;
    pc_easy_coh_easy_rule(i,:) = correct/total*100;
    clear temp correct incorrect total;
    
    temp = d.correct(1,find(stim_mat(:,5,i)==1 & stim_mat(:,8,i)==2));  % d.correct for easy coherence and hard rule trials
    rt_easy_coh_hard_rule_correcttrials = rt_easy_coh_hard_rule(find(temp==1)); % rt for correct trials
    correct = sum(temp == 1);
    incorrect = sum(temp == 0);
    total = correct+incorrect;
    pc_easy_coh_hard_rule(i,:) = correct/total*100;
    clear temp correct incorrect total;
    
    temp = d.correct(1,find(stim_mat(:,5,i)==2 & stim_mat(:,8,i)==1)); % d.correct for hard coherence and easy rule trials
    rt_hard_coh_easy_rule_correcttrials = rt_hard_coh_easy_rule(find(temp==1)); % rt for correct trials
    correct = sum(temp == 1);
    incorrect = sum(temp == 0);
    total = correct+incorrect;
    pc_hard_coh_easy_rule(i,:) = correct/total*100;
    clear temp correct incorrect total;
    
    temp = d.correct(1,find(stim_mat(:,5,i)==2 & stim_mat(:,8,i)==2)); % d.correct for hard coherence and hard rule trials
    rt_hard_coh_hard_rule_correcttrials = rt_hard_coh_hard_rule(find(temp==1)); % rt for correct trials
    correct = sum(temp == 1);
    incorrect = sum(temp == 0);
    total = correct+incorrect;
    pc_hard_coh_hard_rule(i,:) = correct/total*100;
    clear temp correct incorrect total;
    
    % average those
    av_rt_easy_coh_corr(i,:) = nanmean(rt_easy_coh_correcttrials);
    sem_rt_easy_coh_corr(i,:) = nansem(rt_easy_coh_correcttrials);
    av_rt_hard_coh_corr(i,:) = nanmean(rt_hard_coh_correcttrials);
    sem_rt_hard_coh_corr(i,:) = nansem(rt_hard_coh_correcttrials);
    av_rt_easy_rule_corr(i,:) = nanmean(rt_easy_rule_correcttrials);
    sem_rt_easy_rule_corr(i,:) = nansem(rt_easy_rule_correcttrials);
    av_rt_hard_rule_corr(i,:) = nanmean(rt_hard_rule_correcttrials);
    sem_rt_hard_rule_corr(i,:) = nansem(rt_hard_rule_correcttrials);
    av_rt_easy_coh_easy_rule_corr(i,:) = nanmean(rt_easy_coh_easy_rule_correcttrials);
    sem_rt_easy_coh_easy_rule_corr(i,:) = nansem(rt_easy_coh_easy_rule_correcttrials);
    av_rt_easy_coh_hard_rule_corr(i,:) = nanmean(rt_easy_coh_hard_rule_correcttrials);
    sem_rt_easy_coh_hard_rule_corr(i,:) = nansem(rt_easy_coh_hard_rule_correcttrials);
    av_rt_hard_coh_easy_rule_corr(i,:) = nanmean(rt_hard_coh_easy_rule_correcttrials);
    sem_rt_hard_coh_easy_rule_corr(i,:) = nansem(rt_hard_coh_easy_rule_correcttrials);
    av_rt_hard_coh_hard_rule_corr(i,:) = nanmean(rt_hard_coh_hard_rule_correcttrials);
    sem_rt_hard_coh_hard_rule_corr(i,:) = nansem(rt_hard_coh_hard_rule_correcttrials);
    
end

% so now each variable is a list of the nanmean of its thing for each block

av_all_ind = [nanmean(av_rt_easy_coh), nanmean(av_rt_hard_coh), nanmean(av_rt_easy_rule), nanmean(av_rt_hard_rule)];
pc_ind = [nanmean(pc_easy_coh), nanmean(pc_hard_coh), nanmean(pc_easy_rule), nanmean(pc_hard_rule)];
av_corr_ind = [nanmean(av_rt_easy_coh_corr), nanmean(av_rt_hard_coh_corr), nanmean(av_rt_easy_rule_corr), nanmean(av_rt_hard_rule_corr)];
av_all_comb = [nanmean(av_rt_easy_coh_easy_rule), nanmean(av_rt_easy_coh_hard_rule), nanmean(av_rt_hard_coh_easy_rule), nanmean(av_rt_hard_coh_hard_rule)];
pc_comb = [nanmean(pc_easy_coh_easy_rule), nanmean(pc_easy_coh_hard_rule), nanmean(pc_hard_coh_easy_rule), nanmean(pc_hard_coh_hard_rule)];
av_corr_comb = [nanmean(av_rt_easy_coh_easy_rule_corr), nanmean(av_rt_easy_coh_hard_rule_corr), nanmean(av_rt_hard_coh_easy_rule_corr), nanmean(av_rt_hard_coh_hard_rule_corr)];
data_summary = [av_all_ind; pc_ind; av_corr_ind; av_all_comb; pc_comb; av_corr_comb];

sem_summary.easycohcorr = sem_rt_easy_coh_corr;
sem_summary.hardcohcorr = sem_rt_hard_coh_corr;
sem_summary.easyrulecorr = sem_rt_easy_rule_corr;
sem_summary.hardrulecorr = sem_rt_hard_rule_corr;
sem_summary.easycoh_easyrulecorr = sem_rt_easy_coh_easy_rule_corr;
sem_summary.easycoh_hardrulecorr = sem_rt_easy_coh_hard_rule_corr;
sem_summary.hardcoh_easyrulecorr = sem_rt_hard_coh_easy_rule_corr;
sem_summary.hardcoh_hardrulecorr = sem_rt_hard_coh_hard_rule_corr;

end

%% get sem without nans
function semval = nansem(vector_data)
% Recall that s.e.m. = std(x)/sqrt(length(x));
nonan_std = nanstd(vector_data);
nonan_len = length(vector_data(~isnan(vector_data)));
% Plug in values
semval = nonan_std / sqrt(nonan_len);
end
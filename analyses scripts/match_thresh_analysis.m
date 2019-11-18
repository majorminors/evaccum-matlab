%% matching threshold analysis

% note: you're overwriting variables when you produce your plots. You
%       should address this at some point.

% produces:

% 'overview', both for 'low_coh' and 'hi_coh' which is three rows:
%   1) point condition per trial
%   2) correct/incorrect (1/0) per trial
%   3) matching angle (0-90 degrees) per trial
%   note: currently coded such that block one is low coherence and block
%   two is high coherence

% 'summary' which is three rows:
%   1) point condition
%   2) matching angle
%   3) percent correct for low coherence trials
%   4) percent correct for high coherence trials

%  will compute a sinusoid for the data, save parameters as 'sineparams'
%  and produce a figure
%    figure is x(top) - matching angle
%              y      - percent correct

% will save all vars and figs into a folder in data named after the filename which
% starts with the participant number

%% set up

close all;
clearvars;
clc;

% enter filename
data_file = 'S20_EvAccum_matching_threshold_test';
vars = {'p','d'}; % which variables do you need
% enter your thresholds
low_threshold_pc = 0.6;
high_threshold_pc = 0.8; % not using this currently, but we'll have it pop up anyway

% directory mapping
rootdir = 'C:\Users\doria\Google Drive\04 Research\05 Evidence Accumulation\01 EvAccum Code';
addpath(genpath(fullfile(rootdir, 'tools'))) % add sigm_fit tool to path for analysis
datadir = fullfile(rootdir, 'data');
backupdir = fullfile(datadir, 'backup');
if ~exist(fullfile(datadir, data_file),'dir') % make a folder for the paricipant data in 'datadir' if none exist
    mkdir(fullfile(datadir, data_file));
end
if ~exist(backupdir,'dir')
    mkdir(backupdir);
end


% load file
load(fullfile(datadir,[data_file '.mat']),vars{:});
save(fullfile(backupdir, [data_file '_backup'])); % backup subject data before we mess with it

% create matrix specifying stimulus conditions per trial:
%    1)  cue direction (1-4) - evenly allocates trials to cues
%    2)  cue direction in degrees for each trial - evenly adds cue
%        directions to trials in a similar manner to (1)
%    3)  coherence condition (1 = match cue 1, 2 = match cue 2) - evenly allocates
%        trials as either a match cue 1 or match cue 2
%    4)  point condition (1-10) - each repeated
%        p.num_trials_per_block/p.num_points times
%    5)  matching points allocated to point conditions
%    6)  whether trial should add or subtract degrees from cues test (1 =
%        add, 2 = subtract) - currently 2 of each per cue (since four reps
%        of point conditions per cue)
%    7)  coherence direction from cued direction in degrees - calculated from cue direction and matching
%        point
%    8)  coherence direction in degrees incorporating match or non match

% p.rule_points = union([0:5:20],[70:5:90]); % *p.stim_mat* - length(p.rule_points) must == p.num_points

% for fixing S03 and S04 who for non-matches had stimuli presented as
%   matches
% for i = 1:2
%     correct_these(i,:) = d.stim_mat_all(:,3,i)' == 2;
%     for iCorr = 1:160
%         if correct_these(i,iCorr)
%             if strcmp(d.resp_key_name(i,iCorr), d.incorrect_resp(i,iCorr))
%                 d.correct(i,iCorr) = 1;
%             elseif strcmp(d.resp_key_name(i,iCorr), d.correct_resp(i,iCorr))
%                 d.correct(i,iCorr) = 0;
%             end
%         end
%     end
% end
% clear i;


overview_hi_coh = [d.stim_mat_all(:,4,1)';d.correct(1,:);d.stim_mat_all(:,5,1)';d.rt(1,:)];
% note: overview_hi_coh = [point condition; correct/incorrect; matching angle; rts]

overview_low_coh = [d.stim_mat_all(:,4,2)';d.correct(2,:);d.stim_mat_all(:,5,2)';d.rt(2,:)];
% note: overview_low_coh = [point condition; correct/incorrect; matching angle, rts]

rts_low_coh = NaN(160,10);
rts_hi_coh = NaN(160,10);
for i = 1:p.num_points
    corr_per_point_low_coh(i,:) = length(overview_low_coh(3,find(overview_low_coh(1,:)==i & overview_low_coh(2,:)==1)));
    incorr_per_point_low_coh(i,:) = length(overview_low_coh(3,find(overview_low_coh(1,:)==i & overview_low_coh(2,:)==0)));
    temp = overview_low_coh(1,:) == i & overview_low_coh(2,:)==1;
    rts_low_coh(temp,i) = overview_low_coh(4,temp);
    corr_per_point_hi_coh(i,:) = length(overview_hi_coh(3,find(overview_hi_coh(1,:)==i & overview_hi_coh(2,:)==1)));
    incorr_per_point_hi_coh(i,:) = length(overview_hi_coh(3,find(overview_hi_coh(1,:)==i & overview_hi_coh(2,:)==0)));
    temp = overview_hi_coh(1,:) == i & overview_hi_coh(2,:)==1;
    rts_hi_coh(temp,i) = overview_hi_coh(4,temp);
end
pc_low_coh = corr_per_point_low_coh./(corr_per_point_low_coh+incorr_per_point_low_coh)*100;
pc_hi_coh = corr_per_point_hi_coh./(corr_per_point_hi_coh+incorr_per_point_hi_coh)*100;
av_corr_rts_low = nanmean(rts_low_coh);
av_corr_rts_hi = nanmean(rts_hi_coh);
clear ans temp i corr_per_point_low_coh corr_per_point_hi_coh;

summary = [1:p.num_points;p.rule_points;pc_low_coh';av_corr_rts_low;pc_hi_coh';av_corr_rts_hi];
% summary is six rows:
%   1) point condition
%   2) matching angle
%   3) percent correct (low coherence)
%   4) mean rt for correct trials (low coherence)
%   5) percent correct (high coherence)
%   6) mean rt for incorrect trials (high coherence)

save(fullfile(datadir, data_file, [data_file '_summary']), 'summary');

% prep data for psignifit
% Note: the psignifit tool only goes low to high, so if as in this case the
%       values of your stimulus level goes high to low, then you can flip
%       the intensity values (i.e. lie to matlab) or invert them to their
%       negative.
% block/test 1 first
temp=[-d.stim_mat_all(:,5,1) d.correct(1,:)']; % pull required data (matching angle and correct) into temp
temp = sortrows(temp,1); % sort by matching angle
test1 =[]; % open a variable for data
for i = 1:numel(p.rule_points) % for 1:number of matching angle points
    
    idx = temp(:,1)==-p.rule_points(i);
    test1(i,1)=-p.rule_points(i);
    test1(i,2)=sum(temp(idx,2));
    test1(i,3)=sum(idx);
    
    % makes a structure that looks like
    % | matching difficulty point | number correct | number of trials |
    % and builds row-wise each iteration
    
end
%test1(:,1) = flip(test1(:,1)); % need to flip the points (lie to matlab) so the psychcurve tool works

% then block/test 2
temp=[-d.stim_mat_all(:,5,2) d.correct(2,:)']; % pull required data (matching angle and correct) into temp
temp = sortrows(temp,1); % sort by matching angle
test2 =[]; % open a variable for data
for i = 1:numel(p.rule_points) % for 1:number of matching angle points
    
    idx = temp(:,1)==-p.rule_points(i);
    test2(i,1)=-p.rule_points(i);
    test2(i,2)=sum(temp(idx,2));
    test2(i,3)=sum(idx);
    
    % makes a structure that looks like
    % | matching difficulty point | number correct | number of trials |
    % and builds row-wise each iteration
    
end
%test2(:,1) = flip(test2(:,1)); % need to flip the points (lie to matlab) so the psychcurve tool works

% make a sigmoid and put it on a figure
sigmoid_low = figure('visible','on');
[plotresult, plotline, plotdata] = psychcurve(test2);
hold on
% find the x value for proportion correct:
[~, low_threshold_idx] = min(abs(plotline.YData-low_threshold_pc(1))); % first find the index of the value closest to the threshold pc
low_threshold = plotline.XData(low_threshold_idx); % then find the value using the index
high_threshold = -90-low_threshold; % this is the inverse of the low_threshold value - currently I've flipped everything into negative so it is oriented correctly in the psignifit tools
% add plot lines at the threshold value on y:
plot([-90 -0], [low_threshold_pc low_threshold_pc], '-', 'Color',[1 0 0])
[~, high_pc_idx] = min(abs(plotline.XData-high_threshold(1))); % then find the index of the value closest to the high_threshold
high_threshold_pc = plotline.YData(high_pc_idx); % then find the value using the index and print that in the command window
plot([-90 -0], [high_threshold_pc high_threshold_pc], '-', 'Color',[0 1 0])
% add plot lines at the threshold value on x:
plot([low_threshold low_threshold], [0.3 1], '-', 'Color',[1 0 0])
plot([high_threshold high_threshold], [0.3 1], '-', 'Color',[0 1 0])
savefig(fullfile(datadir, data_file, [data_file '_lowcohsigmoid']));
hold off
% diplay rts on a figure
rts_low = figure('visible','off');
plot(summary(2,:),summary(4,:),'ro:');
savefig(fullfile(datadir, data_file, [data_file '_low_coh_rts']));
% make a sigmoid and put it on a figure
sigmoid_hi = figure('visible','off');
[plotresult, plotline, plotdata] = psychcurve(test1);
hold on
[~, low_pc_idx] = min(abs(plotline.XData-low_threshold(1))); % then find the index of the value closest to the high_threshold
low_threshold_pc = plotline.YData(low_pc_idx); % then find the value using the index and print that in the command window
[~, high_pc_idx] = min(abs(plotline.XData-high_threshold(1))); % then find the index of the value closest to the high_threshold
high_threshold_pc = plotline.YData(high_pc_idx); % then find the value using the index and print that in the command window
plot([-90 -0], [low_threshold_pc low_threshold_pc], '-', 'Color',[1 0 0])
plot([-90 -0], [low_threshold_pc low_threshold_pc], '-', 'Color',[1 0 0])
plot([low_threshold low_threshold], [0.3 1], '-', 'Color',[1 0 0])
plot([high_threshold high_threshold], [0.3 1], '-', 'Color',[0 1 0])
savefig(fullfile(datadir, data_file, [data_file '_hicohsigmoid']));
hold off
% diplay rts on a figure
rts_hi = figure('visible','off');
plot(summary(2,:),summary(6,:),'ro:');
savefig(fullfile(datadir, data_file, [data_file '_hi_coh_rts']));

low_threshold = abs(low_threshold) % print that in the command window
high_threshold = abs(high_threshold)
% load those figures into variables
% sigmoid_low=hgload(fullfile(datadir, data_file, [data_file '_lowcohsigmoid.fig']));
% rts_low=hgload(fullfile(datadir, data_file, [data_file '_low_coh_rts.fig']));
% sigmoid_hi=hgload(fullfile(datadir, data_file, [data_file '_hicohsigmoid.fig']));
% rts_hi=hgload(fullfile(datadir, data_file, [data_file '_hi_coh_rts.fig']));

% prepare subplots
figure
visualise(1)=subplot(2,2,1);
visualise(2)=subplot(2,2,2);
visualise(3)=subplot(2,2,3);
visualise(4)=subplot(2,2,4);
% paste our figures on the subplots
copyobj(allchild(get(sigmoid_low,'CurrentAxes')),visualise(1));
copyobj(allchild(get(rts_low,'CurrentAxes')),visualise(2));
copyobj(allchild(get(sigmoid_hi,'CurrentAxes')),visualise(3));
copyobj(allchild(get(rts_hi,'CurrentAxes')),visualise(4));
% add a legend
t(1)=title(visualise(1),'percent correct low coh');
t(2)=title(visualise(2),'reaction time low coh');
t(3)=title(visualise(3),'percent correct hi coh');
t(4)=title(visualise(4),'reaction time hi coh');



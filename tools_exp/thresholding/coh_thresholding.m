function [easy_threshold,hard_threshold,overview,summary] = coh_thresholding(p,d,save_file)
% function [easy_value,hard_value] = coh_thresholding(p,d)
%
% coherence threshold analysis
%
% finds coherence value to achieve a specified percent correct for a
% participant
% 
% requires:
%    various test parameters (expects to find in structure 'p')
%    various saved outputs of test (expects to find in structure 'd')
%    the full path to the saved data (so it can use the same save directory
%       and savename conventions (save_file)
%
% specify in this function your desired percent correct for thresholding (low_threshold_pc,
%   high_threshold_pc)
% 
% produces:
%
% easy_threshold = value to achieve your higher percent correct
% hard_threshold = value to achieve your lower percent correct
%
% 'overview' which is four rows:
%   point condition per trial | correct/incorrect (1/0) per trial | coherence value per trial | rt for each trial
%
% summary is four rows:
%   point condition | coherence value | percent correct | average rt for correct trials
%
%  will compute a sinusoid for the data, save parameters as 'sineparams'
%  and produce a figure
%    figure is x(top) - coherence
%              y      - percent correct
%
% will save figs using the full path save_file

%% set up

% enter your desired thresholds
low_threshold_pc = 0.7; % low coherence/lower percent correct = hard
high_threshold_pc = 0.9; % high coherence/higher percent correct = easy

% create matrix specifying stimulus conditions per trial:
%    1)  cue direction (1-4) - evenly allocates trials to cues
%    2)  cue direction in degrees for each trial - evenly adds cue
%        directions to trials in a similar manner to (1)
%    3)  coherence condition (1 = match cue 1, 2 = match cue 2) - evenly allocates
%        trials as either a match cue 1 or a match cue 2
%    4)  coherence direction in degrees (in relation to cue) (either cue
%        direction, or 180 degrees from cue
%    5)  point condition (1-10) - each repeated
%        p.num_trials_per_block/p.num_points times
%    6)  coherence points allocated to point conditions

%p.coh_points = union([0.1:0.05:0.3], [0.65:0.05:0.85]); % *p.stim_mat* - length(p.coh_points) must == p.num_points

overview = [d.stim_mat_all(:,5,1)';d.correct(1,:);d.stim_mat_all(:,6,1)';d.rt];
% note: overview = [point condition; correct/incorrect; coherence value; rt]

rts = NaN(160,10); 
for i = 1:p.num_points
    corr_per_point(i,:) = length(overview(3,find(overview(1,:)==i & overview(2,:)==1)));
    incorr_per_point(i,:) = length(overview(3,find(overview(1,:)==i & overview(2,:)==0)));
    temp = overview(1,:) == i & overview(2,:)==1;
    rts(temp,i) = overview(4,temp);
end
pc = corr_per_point./(corr_per_point+incorr_per_point)*100;
av_corr_rts = nanmean(rts);
clear ans i temp corr_per_point incorr_per_point;

summary = [1:p.num_points;p.coh_points;pc';av_corr_rts];
% summary is four rows:
%   1) point condition
%   2) coherence value
%   3) percent correct
%   4) average rt for correct trials

% prep data for psignifit
temp=[d.stim_mat_all(:,6) d.correct']; % pull required data (coherence and correct) into temp
temp = sortrows(temp,1); % sort by coherence
temp(find(temp(:,2)==-1),2) = 0; % make invalid responses 0 (incorrect)
data =[]; % open a variable for data
for i = 1:numel(p.coh_points) % for 1:number of coherence points
    
    idx = temp(:,1)==p.coh_points(i);
    data(i,1)=p.coh_points(i);
    data(i,2)=sum(temp(idx,2));
    data(i,3)=sum(idx);
    
    % makes a structure that looks like
    % | coherence point | number correct | number of trials |
    % and builds row-wise each iteration
    
end

[~,low_threshold] = doPsignifit(data,low_threshold_pc);  % new way of getting unscaled thresholds that obviates the need for some of the crazy indexing into the plotline I was doing below
[~,high_threshold] = doPsignifit(data,high_threshold_pc);  % new way of getting unscaled thresholds that obviates the need for some of the crazy indexing into the plotline I was doing below

% make the sigmoid and put it on a figure
sigmoid = figure('visible','on');
psychcurve(data);
hold on
% find the x value for proportion correct:
% [~, low_threshold_idx] = min(abs(plotline.YData-low_threshold_pc(1))); % first find the index of the value closest to the threshold pc
% [~, high_threshold_idx] = min(abs(plotline.YData-high_threshold_pc(1)));
% low_threshold = plotline.XData(low_threshold_idx); % then find the values using the index
% high_threshold = plotline.XData(high_threshold_idx);
% add plot lines at the threshold value on y:
plot([0 1], [low_threshold_pc low_threshold_pc], '-', 'Color',[1 0 0])
plot([0 1], [high_threshold_pc high_threshold_pc], '-', 'Color',[0 1 0])
% add plot lines at the threshold value on x:
plot([low_threshold low_threshold], [0.3 1], '-', 'Color',[1 0 0])
plot([high_threshold high_threshold], [0.3 1], '-', 'Color',[0 1 0])
savefig(sigmoid,[save_file '_sigmoid']);
hold off
% diplay rts on a figure
rts = figure('visible','off');
plot(summary(2,:),summary(4,:),'ro:')
savefig(rts,[save_file '_rts']);
% load those figures into variables
%sigmoid=hgload(fullfile(datadir, save_file, [data_file '_sigmoid.fig']));
%rts=hgload(fullfile(datadir, data_file, [data_file '_rts.fig']));

% prepare subplots
figure
visualise(1)=subplot(1,2,1);
visualise(2)=subplot(1,2,2);
% paste our figures on the subplots
copyobj(allchild(get(sigmoid,'CurrentAxes')),visualise(1));
copyobj(allchild(get(rts,'CurrentAxes')),visualise(2));
% add a legend
t(1)=title(visualise(1),'percent correct');
t(2)=title(visualise(2),'reaction time');

% make the output a little less confusing to understand
easy_threshold = high_threshold;
hard_threshold = low_threshold;

return


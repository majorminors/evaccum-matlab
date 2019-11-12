%% coherence threshold analysis

% produces:

% 'overview' which is four rows:
%   1) point condition per trial
%   2) correct/incorrect (1/0) per trial
%   3) coherence value per trial
%   4) rt for each trial

% summary is four rows:
%   1) point condition
%   2) coherence value
%   3) percent correct
%   4) average rt for correct trials

%  will compute a sinusoid for the data, save parameters as 'sineparams'
%  and produce a figure
%    figure is x(top) - coherence
%              y      - percent correct

% will save all vars and figs into a folder in data named after the filename which
% starts with the participant number

%% set up

close all;
clearvars;
clc;

% enter filename
data_file = 'S20_EvAccum_coherence_threshold_test';
vars = {'p','d'}; % which variables do you need


% directory mapping
rootdir = 'C:\Users\doria\Google Drive\04 Research\05 Evidence Accumulation\01 EvAccum Code';
addpath(genpath(fullfile(rootdir, 'tools'))) % add tools to path for analysis (req. psignifit)
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

save(fullfile(datadir, data_file, [data_file '_summary']), 'summary');

% prep data for psignifit
temp=[d.stim_mat_all(:,6) d.correct']; % pull required data (coherence and correct) into temp
temp = sortrows(temp,1); % sort by coherence
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

set(0,'DefaultFigureVisible','off'); % stop figures from popping up

% make the sigmoid and put it on a figure
psychcurve(data)
savefig(fullfile(datadir, data_file, [data_file '_sigmoid']));
% diplay rts on a figure
plot(summary(2,:),summary(4,:),'ro:')
savefig(fullfile(datadir, data_file, [data_file '_rts']));
% load those figures into variables
sigmoid=hgload(fullfile(datadir, data_file, [data_file '_sigmoid.fig']));
rts=hgload(fullfile(datadir, data_file, [data_file '_rts.fig']));

set(0,'DefaultFigureVisible','on'); % allow figures to pop up again

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


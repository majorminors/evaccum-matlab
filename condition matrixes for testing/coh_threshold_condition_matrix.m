close all;
clearvars;
clc;

% trial settings (*p.stim_mat* = parameter required to calculate stimulus condition matrix)
p.num_trials_per_block = 160; % *p.stim_mat* - must be divisible by p.num_cues && >=15*p.num_points
p.num_cues = 4; % *p.stim_mat*
p.cue_directions = 45:90:315; % *p.stim_mat* - refers to the direction of the upward arrow of the doublesided arrow cue in stimdir
p.num_points = 10; % *p.stim_mat* - number of points to test participants on for each test
p.coh_points = 0.05:0.1:0.95; % *p.stim_mat* - length(p.coh_points) must == p.num_points
if length(p.coh_points) ~= p.num_points; error('number of coherence points is not equal to the number of testing points you specified'); end % check length(p.coh_points) == p.num_points, or error

% create matrix specifying stimulus conditions per trial:
%    1)  cue direction (1-4) - evenly allocates trials to cues
%    2)  blue cue direction in degrees for each trial - evenly adds cue
%        directions to trials in a similar manner to (1)
%    3)  coherence condition (1 = match blue, 2 = match orange) - evenly allocates
%        trials
%    4)  coherence direction in degrees (in relation to cue) (either blue
%        direction, or 180 degrees from blue
%    5)  point condition (1-10) - each repeated
%        p.num_trials_per_block/p.num_points times
%    6)  coherence points allocated to point conditions

% create matrix
p.stim_mat = zeros(p.num_trials_per_block,6);
% create columns
p.stim_mat(:,1) = sort(repmat(1:p.num_cues,[1,p.num_trials_per_block/p.num_cues]));
p.stim_mat(:,2) = p.cue_directions(p.stim_mat(:,1));
p.stim_mat(:,3) = repmat(sort(repmat(1:2,[1,p.num_trials_per_block/p.num_cues/2])),[1,p.num_cues]);
temp = p.stim_mat(:,2)+180; % get 180 degrees on each cue
temp(temp>360) = temp(temp>360)-360; % make any over 360 wrap around from 0 again
p.stim_mat(:,4) = times(p.stim_mat(:,3)-1,temp); % make 1s and 2s from (4) into 0s and 1s, then times by temp (so becomes 0s and cue+180degs)
temp = (p.stim_mat(:,4)==0); % index all posns in (5) with 0s into 'temp'
p.stim_mat(temp,4) = p.stim_mat(temp,2); % insert cue directions from (2) into (5) where there are 0s (using 'temp' index)
p.stim_mat(:,5) = repmat(1:p.num_points,[1,p.num_trials_per_block/p.num_points]);
p.stim_mat(:,6) = repmat(p.coh_points,[1,p.num_trials_per_block/p.num_points]);
% clear floating variables
clear temp;

% shuffle the trial condition order for each block into 'd.stim_mat_all' then re-sort by cue
for block=1:p.num_blocks
    d.stim_mat_all(:,:,block) = p.stim_mat(Shuffle(1:p.num_trials_per_block),:);
    d.stim_mat_all(:,:,block) = sortrows(d.stim_mat_all(:,:,block),1);
end
% clear floating variables
clear block;

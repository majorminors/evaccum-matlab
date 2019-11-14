close all;
clearvars;
clc;

% trial settings (*p.stim_mat* = parameter required to calculate stimulus condition matrix)
p.num_trials_per_block = 160; % *p.stim_mat* - must be divisible by p.num_cues && >= 15*p.num_points
p.num_cues = 4; % *p.stim_mat*
p.cue_directions = 45:90:315; % *p.stim_mat* - refers to the direction of the upward arrow of the doublesided arrow cue in stimdir
p.num_points = 10; % *p.stim_mat* - number of points to test participants on for each test
p.rule_points = union([0:5:20],[70:5:90]); % *p.stim_mat* - length(p.rule_points) must == p.num_points
if length(p.rule_points) ~= p.num_points; error('number of coherence points is not equal to the number of testing points you specified'); end % check length(p.rule_points) == p.num_points, or error

% create matrix specifying stimulus conditions per trial:
%    1)  cue direction (1-4) - evenly allocates trials to cues
%    2)  blue cue direction in degrees for each trial - evenly adds cue
%        directions to trials in a similar manner to (1)
%    3)  coherence condition (1 = match blue, 2 = match orange) - evenly allocates
%        trials
%    4)  point condition (1-10) - each repeated
%        p.num_trials_per_block/p.num_points times
%    5)  matching points allocated to point conditions
%    6)  whether trial should add or subtract degrees from cues test (1 =
%        add, 2 = subtract) - currently 2 of each per cue (since four reps
%        of point conditions per cue)
%    7)  coherence direction from cued direction in degrees - calculated from cue direction and matching
%        point
%    8)  coherence direction in degrees incorporating blue match or orange match

% create matrix
p.stim_mat = zeros(p.num_trials_per_block,8);
% create columns
p.stim_mat(:,1) = sort(repmat(1:p.num_cues,[1,p.num_trials_per_block/p.num_cues]));
p.stim_mat(:,2) = p.cue_directions(p.stim_mat(:,1));
p.stim_mat(:,3) = repmat(sort(repmat(1:2,[1,p.num_trials_per_block/p.num_cues/2])),[1,p.num_cues]);
p.stim_mat(:,4) = repmat(1:p.num_points,[1,p.num_trials_per_block/p.num_points]);
p.stim_mat(:,5) = repmat(p.rule_points,[1,p.num_trials_per_block/p.num_points]);
p.stim_mat(:,6) = repmat(sort(repmat(1:2,[1,p.num_trials_per_block/p.num_cues/4])),[1,p.num_cues*2]);
add = (p.stim_mat(:,6)==1); % index all posns in (6) with 1s into 'add'
subtract = (p.stim_mat(:,6)==2); % index all posns in (6) with 2s into 'subtract'
p.stim_mat(add,7) = p.stim_mat(add,2)+p.stim_mat(add,5); % insert addition of matching point to cue direction into (7)
p.stim_mat(subtract,7) = p.stim_mat(subtract,2)-p.stim_mat(subtract,5); % insert subtraction of matching point to cue direction into (7)
temp = p.stim_mat(:,7)>360; % get all places where (7) > 360 degrees
p.stim_mat(temp,7) = p.stim_mat(temp,7)-360; % make any in (7) over 360 wrap around from 0 again
temp = p.stim_mat(:,7)<0; % get all places where (7) < 0 degrees
p.stim_mat(temp,7) = p.stim_mat(temp,7)+360; % make any in (7) under 0 wrap around from 360 again
temp = p.stim_mat(:,7)+180; % get 180 degrees on each cue
temp(temp>360) = temp(temp>360)-360; % make any over 360 wrap around from 0 again
p.stim_mat(:,8) = times(p.stim_mat(:,3)-1,temp); % make 1s and 2s from (3) into 0s and 1s, then times by temp (so becomes 0s and cue+180degs)
temp = (p.stim_mat(:,8)==0); % index all posns in (8) with 0s into 'temp'
p.stim_mat(temp,8) = p.stim_mat(temp,7); % insert coherence directions from (7) into (8) where there are 0s (using 'temp' index)
% clear floating variables
clear temp add subtract;

% shuffle the trial condition order for each block into 'd.stim_mat_all' then re-sort by cue
for block=1:p.num_blocks
    d.stim_mat_all(:,:,block) = p.stim_mat(Shuffle(1:p.num_trials_per_block),:);
    d.stim_mat_all(:,:,block) = sortrows(d.stim_mat_all(:,:,block),1);
end
% clear floating variables
clear block;

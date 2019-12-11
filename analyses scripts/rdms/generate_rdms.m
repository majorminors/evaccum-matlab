%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);

rootdir = '\\cbsu\data\Group\Woolgar-Lab\projects\EvAccum'; % root directory

%% do stimulus matrix

% test rule entry values (entered as prompted in script)
d.easy_rule = 80;
d.hard_rule = 10;

% how many blocks
p.num_blocks = 20;

% trial settings (*p.stim_mat* = parameter required to calculate stimulus condition matrix)
p.num_trials_per_block = 64; % *p.stim_mat*
p.num_cues = 4; % *p.stim_mat*
p.num_motion_coherence = 8; % *p.stim_mat*
p.cue_directions = 45:90:315; % *p.stim_mat* - refers to the direction of the upward arrow of the doublesided arrow cue in stimdir
p.dot_motion_directions = union([p.cue_directions+d.easy_rule],[p.cue_directions+d.hard_rule]); % *p.stim_mat* - adds easy rule and hard rule to each cue, then puts them in a vector sorted low to high

% create matrix specifying stimulus conditions per trial:
%  1)  cue direction (1-4)
%  2)  cue direction in degrees
%  3)  dot motion direction condition (1-8)
%  4)  dot motion direction in degrees
%  5)  coherence difficulty (1 = easy, 2 = hard)
%  6)  matching distance from cue direction (absolute value) - used to calc match and match difficulty
%  7)  match blue (1) or match orange (2)
%  8)  matching difficulty (1 = easy, 2 = difficult)
%  9)  gives you a unique number for each trial condition
% 10)  gives a number based on 9 for meg triggers
% note: if you try to test with two matching angles that are the same, you
%       will get an error, so make them different by at least 1 degree. this is
%       because we use 'union()' to calc matching distance from cue direction.

% create matrix
p.stim_mat = zeros(p.num_trials_per_block,10);

% create columns
p.stim_mat(:,1) = sort(repmat(1:p.num_cues,[1,p.num_trials_per_block/p.num_cues]));
p.stim_mat(:,2) = p.cue_directions(p.stim_mat(:,1));
p.stim_mat(:,3) = repmat(sort(repmat(1:p.num_motion_coherence,[1,p.num_trials_per_block/p.num_cues/p.num_motion_coherence])),[1,p.num_cues]);
p.stim_mat(:,4) = p.dot_motion_directions(p.stim_mat(:,3));
p.stim_mat(:,5) = repmat(sort(repmat(1:2,[1,p.num_trials_per_block/p.num_cues/p.num_motion_coherence/2])),[1,p.num_cues*p.num_motion_coherence]);
dist = abs([p.stim_mat(:,4)-360,p.stim_mat(:,4),p.stim_mat(:,4)+360]-repmat(p.stim_mat(:,2),[1,3]));
p.stim_mat(:,6) = min(dist,[],2);
p.stim_mat(:,7) = (p.stim_mat(:,6)>90)+1;
p.stim_mat(:,8) = ~((p.stim_mat(:,6)==min(p.stim_mat(:,6)))|(p.stim_mat(:,6)==max(p.stim_mat(:,6))))+1;
p.stim_mat(:,9) = 1:length(p.stim_mat(:,9));
p.stim_mat(:,10) = p.stim_mat(:,9)+[0:2:127]'+5;

% clear floating variables
clear dist;

% shuffle the trial condition order for each block into 'd.stim_mat_all' then re-sort by cue
% for block=1:p.num_blocks
%     d.stim_mat_all(:,:,block) = p.stim_mat(Shuffle(1:p.num_trials_per_block),:);
%     d.stim_mat_all(:,:,block) = sortrows(d.stim_mat_all(:,:,block),1);
% end
% clear floating variables
clear block;

%% RDM - stimulus

mat1 = repmat(p.stim_mat(:,4),1,64);

mat2 = repmat(p.stim_mat(:,4)',64,1);

rdm_stim = mat1~=mat2;

clear mat1 mat2

fig.stim = imagesc(rdm_stim);

%% RDM - rule

mat1 = repmat(p.stim_mat(:,6),1,64);

mat2 = repmat(p.stim_mat(:,6)',64,1);

rdm_rule = mat1~=mat2;

clear mat1 mat2

fig.rule = imagesc(rdm_rule);

%% RDM - cue

mat1 = repmat(p.stim_mat(:,1),1,64);

mat2 = repmat(p.stim_mat(:,1)',64,1);

rdm_cue = mat1~=mat2;

clear mat1 mat2

fig.cue = imagesc(rdm_cue);

%% RDM - response

mat1 = repmat(p.stim_mat(:,7),1,64);

mat2 = repmat(p.stim_mat(:,7)',64,1);

rdm_resp = mat1~=mat2;

clear mat1 mat2

fig.resp = imagesc(rdm_resp);

%% show all figures

figure
visualise(1)=subplot(2,2,1);
visualise(2)=subplot(2,2,2);
visualise(3)=subplot(2,2,3);
visualise(4)=subplot(2,2,4);
% paste our figures on the subplots
copyobj(imagesc(rdm_stim),visualise(1));
copyobj(imagesc(rdm_cue),visualise(2));
copyobj(imagesc(rdm_rule),visualise(3));
copyobj(imagesc(rdm_resp),visualise(4));

linkaxes([visualise],'xy')
% add a legend
t(1)=title(visualise(1),'stim rdm');
t(2)=title(visualise(2),'cue rdm');
t(3)=title(visualise(3),'rule rdm');
t(4)=title(visualise(4),'resp rdm');
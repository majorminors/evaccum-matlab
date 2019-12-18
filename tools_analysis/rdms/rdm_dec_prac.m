%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);

rootdir = '\\cbsu\data\Group\Woolgar-Lab\projects\Dorian\EvAccum'; % root directory

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
temp = p.stim_mat(:,4);
temp(temp > 360) = temp(temp > 360) - 360;
p.stim_mat(:,4) = temp; clear temp;
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

%%

mat1 = repmat(p.stim_mat(:,1),1,64);

mat2 = repmat(p.stim_mat(:,1)',64,1);

for i1 = 1:length(mat1(:,1))
    for i2 = 1:length(mat2(1,:))
        if mat1(i1,i2) == mat2(i1,i2)
            rdm_dec2(i1,i2) = 0;
        elseif mat1(i1,i2) == 1 || mat1(i1,i2) == 3 && mat2(i1,i2) == 2 || mat2(i1,i2) == 4
            rdm_dec2(i1,i2) = 1;
        elseif mat2(i1,i2) == 1 || mat2(i1,i2) == 3 && mat1(i1,i2) == 2 || mat1(i1,i2) == 4
            rdm_dec2(i1,i2) = 1;
        elseif mat1(i1,i2) == 1 && mat2(i1,i2) == 3 || mat2(i1,i2) == 1 && mat1(i1,i2) == 3
            rdm_dec2(i1,i2) = 2;
        elseif mat1(i1,i2) == 2 && mat2(i1,i2) == 4 || mat2(i1,i2) == 2 && mat1(i1,i2) == 4
            rdm_dec2(i1,i2) = 2;
        else
            error('didnt catch em all');
        end
    end
end
    

for i1 = 1:length(mat1(:,1))
    for i2 = 1:length(mat1(1,:))
        if mat1(i1,i2) == mat2(i1,i2)
            rdm_dec(i1,i2) = 0;
        elseif mat1(i1,i2) == 1
            if mat2(i1,i2) == 2 || 4
                rdm_dec(i1,i2) = 1;
            elseif mat2(i1,i2) == 3
                rdm_dec(i1,i2) = 2;
            end
        elseif mat1(i1,i2) == 2
            if mat2(i1,i2) == 1 || 3
                rdm_dec(i1,i2) = 1;
            elseif mat2(i1,i2) == 4
                rdm_dec(i1,i2) = 2;
            end
        elseif mat1(i1,i2) == 3
            if mat2(i1,i2) == 2 || 4
                rdm_dec(i1,i2) = 1;
            elseif mat2(i1,i2) == 1
                rdm_dec(i1,i2) = 2;
            end
        elseif mat1(i1,i2) == 4
            if mat2(i1,i2) == 1 || 3
                rdm_dec(i1,i2) = 1;
            elseif mat2(i1,i2) == 2
                rdm_dec(i1,i2) = 2;
            end
        else
            error('didnt catch em all');
        end
    end
end
    
    

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

%% RDM - stimulus distance

mat1 = repmat(p.stim_mat(:,4),1,64); % repeat column 3 (cue direction code) for trial columns

mat2 = repmat(p.stim_mat(:,4)',64,1); % repeat column 3 transposed for trial rows

mat3 = mod(mat1-mat2,360); % normalised difference between angles

rdm_stim = min(360-mat3,mat3); % the smallest of the two differences between angles

clear mat1 mat2 mat3

imagesc(rdm_stim);
savefig('rdm_stim')

% then you need to do two more:
% 1. decisions made - so 1-4 cue one is different from 3-6 cue 2 and 5-8 cue 3 and 7-2 cue 4
% 2. decisions made detail - so 1-4 cue one is LESS different to 3-6 cue 2 THAN 5-8 cue 3

%% RDM - coherence

mat1 = repmat(p.stim_mat(:,5),1,64);

mat2 = repmat(p.stim_mat(:,5)',64,1);

rdm_coh = mat1~=mat2; % if not equal, 1 (dissimilar) else 0

clear mat1 mat2

imagesc(rdm_coh);
savefig('rdm_coh')

%% RDM - cue

mat1 = repmat(p.stim_mat(:,1),1,64); % repeat cue code for trials as columns

mat2 = repmat(p.stim_mat(:,1)',64,1); % repeat cue code for trials as rows

rdm_cue = mat1~=mat2; % if they aren't equal, then 1 (dissimilar) otherwise 0

imagesc(rdm_cue);
savefig('rdm_cue');

%% RDM - decision boundary

rdm_decbdry = abs(mat1-mat2); % minus one from the other to get dissimilarity

rdm_decbdry(rdm_decbdry == 3) = 1;

clear mat1 mat2

imagesc(rdm_decbdry);
savefig('rdm_ruledec');

%% RDM - rule difficulty

mat1 = repmat(p.stim_mat(:,8),1,64); % repeat difficulty for trial cols

mat2 = repmat(p.stim_mat(:,8)',64,1); % repeat difficulty for trial rows

rdm_rulediff = mat1~=mat2; % if not equal, 1 (dissimilar) else 0

clear mat1 mat2

imagesc(rdm_rulediff);
savefig('rdm_rulediff');

%% RDM - response

mat1 = repmat(p.stim_mat(:,7),1,64); % repeat for trial columns

mat2 = repmat(p.stim_mat(:,7)',64,1); % repeat for trial rows

rdm_resp = mat1~=mat2; % if not equal, 1 (dissimilar) else 0

clear mat1 mat2

imagesc(rdm_resp);
savefig('rdm_resp');

%% RDM - decision

rdm_dec = rdm_decbdry;

for i1 = 1:length(rdm_dec(:,1))
    for i2 = 1:length(rdm_dec(1,:))
        if rdm_resp(i1,i2) == 1 && rdm_dec(i1,i2) == 0
            rdm_dec(i1,i2) = 2;
        elseif rdm_resp(i1,i2) == 1 &&  rdm_dec(i1,i2) == 2
            rdm_dec(i1,i2) = 0;
        end
    end
end

rdm_dec_equal = rdm_decbdry;

for i1 = 1:length(rdm_dec_equal(:,1))
    for i2 = 1:length(rdm_dec_equal(1,:))
        if rdm_resp(i1,i2) == 1 && rdm_dec_equal(i1,i2) == 0
            rdm_dec_equal(i1,i2) = 2;
        elseif rdm_resp(i1,i2) == 1 &&  rdm_dec_equal(i1,i2) == 2
            rdm_dec_equal(i1,i2) = 1;
        end
    end
end

imagesc(rdm_dec);
savefig('rdm_dec');

imagesc(rdm_dec_equal);
savefig('rdm_dec_equal');

%% show all figures

% create subplot and map rdms to it
figure
visualise(1)=subplot(4,2,1);
imagesc(rdm_stim);
visualise(2)=subplot(4,2,2);
imagesc(rdm_coh);
visualise(3)=subplot(4,2,3);
imagesc(rdm_cue);
visualise(4)=subplot(4,2,4);
imagesc(rdm_decbdry);
visualise(5)=subplot(4,2,5);
imagesc(rdm_rulediff);
visualise(6)=subplot(4,2,6);
imagesc(rdm_resp);
visualise(7)=subplot(4,2,7);
imagesc(rdm_dec);
visualise(8)=subplot(4,2,8);
imagesc(rdm_dec_equal);

% add a legend
t(1)=title(visualise(1),'stimulus');
t(5)=title(visualise(2),'coherence difficulty');
t(2)=title(visualise(3),'cue');
t(3)=title(visualise(4),'decision boundary');
t(4)=title(visualise(5),'rule difficulty');
t(5)=title(visualise(6),'response');
t(5)=title(visualise(7),'decision');
t(5)=title(visualise(8),'decision treated equal');
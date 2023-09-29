%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';
toolsdir = fullfile(rootdir,'tools_analysis');
modeldir = fullfile(toolsdir,'rdms');
savefigs = 0;

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

easy_coh = (find(p.stim_mat(:,5) == 1));
hard_coh = (find(p.stim_mat(:,5) == 2));
easy_dec = (find(p.stim_mat(:,8) == 1));
hard_dec = (find(p.stim_mat(:,8) == 2));

eced = easy_coh(ismember(easy_coh,easy_dec));
echd = easy_coh(ismember(easy_coh,hard_dec));
hced = hard_coh(ismember(hard_coh,easy_dec));
hchd = hard_coh(ismember(hard_coh,hard_dec));

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
% minimum angle between cues

mat1 = repmat(p.stim_mat(:,4),1,64); % repeat column 3 (cue direction code) for trial columns

mat2 = repmat(p.stim_mat(:,4)',64,1); % repeat column 3 transposed for trial rows

mat3 = mod(mat1-mat2,360); % normalised difference between angles

rdm_stim = min(360-mat3,mat3); % the smallest of the two differences between angles

clear mat1 mat2 mat3

makeModel(rdm_stim,'rdm_stim',modeldir,savefigs)

%% RDM - coherence
% coherence level same or different

mat1 = repmat(p.stim_mat(:,5),1,64); % repeat coherence for trial cols

mat2 = repmat(p.stim_mat(:,5)',64,1); % repeat coherence for trial rows

rdm_coh = mat1~=mat2; % if not equal, 1 (dissimilar) else 0

clear mat1 mat2

makeModel(rdm_coh,'rdm_coh',modeldir,savefigs)

%% RDM - cue
% cue as same or different

mat1 = repmat(p.stim_mat(:,1),1,64); % repeat cue code for trials as columns

mat2 = repmat(p.stim_mat(:,1)',64,1); % repeat cue code for trials as rows

rdm_cue_simple = mat1~=mat2; % if they aren't equal, then 1 (dissimilar) otherwise 0

makeModel(rdm_cue_simple,'rdm_cue_simple',modeldir,savefigs)

%% RDM - cue detail
% same and opposite cues different, and medial are equal

mat1 = repmat(p.stim_mat(:,1),1,64); % repeat cue code for trials as columns

mat2 = repmat(p.stim_mat(:,1)',64,1); % repeat cue code for trials as rows

rdm_cue_detail = abs(mat1-mat2); % minus one from the other to get dissimilarity

rdm_cue_detail(rdm_cue_detail == 3) = 1; % the furthest distant are actually closer again

clear mat1 mat2

makeModel(rdm_cue_detail,'rdm_cue_detail',modeldir,savefigs)

%% RDM - decision boundary
% same and opposite cues same, and medial are equally different

mat1 = repmat(p.stim_mat(:,1),1,64); % repeat cue code for trials as columns

mat2 = repmat(p.stim_mat(:,1)',64,1); % repeat cue code for trials as rows

rdm_decbdry = abs(mat1-mat2); % minus one from the other to get dissimilarity

rdm_decbdry(rdm_decbdry == 3) = 1;
rdm_decbdry(rdm_decbdry == 2) = 0;

clear mat1 mat2

makeModel(rdm_decbdry,'rdm_decbdry',modeldir,savefigs)

%% RDM - rule difficulty
% rule difficulty same or different
mat1 = repmat(p.stim_mat(:,8),1,64); % repeat difficulty for trial cols

mat2 = repmat(p.stim_mat(:,8)',64,1); % repeat difficulty for trial rows

rdm_rulediff = mat1~=mat2; % if not equal, 1 (dissimilar) else 0

clear mat1 mat2

makeModel(rdm_rulediff,'rdm_rulediff',modeldir,savefigs)

%% RDM - response
% response button same or different

mat1 = repmat(p.stim_mat(:,7),1,64); % repeat for trial columns

mat2 = repmat(p.stim_mat(:,7)',64,1); % repeat for trial rows

rdm_resp = mat1~=mat2; % if not equal, 1 (dissimilar) else 0

clear mat1 mat2

makeModel(rdm_resp,'rdm_resp',modeldir,savefigs)

%% RDM - decision
% same or different side of the decision boundary

% then you need to do two more:
% 1. decisions made - so 1-4 cue one is different from 3-6 cue 2 and 5-8 cue 3 and 7-2 cue 4
% 2. decisions made detail - so 1-4 cue one is LESS different to 3-6 cue 2 THAN 5-8 cue 3

rdm_dec = rdm_cue_detail; % pull in the cue information (same, medial, or opposite) so we can use the response as a proxy for all decisions on one half of the dec bndry

for i1 = 1:length(rdm_dec(:,1))
    for i2 = 1:length(rdm_dec(1,:))
        if rdm_resp(i1,i2) == 1 && rdm_dec(i1,i2) == 0 % when the cue is the same but the response (proxy for direction) is on the other side of the decision boundary
            rdm_dec(i1,i2) = 2; % code as opposite
        elseif rdm_resp(i1,i2) == 1 &&  rdm_dec(i1,i2) == 2 % when the cue is opposite and the response (proxy for direction) is on the other side of the decision boundary
            rdm_dec(i1,i2) = 0; % code as the same
        elseif rdm_dec(i1,i2) == 1
                rdm_dec(i1,i2) = 1; % otherwise make no assumption about distance (should be NaNs, but code as number for visual)
        end
    end
end

makeModel(rdm_dec,'rdm_dec',modeldir,savefigs)

rdm_dec_diff = rdm_cue_simple; % pull in the cue information (same or different) so we can use the response as a proxy for all decisions on one half of the dec bndry
rdm_dec_diff = rdm_dec_diff-1; % make different trials 0 and same trials -1
for i1 = 1:length(rdm_dec_diff(:,1))
    for i2 = 1:length(rdm_dec_diff(1,:))
        if rdm_dec_diff(i1,i2) == -1 && rdm_resp(i1,i2) == 1 % if the response is a 1 on same trials, make it a 1
            rdm_dec_diff(i1,i2) = 1;
        end
    end
end

makeModel(rdm_dec_diff,'rdm_dec_diff',modeldir,savefigs)

%% extract specific trials

rdm_dec_eced = rdm_dec(eced,eced);
rdm_dec_eced_tri = triu(rdm_dec_eced,1); % get upper triangle above the first diagonal
%corrcoef(x)
rdm_dec_echd = rdm_dec(echd,echd);
rdm_dec_echd_tri = triu(rdm_dec_echd,1); % get upper triangle above the first diagonal
rdm_dec_hced = rdm_dec(hced,hced);
rdm_dec_hced_tri = triu(rdm_dec_hced,1); % get upper triangle above the first diagonal
rdm_dec_hchd = rdm_dec(hchd,hchd);
rdm_dec_hchd_tri = triu(rdm_dec_hchd,1); % get upper triangle above the first diagonal

rdm_stim_ec = rdm_stim(easy_coh,easy_coh);
makeModel(rdm_stim_ec,'rdm_stim_ec',modeldir,savefigs);
rdm_stim_hc = rdm_stim(hard_coh,hard_coh);
makeModel(rdm_stim_hc,'rdm_stim_hc',modeldir,savefigs);
rdm_stim_ed = rdm_stim(easy_dec,easy_dec);
makeModel(rdm_stim_ed,'rdm_stim_ed',modeldir,savefigs);
rdm_stim_hd = rdm_stim(hard_dec,hard_dec);
makeModel(rdm_stim_hd,'rdm_stim_hd',modeldir,savefigs);
rdm_coh_ec = rdm_coh(easy_coh,easy_coh);
makeModel(rdm_coh_ec,'rdm_coh_ec',modeldir,savefigs);
rdm_coh_hc = rdm_coh(hard_dec,hard_dec);
makeModel(rdm_coh_hc,'rdm_coh_hc',modeldir,savefigs);
rdm_coh_ed = rdm_coh(easy_dec,easy_dec);
makeModel(rdm_coh_ed,'rdm_coh_ed',modeldir,savefigs);
rdm_coh_hd = rdm_coh(hard_dec,hard_dec);
makeModel(rdm_coh_hd,'rdm_coh_hd',modeldir,savefigs);
rdm_cue_simple_ec = rdm_cue_simple(easy_coh,easy_coh);
makeModel(rdm_cue_simple_ec,'rdm_cue_simple_ec',modeldir,savefigs);
rdm_cue_simple_hc = rdm_cue_simple(hard_dec,hard_dec);
makeModel(rdm_cue_simple_hc,'rdm_cue_simple_hc',modeldir,savefigs);
rdm_cue_simple_ed = rdm_cue_simple(easy_dec,easy_dec);
makeModel(rdm_cue_simple_ed,'rdm_cue_simple_ed',modeldir,savefigs);
rdm_cue_simple_hd = rdm_cue_simple(hard_dec,hard_dec);
makeModel(rdm_cue_simple_hd,'rdm_cue_simple_hd',modeldir,savefigs);
rdm_cue_detail_ec = rdm_cue_detail(easy_coh,easy_coh);
makeModel(rdm_cue_detail_ec,'rdm_cue_detail_ec',modeldir,savefigs);
rdm_cue_detail_hc = rdm_cue_detail(hard_dec,hard_dec);
makeModel(rdm_cue_detail_hc,'rdm_cue_detail_hc',modeldir,savefigs);
rdm_cue_detail_ed = rdm_cue_detail(easy_dec,easy_dec);
makeModel(rdm_cue_detail_ed,'rdm_cue_detail_ed',modeldir,savefigs);
rdm_cue_detail_hd = rdm_cue_detail(hard_dec,hard_dec);
makeModel(rdm_cue_detail_hd,'rdm_cue_detail_hd',modeldir,savefigs);
rdm_decbdry_ec = rdm_decbdry(easy_coh,easy_coh);
makeModel(rdm_decbdry_ec,'rdm_decbdry_ec',modeldir,savefigs);
rdm_decbdry_hc = rdm_decbdry(hard_dec,hard_dec);
makeModel(rdm_decbdry_hc,'rdm_decbdry_hc',modeldir,savefigs);
rdm_decbdry_ed = rdm_decbdry(easy_dec,easy_dec);
makeModel(rdm_decbdry_ed,'rdm_decbdry_ed',modeldir,savefigs);
rdm_decbdry_hd = rdm_decbdry(hard_dec,hard_dec);
makeModel(rdm_decbdry_hd,'rdm_decbdry_hd',modeldir,savefigs);
rdm_rulediff_ec = rdm_rulediff(easy_coh,easy_coh);
makeModel(rdm_rulediff_ec,'rdm_rulediff_ec',modeldir,savefigs);
rdm_rulediff_hc = rdm_rulediff(hard_dec,hard_dec);
makeModel(rdm_rulediff_hc,'rdm_rulediff_hc',modeldir,savefigs);
rdm_rulediff_ed = rdm_rulediff(easy_dec,easy_dec);
makeModel(rdm_rulediff_ed,'rdm_rulediff_ed',modeldir,savefigs);
rdm_rulediff_hd = rdm_rulediff(hard_dec,hard_dec);
makeModel(rdm_rulediff_hd,'rdm_rulediff_hd',modeldir,savefigs);
rdm_resp_ec = rdm_resp(easy_coh,easy_coh);
makeModel(rdm_resp_ec,'rdm_resp_ec',modeldir,savefigs);
rdm_resp_hc = rdm_resp(hard_dec,hard_dec);
makeModel(rdm_resp_hc,'rdm_resp_hc',modeldir,savefigs);
rdm_resp_ed = rdm_resp(easy_dec,easy_dec);
makeModel(rdm_resp_ed,'rdm_resp_ed',modeldir,savefigs);
rdm_resp_hd = rdm_resp(hard_dec,hard_dec);
makeModel(rdm_resp_hd,'rdm_resp_hd',modeldir,savefigs);
rdm_dec_ec = rdm_dec(easy_coh,easy_coh);
makeModel(rdm_dec_ec,'rdm_dec_ec',modeldir,savefigs);
rdm_dec_hc = rdm_dec(hard_dec,hard_dec);
makeModel(rdm_dec_hc,'rdm_dec_hc',modeldir,savefigs);
rdm_dec_ed = rdm_dec(easy_dec,easy_dec);
makeModel(rdm_dec_ed,'rdm_dec_ed',modeldir,savefigs);
rdm_dec_hd = rdm_dec(hard_dec,hard_dec);
makeModel(rdm_dec_hd,'rdm_dec_hd',modeldir,savefigs);
rdm_dec_diff_ec = rdm_dec_diff(easy_coh,easy_coh);
makeModel(rdm_dec_diff_ec,'rdm_dec_diff_ec',modeldir,savefigs);
rdm_dec_diff_hc = rdm_dec_diff(hard_dec,hard_dec);
makeModel(rdm_dec_diff_hc,'rdm_dec_diff_hc',modeldir,savefigs);
rdm_dec_diff_ed = rdm_dec_diff(easy_dec,easy_dec);
makeModel(rdm_dec_diff_ed,'rdm_dec_diff_ed',modeldir,savefigs);
rdm_dec_diff_hd = rdm_dec_diff(hard_dec,hard_dec);
makeModel(rdm_dec_diff_hd,'rdm_dec_diff_hd',modeldir,savefigs);

figure
visualise(1)=subplot(2,2,1);
imagesc(rdm_dec_eced);
visualise(2)=subplot(2,2,2);
imagesc(rdm_dec_echd);
visualise(3)=subplot(2,2,3);
imagesc(rdm_dec_hced);
visualise(4)=subplot(2,2,4);
imagesc(rdm_dec_hchd);
t(1)=title(visualise(1),'eced');
t(2)=title(visualise(2),'echd');
t(3)=title(visualise(3),'hced');
t(4)=title(visualise(4),'hchd');

%% show all figures

% create subplot and map rdms to it
figure
visualise(1)=subplot(2,3,1);
imagesc(rdm_stim);
visualise(2)=subplot(2,3,2);
imagesc(rdm_cue_detail);
visualise(3)=subplot(2,3,3);
imagesc(rdm_decbdry);
visualise(4)=subplot(2,3,4);
imagesc(rdm_resp);
visualise(5)=subplot(2,3,5);
imagesc(rdm_dec);
visualise(6)=subplot(2,3,6);
imagesc(rdm_dec_diff);


% add a legend
t(1)=title(visualise(1),'stimulus');
t(5)=title(visualise(2),'cue');
t(2)=title(visualise(3),'decision boundary');
t(3)=title(visualise(4),'button press');
t(4)=title(visualise(5),'motion classification detail');
t(5)=title(visualise(6),'motion classification simple');

%% subfunctions

function makeModel(model,savename,savefolder,savefigs)

imagesc(model);
if savefigs
    print([savefolder filesep savename '.png'], '-dpng');
    savefig(savename)
end
csvwrite([savefolder filesep savename '.csv'], model);

return
end

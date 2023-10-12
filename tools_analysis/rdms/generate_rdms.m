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

trialId = p.stim_mat(:,9);
save([modeldir filesep 'trialIds.mat'],'trialId');

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
% same or different side of the decision boundary, but perpendicular decision boundaries are NaN (no prediction) 

rdm_dec_detail_null = rdm_cue_detail; % pull in the cue information (same, medial, or opposite) so we can use the response as a proxy for all decisions on one half of the dec bndry

for i1 = 1:length(rdm_dec_detail_null(:,1))
    for i2 = 1:length(rdm_dec_detail_null(1,:))
        if rdm_resp(i1,i2) == 1 && rdm_dec_detail_null(i1,i2) == 0 % when the cue is the same but the response (proxy for direction) is on the other side of the decision boundary
            rdm_dec_detail_null(i1,i2) = 2; % code as opposite
        elseif rdm_resp(i1,i2) == 1 &&  rdm_dec_detail_null(i1,i2) == 2 % when the cue is opposite and the response (proxy for direction) is on the other side of the decision boundary
            rdm_dec_detail_null(i1,i2) = 0; % code as the same
        elseif rdm_dec_detail_null(i1,i2) == 1
                rdm_dec_detail_null(i1,i2) = NaN; % otherwise make no assumption about distance
        end
    end
end

makeModel(rdm_dec_detail_null,'rdm_dec_detail_null',modeldir,savefigs)

%% RDM - decision 2
% same or different side of the decision boundary, and perpendicular
% decision boundaries are equally different from both

rdm_dec_detail_pred = rdm_cue_detail; % pull in the cue information (same, medial, or opposite) so we can use the response as a proxy for all decisions on one half of the dec bndry

for i1 = 1:length(rdm_dec_detail_pred(:,1))
    for i2 = 1:length(rdm_dec_detail_pred(1,:))
        if rdm_resp(i1,i2) == 1 && rdm_dec_detail_pred(i1,i2) == 0 % when the cue is the same but the response (proxy for direction) is on the other side of the decision boundary
            rdm_dec_detail_pred(i1,i2) = 2; % code as opposite
        elseif rdm_resp(i1,i2) == 1 &&  rdm_dec_detail_pred(i1,i2) == 2 % when the cue is opposite and the response (proxy for direction) is on the other side of the decision boundary
            rdm_dec_detail_pred(i1,i2) = 0; % code as the same
        elseif rdm_dec_detail_pred(i1,i2) == 1
                rdm_dec_detail_pred(i1,i2) = 1; % otherwise assume they are equally different
        end
    end
end

makeModel(rdm_dec_detail_pred,'rdm_dec_detail_pred',modeldir,savefigs)

%% RDM - decision 3
% same or different side of the decision boundary when cue is the same, but
% perpendicular decision boundaries and the same decision boundary but
% associated with the opposite cue are different

rdm_dec_simple = rdm_cue_simple; % pull in the cue information (same or different)
rdm_dec_simple = rdm_dec_simple-1; % make different trials 0 and same trials -1
for i1 = 1:length(rdm_dec_simple(:,1))
    for i2 = 1:length(rdm_dec_simple(1,:))
        if rdm_dec_simple(i1,i2) == -1 && rdm_resp(i1,i2) == 1 % if the response is a 1 on same trials, make it a 1
            rdm_dec_simple(i1,i2) = 1;
        end
    end
end

% but cosmo doesn't like that the diagonals aren't zero, so let's shift it
% all up one, so -1 is 0, opposite is 2, and equally different is 1
rdm_dec_simple = rdm_dec_simple+1;

makeModel(rdm_dec_simple,'rdm_dec_simple',modeldir,savefigs)

%% extract specific trials
% we don't need to do this: we have a function that generates a trialwise rdm for the trials a subject actually did
% might be useful for figures though

%rdm_stim_eced = rdm_stim(eced,eced);
%rdm_stim_eced_tri = triu(rdm_stim_eced,1); % get upper triangle above the first diagonal
%%corrcoef(x)
%rdm_stim_echd = rdm_stim(echd,echd);
%rdm_stim_echd_tri = triu(rdm_stim_echd,1); % get upper triangle above the first diagonal
%rdm_stim_hced = rdm_stim(hced,hced);
%rdm_stim_hced_tri = triu(rdm_stim_hced,1); % get upper triangle above the first diagonal
%rdm_stim_hchd = rdm_stim(hchd,hchd);
%rdm_stim_hchd_tri = triu(rdm_stim_hchd,1); % get upper triangle above the first diagonal

%rdm_stim_ec = rdm_stim(easy_coh,easy_coh);
%makeModel(rdm_stim_ec,'rdm_stim_ec',modeldir,savefigs);
%rdm_stim_hc = rdm_stim(hard_coh,hard_coh);
%makeModel(rdm_stim_hc,'rdm_stim_hc',modeldir,savefigs);
%rdm_stim_ed = rdm_stim(easy_dec,easy_dec);
%makeModel(rdm_stim_ed,'rdm_stim_ed',modeldir,savefigs);
%rdm_stim_hd = rdm_stim(hard_dec,hard_dec);
%makeModel(rdm_stim_hd,'rdm_stim_hd',modeldir,savefigs);

%figure
%visualise(1)=subplot(2,2,1);
%imagesc(rdm_dec_eced);
%visualise(2)=subplot(2,2,2);
%imagesc(rdm_dec_echd);
%visualise(3)=subplot(2,2,3);
%imagesc(rdm_dec_hced);
%visualise(4)=subplot(2,2,4);
%imagesc(rdm_dec_hchd);
%t(1)=title(visualise(1),'eced');
%t(2)=title(visualise(2),'echd');
%t(3)=title(visualise(3),'hced');
%t(4)=title(visualise(4),'hchd');

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
imagesc(rdm_dec_detail_pred);
visualise(6)=subplot(2,3,6);
imagesc(rdm_dec_simple);


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

%% set up

close all;
clearvars;
clc;

% enter filename
datafilenames = ["S01_EvAccum_stitched", "S02_EvAccum", "S03_EvAccum", "S04_EvAccum", "S05_EvAccum",...
    "S06_EvAccum", "S07_EvAccum", "S08_EvAccum", "S09_EvAccum", "S10_EvAccum", "S11_EvAccum",...
    "S12_EvAccum", "S13_EvAccum", "S14_EvAccum", "S15_EvAccum"];
vars = {'d'}; % which variables do you need
rootdir = 'C:\Users\doria\Google Drive\04 Research\05 Evidence Accumulation\01 EvAccum Code';
datadir = fullfile(rootdir, 'data');

for i = 1:length(datafilenames)
    load(fullfile(datadir,strcat(datafilenames(i),'.mat')),vars{:});
    rts.all(:,:,i) = d.rt;
    rts.max_all(:,i) = max(d.rt); % get max rts for each block, every participant
    rts.min_all(:,i) = min(d.rt(d.rt>0)); % same for min rts (excluding zero) - haven't figured out how to get this for every block. Is just getting one value per subject
end

rts.max_sub = max(rts.max_all); % for each subject
rts.min_sub = min(rts.min_all);
rts.max_overall = max(rts.max_sub); % overall
rts.min_overall = min(rts.min_sub);

ii = 0;
for i = 1:0.1:1.5 % for rts higher these values
    ii = ii+1;
    pc(ii,1) = i; % record which rt range
    pc(ii,2) = numel(rts.all(rts.all>i))/numel(rts.all); % calculate proportion of responses higher than the value
end
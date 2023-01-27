%% set up

close all;
clearvars;
clc;

fprintf('setting up %s\n', mfilename);
p = struct(); % est structure for parameter values
d = struct(); % est structure for data
t = struct(); % another structure for temp vars

% extract preprocessed sensor data into cosmo friendly format
% set up some directory paths
rootdir    = '/group/woolgar-lab/projects/Dorian/EvAccum/';
datadir = fullfile(rootdir,'data','meg_pilot_1');
decodingdatadir = fullfile(datadir,'decoding'); % find or make directory to output decoding results
p.filename = fullfile(decodingdatadir, 'extracted_data');

% add tools
addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7487')); %spm
addpath('/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/meg_data'); %where Ale's extract_chans function lives
addpath(genpath('/group/woolgar-lab/projects/Tijl/MD_dtb/Data_and_Analysis/Toolboxes/CoSMoMVPA-master/')) %cosmo
addpath(genpath('/group/woolgar-lab/projects/Tijl/MD_dtb/Data_and_Analysis/Analysis/Codes/')); %some of Tijl's functions in preference
addpath(genpath('/group/woolgar-lab/projects/Tijl/MD_dtb/Data_and_Analysis/Toolboxes/libsvm3.17/'));

% load file
fprintf('loading file %s\n', p.filename);
load(p.filename)

%% lets play with the data

% d.behavioural contains (columnwise):
    %  1)  cue direction (1-4)
    %  2)  dot motion direction condition (1-8)
    %  3)  coherence difficulty (1 = easy, 2 = hard)
    %  4)  matching difficulty (1 = easy, 2 = difficult)
    %  5)  unique number for each trial condition (1-64)
    %  6)  gives a number based on 9 for meg triggers
    %  7)  conditions (HcHr, HcLr, LcHr, LcLr)
    %  8)  accuracy (0, 1, or -1 for missed trials)
    %  9)  rts from behavioural data
    % 10)  rts from MEG
    % 11)  something to tag which run this sequence of data belongs to
    % 12)  we later add a row of unique numbers to tag each trial
% first, simplest thing: decode coherence direction through time. Ignore
% cues and coherence levels. This is col 2 of the behavioural data

%all_conditions = d.behavioural(:,2)'; % coherence directions

%all_conditions = d.behavioural(:,3)'; % coherence difficulty

%all_conditions = d.behavioural(:,4)'; % matching difficulty

%t.acc_idx = find(d.behavioural(:,8)>=0);
%all_conditions = d.behavioural(t.acc_idx,8)'; % participant accuracy
all_conditions = d.behavioural(:,8)';

rel_conditions = all_conditions(d.trialinfo(1,:));% filter down to just the trials that we have MEG data for (in this case, everything)
ds.sa.targets = rel_conditions';

% now asign chunks
% simplest: use each run as a chunk
ds.sa.chunks = d.trialinfo(2,:)';


%% now actually do some decoding

fprintf('starting decoding\n');

ma={};
% pick your classifier
%ma.classifier = @cosmo_classify_libsvm; % slower
ma.classifier = @cosmo_classify_lda; % faster - better for debugging

% choose partitions
ma.partitions = cosmo_nfold_partitioner(ds.sa.chunks); % put chunks in here

% normalisation
ma.normalization='zscore'; % estimate normalisation for train and test sets using the training data. 1: zscore each channel. 2: zscore each trial. default=1

ma.output = 'winner_predictions'; % gets the predicted label for each sample. to get the predicted label for each fold, use 'fold_predictions'

measure=@cosmo_crossvalidation_measure;
nbrhood=cosmo_interval_neighborhood(ds,'time','radius',0); % look in neighbourhoods of 0 radius wide timepoints
res = cosmo_searchlight(ds, nbrhood, measure, ma);

for i = 1:size(res.samples,2) % for each timepoint
    accuracy(i)=sum(numel(find(res.samples(:,i)==res.sa.targets)))/size(res.samples,1); % accuracy = sum of correct labels/total number of predictions
end

plot(accuracy)

%res=cosmo_crossvalidation_measure(ds,ma); % put in your ds, and args
% not doing something for each timepoint - neightbourhood?

% not sure how it would deal with 8 conditions

% need more detail in output - winner predictions (label for the winning

% condition of each trial across cross validation folds)

% cosmo_interval_neighbourhood = what neighbourhood to look at, (timepoint)
% -> radius of 1

% LDA is faster so use that as classifier when debugging
% normalise first because you want to normalise training and testing
% seperately so theres no information leaking across

% check when it outputs the labels, that its using all of them, because it
% might default to two - if no way to do more, then can just do all
% pairwise comparisons
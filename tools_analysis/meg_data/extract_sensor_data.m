% AW 5/8/20
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
megdatadir = fullfile(rootdir, 'data/meg_pilot_1/megdata/');
behavdatadir   = fullfile(rootdir,'data/meg_pilot_1/behavioural/');
behavfile = fullfile(behavdatadir,'S03/subj_3_MEGRTs.mat');

% add tools
addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7487')); %spm
addpath('/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/meg_data'); %where Ale's extract_chans function lives
addpath(genpath('/group/woolgar-lab/projects/Tijl/MD_dtb/Data_and_Analysis/Toolboxes/CoSMoMVPA-master/')) %cosmo
addpath(genpath('/group/woolgar-lab/projects/Tijl/MD_dtb/Data_and_Analysis/Analysis/Codes/')); %some of Tijl's functions in preference
addpath(genpath('/group/woolgar-lab/projects/Tijl/MD_dtb/Data_and_Analysis/Toolboxes/libsvm3.17/'));

% extract the sensor data
p.filename='/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/megdata/subj_3/MEEG/Preprocess/SL_subj_3.mat';
[EEG,MEGMAG,MEGPLANAR,conditions,chanlabels,badchans, trialinfo] = extract_chans_withtrialnums(p.filename);

% now we'll pull the behaviour data
    % a reminder that it looks like this:
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
load(behavfile,'MEG_RT'); % load the behavioural data
d.behavioral = double(MEG_RT); clear MEG_RT; % convert it to double and put it into a better variable
% we can index into this with the first row of trialinfo to get at D

%% lets set up cosmo
% now we have the sensor data from three sources (EEG + 2xMEG) in the
% format: channel * timepoints * trials matrices
% cosmo wants it as trials * (timepoints*channels)

% if we want to play with all three sources, we need to stack the
% channels so let's check they're sized properly
if size(EEG,3) ~= size(MEGMAG,3) || size(EEG,3) ~= size(MEGPLANAR,3)
    error('it appears your EEG and MEG data dont have the same number of trials');
elseif size(EEG,2) ~= size(MEGMAG,2) || size(EEG,2) ~= size(MEGPLANAR,2)
    error('it appears your EEG and MEG data dont have the same number of timepoints');
end

% now we'll convert them into the cosmo format
t.num_trials = size(EEG,3); % we'll use EEG to size, but could be any source
t.num_timepoints = size(EEG,2); % since they should be the same size
% stack the three sets of channels for a reshape
t.combined = []; % init this bad boy
for itrial = 1:t.num_trials
    t.combined(:,:,itrial) = [EEG(:,:,itrial);MEGMAG(:,:,itrial);MEGPLANAR(:,:,itrial)];
end
ds.samples = reshape(t.combined,t.num_trials,[]); % here we slide the trials (D3) into the rows, and we free the columns so matlab can stack D1 and D2 into them
% note - now the sources are combined, if we normalise we need to be careful not to do that over the different kinds of sensor

% now we want some feature attributes (i.e. information about what we just did)
% ds.fa.chan = for each column (feature), what channel
% ds.fa.time = for each column (feature), what timepoint

% now for the sample attributes
% ds.sa.rep = repetition of each condition per chunk (can all be ones?)
% these following two, we'll do later but here's a description
    % ds.sa.chunks = for each row (sample), what chunk -> do I want to randomly assign things to train and test set, or deliberately assign irrelevant features so that they're balanced out (equal cooccurance of conditions in train and test set)
    % ds.sa.targets = for each row (sample), what conditions

% and finally your dataset attributes
% ds.a = information about the whole dataset

%% lets play with the data
% first, simplest thing: decode coherence direction through time. Ignore
% cues and coherence levels. This is col 2 of the behavioural data
all_conditions = d.behavioural(:,2)'; % take the directions
rel_conditions = all_conditions(trialinfo(1,:));% filter down to just the trials that we have MEG data for (in this case,
% everything)
ds.sa.targets = rel_conditions;

% now asign chunks
% simplest: use each run as a chunk
ds.sa.chunks = trialinfo(2,:)';


%% now actually do some decoding
ma={};
ma.classifier = @cosmo_classify_libsvm; % this is where you pick your classifier
ma.partitions = cosmo_nfold_partitioner(ds_s);
res=cosmo_crossvalidation_measure(ds_s,ma);




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
d.behavioural = double(MEG_RT); clear MEG_RT; % convert it to double and put it into a better variable
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
% let's start by getting some information
t.num_trials = size(EEG,3); % we'll use EEG to size, but could be any source
t.num_timepoints = size(EEG,2); % since they should be the same size
t.num_chans = size(EEG(:,1,1))+size(MEGMAG(:,1,1))+size(MEGPLANAR(:,1,1)); % delivers [sum of channels, sum of sources (accidentally but could be useful)]
% stack the three sets of channels for a reshape
t.combined = []; % init this bad boy
for itrial = 1:t.num_trials
    t.combined(:,:,itrial) = [EEG(:,:,itrial);MEGMAG(:,:,itrial);MEGPLANAR(:,:,itrial)];
end
ds.samples = reshape(t.combined,t.num_trials,[]); % here we slide the trials (D3) into the rows, and we free the columns so matlab can stack D1 and D2 into them
% note - now the sources are combined, if we normalise we need to be careful not to do that over the different kinds of sensor

% now we want some feature attributes (i.e. information about what we just did)
% ds.fa.time = for each column (feature), what timepoint
% ds.fa.chan = for each column (feature), what channel
% there are also ds.fa.freq for frequencies, but Lydia hasn't mentioned these so we'll leave it
t.timecol = [];
t.chancol = [];
t.thesenames = [];
t.eeg_chan_nums = 1:length(EEG(:,1,1)); % I'm not actually sure this is right, but when we find out we can alter this
t.megmag_chan_nums = 1:length(MEGMAG(:,1,1));
t.megplanar_chan_nums = 1:length(MEGPLANAR(:,1,1));
for itimepoint = 1:t.num_timepoints % seems easiest to do this by timepoint (although slowish)
    % so we'll do the timepoint for all channels first
    t.thistime(1:t.num_chans) = itimepoint; % create a vector of itimepoint as long as the total number of channels
    t.timecol = [t.timecol,t.thistime]; % stack it
    
    % now we'll do the channels for this timepoint
    % I'm doing an index for each type here, though we could also do an
    % index over all three and use the channel names to figure them out
    t.thesechans = [t.eeg_chan_nums, t.megmag_chan_nums, t.megplanar_chan_nums]; % stack the channel numbers
    t.chancol = [t.chancol, t.thesechans]; % stack that stack!
    
    % but since I am doing an index for each typ, we're going to get some
    % information here about the type of sensor to help us
    t.eegnames(1:max(t.eeg_chan_nums)) = "EEG";
    t.megmagnames(1:max(t.megmag_chan_nums)) = "MEGMAG";
    t.megplanarnames(1:max(t.megplanar_chan_nums)) = "MEGPLANAR";
    t.thesenames = [t.thesenames,t.eegnames,t.megmagnames,t.megplanarnames];
end
ds.fa.time = t.timecol;
ds.fa.chan = t.chancol; % so this is now the index for each kind of sensor
ds.fa.chanstring = t.thesenames; % and this is the type of sensor

% now for the sample attributes
% ds.sa.rep = repetition of each condition per chunk (can all be ones?)
% these following two, we'll do later but here's a description
    % ds.sa.chunks = for each row (sample), what chunk -> do I want to randomly assign things to train and test set, or deliberately assign irrelevant features so that they're balanced out (equal cooccurance of conditions in train and test set)
    % ds.sa.targets = for each row (sample), what conditions

% and finally your dataset attributes
% ds.a = information about the whole dataset
% for MEG, this is the names of the feature attributes
% I suppose I'll just identify the sources for now?
% ds.a.fdim.values{1} = channel names
ds.a.fdim.values{1} = [string(chanlabels{1}),string(chanlabels{2}),string(chanlabels{3})]; % convert chan labels to string array in order EEG, MEGMAG, MEGPLANAR
%because right now the index values in ds.fa.chans isn't helpful
% ds.a.fdim.values{2} = timepoint names, what epoch are you using in secs or ms
ds.a.fdim.values{2} = -1500:1:2500; % our epoch is -1.5s to 2.5s
% ds.a.fdim.values{3} = frequencies names
% again, Lydia didn't mention these
% then add the labels we're using (i guess?)
ds.a.fdim.labels{1} = 'chan';
ds.a.fdim.labels{2} = 'time';


%% lets play with the data
% first, simplest thing: decode coherence direction through time. Ignore
% cues and coherence levels. This is col 2 of the behavioural data
all_conditions = d.behavioural(:,2)'; % take the directions
rel_conditions = all_conditions(trialinfo(1,:));% filter down to just the trials that we have MEG data for (in this case,
% everything)
ds.sa.targets = rel_conditions';

% now asign chunks
% simplest: use each run as a chunk
ds.sa.chunks = trialinfo(2,:)';


%% now actually do some decoding
ma={};
ma.classifier = @cosmo_classify_libsvm; % this is where you pick your classifier
ma.partitions = cosmo_nfold_partitioner(ds.sa.chunks); % put chunks in here
res=cosmo_crossvalidation_measure(ds,ma); % put in your ds, and args




% AW 5/8/20

% extract preprocessed sensor data into cosmo friendly format
rootdir    = '/group/woolgar-lab/projects/Dorian/EvAccum/';
megdatadir = fullfile(rootdir, 'data/meg_pilot_1/megdata/');
behavdatadir   = fullfile(rootdir,'data/meg_pilot_1/behavioural/');
behavfile = fullfile(behavdatadir,'S03/subj_3_MEGRTs.mat');

addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7487')); %spm
addpath('/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/meg_data'); %where Ale's extract_chans function lives
addpath(genpath('/group/woolgar-lab/projects/Tijl/MD_dtb/Data_and_Analysis/Toolboxes/CoSMoMVPA-master/')) %cosmo
addpath(genpath('/group/woolgar-lab/projects/Tijl/MD_dtb/Data_and_Analysis/Analysis/Codes/')); %some of Tijl's functions in preference
addpath(genpath('/group/woolgar-lab/projects/Tijl/MD_dtb/Data_and_Analysis/Toolboxes/libsvm3.17/'));

%filename = '/imaging/at07/Matlab/Projects/Dorian/evaccum2/data/meg_pilot_1/megdata/subj_3/MEEG/Preprocess/SL_subj_3.mat';
%filename='/imaging/aw02/tempEA/SL_subj_3.mat';
filename='/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/megdata/subj_3/MEEG/Preprocess/SL_subj_3.mat';
[EEG,MEGMAG,MEGPLANAR,conditions,chanlabels,badchans, trialinfo] = extract_chans_withtrialnums(filename);
    % re-composed matrix looks like:
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
condnums = double(MEG_RT(:,5)); % pull the condition numbers and convert them from string array to number array
% so now we have the condition numbers (condnums) and the trials are tagged
% with the value needed to index into condnums (first row of trialinfo)

% Yields channel * timepoints * trials matrices

% cosmo likes it as trials and timepoints*channels (these are in columns, each column is a MEG channel and the timepoints stack columnwise)
% this is feature attributes (for each column what time and what channel)-ds.fa.chan % ds.fa.time, then ds.sa.chunks = chunks, and finally ds.a - attributes (copy them from the cosmo website) - a reminder of your stuff

% next: reshape into cosmo format: 
% ds.samples = trials * timepoints [Q: where does channel info go?]
% ds.sa.targets = list of conditions (numeric)
% 	we can load the design matrix into a field of sa and then call cols into targets as needed (thanks Lydia)
% will also need chunks - any restrictions on what to hold out?
% ds.sa.rep can all be ones (is the repetitions of each condition per
% chunk)

% ds.samples = reshape(EEG x x x)

% Q: how best to combine over sensors?
% just add them as extra channels in the timepoints - but if we normalise we need to be careful not to do that over the different kinds of sensor
% Q: pseudotrials?
% do I want to randomly assign things to train and test set, or deliberately assign irrelevant features so that they're balanced out (equal cooccurance of conditions in train and test set) 

% first, simplest thing: decode coherence direction through time. Ignore
% cues and coherence levels. This is MEG_RT col 2.
all_conditions = double(MEG_RT(:,2))'; %take everything
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




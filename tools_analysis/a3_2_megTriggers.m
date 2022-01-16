function a3_2_megTriggers(thisSubject)

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab/';
datadir = fullfile(rootdir,'data','meg_pilot_3',thisSubject.id);
origdatadir = fullfile('/megdata/cbu/evaccum',thisSubject.meg_fld,thisSubject.date_meg);
addpath(genpath(fullfile(rootdir,'tools_analysis')))

% set up to use fieldtrip
ft_defaults;
cfg = [];

% tarfld  = fullfile(droot, 'Preprocess'); % output (target) folder
% behav   = fullfile(behavioural_data,[thisSubject.id '_EvAccum.mat']);
% trfile  = fullfile(tarfld, 'run%s_trl.mat');
% megrt   = fullfile(behavioural_data,[thisSubject.id '_MEGRTs.mat']);
% datafld = fullfile(droot, 'MaxfilterOutput');

for runi = 1:numel(thisSubject.meg_labs)

    %rawData = fullfile(datadir,[thisSubject.meg_labs{runi} '_raw.fif']);
    rawData = fullfile(origdatadir,[thisSubject.meg_runs{runi} '.fif']);
    cfg.dataset = rawData;
    thisData = ft_preprocessing(cfg);

end

return
end
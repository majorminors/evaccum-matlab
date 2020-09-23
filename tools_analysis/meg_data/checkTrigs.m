% check on Dorian's triggers

file =     '/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/megdata/subj_3/MEEG/Preprocess/Run5_subj_3_trans.mat';
D = spm_eeg_load(file);
trig_scaled = D(D.indchannel('STI101'),:)/1e6;
unique(trig_scaled)

% semms to have the right values (except for run 5) but sometimes summing
% with response buttons. Adjust code to avoid this - should fix it!
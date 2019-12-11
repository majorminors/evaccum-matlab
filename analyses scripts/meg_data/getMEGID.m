function ID = getMEGID(subject)
% retrieve scans IDs and details for the PD rdk study.
% meg .fif
% windows: \\cbsu\data\Scandata\MEG\cbu\evaccum\megid\megdate\run
% unix: /megdata/cbu/evaccum/megid/megdate/run

switch subject
    
    case {'DM_evaccumpilot_1','DM_evaccumpilot_01'}
        
        ID.meg_fld  = 'meg19_0423';
        ID.meg_runs = {'1_1_raw', '1_2_raw', '1_3_raw', '1_4_raw', '1_5_raw', '1_6_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6'};
        ID.meg_rest = ''; %?
        ID.meg_rlab = ''; %?
        ID.mri      = '';
        ID.date_mri = '';
        
        ID.bad_eeg = [];
        ID.bad_meg = [];
        ID.date_meg= '191125';
        ID.runblks = '6';
        ID.note    = 'button mapping swap after run3, matlab data S#_1';
        ID.emg     = 1;
        ID.ctr     = 0; %?
        ID.usable = 1;
        
        ID.subj  = 'DM_evaccumpilot_02';
        
    case {'DM_evaccumpilot_2','DM_evaccumpilot_02'}
        
        ID.meg_fld  = 'meg19_0429';
        ID.meg_runs = {'2_1_raw', '2_1b_raw', '2_1c_raw'};
        ID.meg_labs = {'Run1','Run2','Run3'};
        ID.meg_rest = ''; %?
        ID.meg_rlab = ''; %?
        ID.mri      = '';
        ID.date_mri = '';
        
        ID.bad_eeg = [];
        ID.bad_meg = [];
        ID.date_meg= '191128';
        ID.runblks = '4';
        ID.note    = 'run1 is 4 blocks (no matlab data), run 2 is 2 blocks (matlab data S#_1), run 3 is 4 blocks (matlab data (S#)';
        ID.emg     = 1;
        ID.ctr     = 0; %?
        ID.usable = 1;
        
        ID.subj  = 'DM_evaccumpilot_02';
        
    case {'DM_evaccumpilot_3','DM_evaccumpilot_03'}
        
        ID.meg_fld  = 'meg19_212';
        ID.meg_runs = {'3_1_raw','3_2_raw','3_3a_raw','3_4_raw','3_5_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5'};
        ID.meg_rest = ''; %?
        ID.meg_rlab = ''; %?
        ID.mri      = '';
        ID.date_mri = '';
        
        ID.bad_eeg = [];
        ID.bad_meg = [];
        ID.date_meg= '191202';
        ID.runblks = '6';
        ID.note    = 'run 3 all except last 2 blocks of 8, run 4 includes the last 8 of 3, run 5 is a single block';
        ID.emg     = 1;
        ID.ctr     = 0; %?
        ID.usable = 1;
        
        ID.subj  = 'DM_evaccumpilot_03';
end






ID.meg_folder = sprintf('/megdata/cbu/evaccum/%s/%s/',ID.meg_fld,ID.date_meg);
ID.mri_folder = sprintf('/mridata/cbu/%s_*',ID.mri);
if ~isfield(ID,'subj'); ID.subj = subject;end
%sanity check
if numel(ID.meg_runs) ~= numel(ID.meg_labs); error('Number of runs and labels does not match');end

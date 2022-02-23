function subjects = importParticipants(subject)

%%% NOTE: we want to enter ID.runid as [<blocks in run 1>;<blocks in run
%%% 2>;<blocks in run n>]

subjects = {};


%% --- S01 carl --- %%

thisSubject.id  = 'S01';
thisSubject.num  = 1;

thisSubject.meg_fld  = 'meg22_018';
thisSubject.date_meg= '220201';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw' 'run5_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4' 'Run5'};
thisSubject.runid = [1:4;4:8;9:12;13:16;17:20]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = '20';
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [1712];

thisSubject.note    = 'shifting in run3; noisy temporal channels in run5';
thisSubject.usable = 1;

thisSubject.mri      = 'CBU210277';
thisSubject.date_mri = '';
        
thisSubject.movement = [-0.5 20.7 -48.0;-1.9 19.6 -46.1;-1.6 19.4 -46.2;-0.9 18.5 -47.4;-2.0 18.0 -47.2]

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S02 charlotte --- %%

thisSubject.id  = 'S02';
thisSubject.num  = 2;

thisSubject.meg_fld  = 'meg22_032';
thisSubject.date_meg= '220215';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;4:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = '16';
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [1712];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = 'none yet';
thisSubject.date_mri = '';
        
thisSubject.movement = [-10.7 7.6 -42.5;-9.3 12 -43.5;-11.1 13.3 -40.4;-10.9 9.4 -44.1;-9.3 14.4 -42.3]

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S03 tabassum --- %%

thisSubject.id  = 'S03';
thisSubject.num  = 3;

thisSubject.meg_fld  = 'meg22_041';
thisSubject.date_meg= '220223';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;4:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = '16';
 
thisSubject.bad_eeg = [32,58];
thisSubject.bad_meg = [1712];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = 'CBU220146';
thisSubject.date_mri = '';
        
thisSubject.movement = [-6.2,1.0,-34.8;-6.2,0.9,-39.6;-7.8,-4.3,-38.2;-6.6,-6.6,-42.1]

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

% ------------------ %


return
end

%%%%%%%%%%%%%%%%%%
%% subfunctions %%
%%%%%%%%%%%%%%%%%%

function thisSubject = doExtras(thisSubject)

thisSubject.meg_folder = sprintf(...
    '/megdata/cbu/evaccum/%s/%s/',...
    thisSubject.meg_fld,...
    thisSubject.date_meg...
    );
thisSubject.mri_folder = sprintf(...
    '/mridata/cbu/%s_*',...
    thisSubject.mri...
    );

if numel(thisSubject.meg_runs) ~= numel(thisSubject.meg_labs); error('Number of runs and labels does not match');end

return
end

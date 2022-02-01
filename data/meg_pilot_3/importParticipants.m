function subjects = importParticipants(subject)

%%% NOTE: we want to enter ID.runid as [<blocks in run 1>;<blocks in run
%%% 2>;<blocks in run n>]

subjects = {};

%% S01 Runhao

thisSubject.id  = 'S01';
thisSubject.num  = 1;

thisSubject.meg_fld  = 'meg21_208';
thisSubject.date_meg= '211118';

thisSubject.meg_runs = {'run1_raw'};
thisSubject.meg_labs = {'Run1'};
thisSubject.runid = [1:3]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = '3';
thisSubject.checkTrigs = 0;
thisSubject.deleteMultiTrigs = 0;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = 'both first and third blocks might not be usable';
thisSubject.usable = 1;

thisSubject.mri      = 'CBU210531';
thisSubject.date_mri = '20211005';
        
thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% S02 Madelena
        
thisSubject.id  = 'S02';
thisSubject.num  = 2;

thisSubject.meg_fld  = 'meg21_209';
thisSubject.date_meg= '211125';

thisSubject.meg_runs = {'run1_raw','run2_raw'}; 
thisSubject.meg_labs = {'Run1', 'Run2'};
thisSubject.runid = [1:6;NaN,8:10,NaN,NaN]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = '10';
thisSubject.checkTrigs = [0,1];
thisSubject.deleteMultiTrigs = [1,0];
thisSubject.reduceTriggers = [1,2]; % if you're cutting blocks off, 1 for cut trials at end, 2 for cut trials at start
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = 'CBU210597';
thisSubject.date_mri = '20211029';

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% S03 Ashley
        
thisSubject.id  = 'S03';
thisSubject.num  = 3;

thisSubject.meg_fld  = 'meg21_213';
thisSubject.date_meg= '211216';

thisSubject.meg_runs = {'Run1_raw','Run2_raw','Run3_raw'}; % run1_raw-1.fif
thisSubject.meg_labs = {'Run1','Run2','Run3'};
thisSubject.runid = [1:6;7:10,NaN,NaN]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = '10';
thisSubject.checkTrigs = [1,1,1];
thisSubject.deleteMultiTrigs = 0;
thisSubject.reduceTriggers = [1,2,0];
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = 'CBU210597';
thisSubject.date_mri = '20211029';

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

return
end

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

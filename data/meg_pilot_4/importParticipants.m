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
thisSubject.runid = [1:4;5:8;9:12;13:16;17:20]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
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
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
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

thisSubject.meg_fld  = 'meg22_042';
thisSubject.date_meg= '220223';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
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

%% --- S04 omer --- %%

thisSubject.id  = 'S04';
thisSubject.num  = 4;

thisSubject.meg_fld  = 'meg22_049';
thisSubject.date_meg= '220228';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = '16';
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-6.5,14.4,-58.4;-4.8,14.3,-54;-5,13,-56.3;-9.6,7.7,-61.2]

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S05 joshua --- %%

thisSubject.id  = 'S05';
thisSubject.num  = 5;

thisSubject.meg_fld  = 'meg22_050';
thisSubject.date_meg= '220301';

thisSubject.meg_runs = {'run1_raw' 'run2_raw'};
thisSubject.meg_labs = {'Run1' 'Run2'};
thisSubject.runid = [1:4;5:8]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = '8';
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [4.5,10.8,-44.6;5.8,10.3,-44.7;5.3,6.9,-46.1]

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S06 joshua --- %%
thisSubject.movement = [6.2,13.3,-44.1;4.8,7.5,-46.9]
thisSubject.meg_runs = {'run3_raw' 'run4_raw' 'run5_raw'};
thisSubject.meg_labs = {'Run3' 'Run4' 'Run5'};
thisSubject.runid = [9:12;13:16;17:20]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = '1';
% ------------------ %
%% --- S07 james --- %%

thisSubject.id  = 'S07';
thisSubject.num  = 7;

thisSubject.meg_fld  = 'meg22_054';
thisSubject.date_meg= '220304';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = '16';
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [0.4,13.3,-39.9;2.4,11.4,-41.5;-1.2,12.6,-43.4;3.7,11.8,-35.5]

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S08 sandra --- %%

thisSubject.id  = 'S08';
thisSubject.num  = 8;

thisSubject.meg_fld  = 'meg22_059';
thisSubject.date_meg= '220308';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = '16';
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-4.3,8.3,-47.3;-4.7,8.1,-48.7;-4.8,7.8,-51.2;-6.9,7.7,-52.6]

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S09 joy --- %%

thisSubject.id  = 'S09';
thisSubject.num  = 9;

thisSubject.meg_fld  = 'meg22_060';
thisSubject.date_meg= '220308';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw' 'run5_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4' 'Run5'};
thisSubject.runid = [1:4;5:7;8,NaN,NaN,NaN;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = '16';
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = 'eye tracking a bit fucked because of goggles use, and also I stopped it after the first run and restarted it with title "Manual Recording Session"';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [0.2,15.9,-44.3;-0.7,15.5,-54.2;-0.3,14.1,-55.3;1.0,14.4,-55.3;-1.3,15.1,-55.7]

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S10 alex shizhang --- %%

thisSubject.id  = 'S10';
thisSubject.num  = 10;

thisSubject.meg_fld  = 'meg22_61';
thisSubject.date_meg= '220309';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = '16';
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = 'first cue and maybe one or two trials are a false start';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-0,-3.9,-46.9;-0.5,-4.3,-47.2;-1,-4.2,-49.2;-0.3,-2.2,-45.9];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S11 eva --- %%

thisSubject.id  = 'S11';
thisSubject.num  = 11;

thisSubject.meg_fld  = 'meg22_061';
thisSubject.date_meg= '220309';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = '16';
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = 'some interference from eye makeup';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-6.6,8.1,-42.8;-8.2,3.2,-44.7;-6.8,2.5,-46.1;-6.7,1.4,-46.6];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S12 Janis --- %%

thisSubject.id  = 'S12';
thisSubject.num  = 12;

thisSubject.meg_fld  = 'meg22_062';
thisSubject.date_meg= '220309';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = '16';
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [5,1,-46.5;4.8,0.7,-47.6;8,-0.2,-46.7;7,0,-46.6];


thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S13 Jess --- %%

thisSubject.id  = 'S13';
thisSubject.num  = 13;

thisSubject.meg_fld  = 'meg22_063';
thisSubject.date_meg= '220310';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = '16';
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [3.1,7.5,-38.3;3.2,7.0,-38.9;3,4.8,-41.4;0.8,3.7,-39.1];


thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S14 Duaa --- %%

thisSubject.id  = 'S14';
thisSubject.num  = 14;

thisSubject.meg_fld  = 'meg22_06';
thisSubject.date_meg= '2203';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = '16';
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-4.6,9.0,-46.3;-5,8.7,-49.8;-1.4,11.3,-55;-0.9,9.1,-50.9;-0.6,9.5,-50.8];


thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

% ------------------ %


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

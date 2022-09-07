function subjects = importParticipants(subject)

%%% NOTE: we want to enter ID.runid as [<blocks in run 1>;<blocks in run
%%% 2>;<blocks in run n>]

subjects = {};


%% --- S01 carl --- %% U01

thisSubject.id  = 'S01';
thisSubject.num  = 1;

thisSubject.meg_fld  = 'meg22_018';
thisSubject.date_meg= '220201';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw' 'run5_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4' 'Run5'};
thisSubject.runid = [1:4;5:8;9:12;13:16;17:20]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 20;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [1712];

thisSubject.note    = 'shifting in run3; noisy temporal channels in run5';
thisSubject.usable = 1;

thisSubject.mri      = 'CBU210277';
thisSubject.date_mri = '20210706';
        
thisSubject.movement = [-0.5 20.7 -48.0;-1.9 19.6 -46.1;-1.6 19.4 -46.2;-0.9 18.5 -47.4;-2.0 18.0 -47.2];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S02 charlotte --- %% U02

thisSubject.id  = 'S02';
thisSubject.num  = 2;

thisSubject.meg_fld  = 'meg22_032';
thisSubject.date_meg= '220215';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [1712];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = 'none yet';
thisSubject.date_mri = '';
        
thisSubject.movement = [-10.7 7.6 -42.5;-9.3 12 -43.5;-11.1 13.3 -40.4;-10.9 9.4 -44.1;-9.3 14.4 -42.3];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S03 tabassum --- %% U03

thisSubject.id  = 'S03';
thisSubject.num  = 3;

thisSubject.meg_fld  = 'meg22_042';
thisSubject.date_meg= '220223';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [32,58];
thisSubject.bad_meg = [1712];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = 'CBU220146';
thisSubject.date_mri = '';
        
thisSubject.movement = [-6.2,1.0,-34.8;-6.2,0.9,-39.6;-7.8,-4.3,-38.2;-6.6,-6.6,-42.1];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S04 omer --- %% U04

thisSubject.id  = 'S04';
thisSubject.num  = 4;

thisSubject.meg_fld  = 'meg22_049';
thisSubject.date_meg= '220228';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-6.5,14.4,-58.4;-4.8,14.3,-54;-5,13,-56.3;-9.6,7.7,-61.2];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S05 joshua --- %% U05

thisSubject.id  = 'S05';
thisSubject.num  = 5;

thisSubject.meg_fld  = 'meg22_050';
thisSubject.date_meg= '220301';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw' 'run5_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4' 'Run5'};
thisSubject.runid = [1:4;5:8;9:12;13:16;17:20]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = [20,8,12]; % total, first file, second file
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = 'run 3, 4, 5 in separate file';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [4.5,10.8,-44.6;5.8,10.3,-44.7;5.3,6.9,-46.1;6.2,13.3,-44.1;4.8,7.5,-46.9];

% S05_2 %
% thisSubject.id  = 'S06';
% thisSubject.num  = 6;
% thisSubject.movement = [6.2,13.3,-44.1;4.8,7.5,-46.9];
% thisSubject.meg_runs = {'run3_raw' 'run4_raw' 'run5_raw'};
% thisSubject.meg_labs = {'Run3' 'Run4' 'Run5'};
% thisSubject.runid = [9:12;13:16;17:20]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
% thisSubject.runblks = 1;


thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S07 james --- %% U06

thisSubject.id  = 'S07';
thisSubject.num  = 7;

thisSubject.meg_fld  = 'meg22_054';
thisSubject.date_meg= '220304';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [0.4,13.3,-39.9;2.4,11.4,-41.5;-1.2,12.6,-43.4;3.7,11.8,-35.5];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S08 sandra --- %% U07

thisSubject.id  = 'S08';
thisSubject.num  = 8;

thisSubject.meg_fld  = 'meg22_059';
thisSubject.date_meg= '220308';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-4.3,8.3,-47.3;-4.7,8.1,-48.7;-4.8,7.8,-51.2;-6.9,7.7,-52.6];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S09 joy --- %% U08

thisSubject.id  = 'S09';
thisSubject.num  = 9;

thisSubject.meg_fld  = 'meg22_060';
thisSubject.date_meg= '220308';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw' 'run5_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4' 'Run5'};
thisSubject.runid = [1:4;5:7,NaN;8,NaN,NaN,NaN;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = 'eye tracking a bit fucked because of goggles use, and also I stopped it after the first run and restarted it with title "Manual Recording Session"';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [0.2,15.9,-44.3;-0.7,15.5,-54.2;-0.3,14.1,-55.3;1.0,14.4,-55.3;-1.3,15.1,-55.7];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S10 alex shizhang --- %% U09

thisSubject.id  = 'S10';
thisSubject.num  = 10;

thisSubject.meg_fld  = 'meg22_61';
thisSubject.date_meg= '220309';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = 'first cue and maybe one or two trials are a false start';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-0,-3.9,-46.9;-0.5,-4.3,-47.2;-1,-4.2,-49.2;-0.3,-2.2,-45.9];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S11 eva --- %% U10

thisSubject.id  = 'S11';
thisSubject.num  = 11;

thisSubject.meg_fld  = 'meg22_061';
thisSubject.date_meg= '220309';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = 'some interference from eye makeup';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-6.6,8.1,-42.8;-8.2,3.2,-44.7;-6.8,2.5,-46.1;-6.7,1.4,-46.6];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S12 Janis --- %% U11

thisSubject.id  = 'S12';
thisSubject.num  = 12;

thisSubject.meg_fld  = 'meg22_062';
thisSubject.date_meg= '220309';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [5,1,-46.5;4.8,0.7,-47.6;8,-0.2,-46.7;7,0,-46.6];


thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S13 Jess --- %% U12

thisSubject.id  = 'S13';
thisSubject.num  = 13;

thisSubject.meg_fld  = 'meg22_063';
thisSubject.date_meg= '220310';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [3.1,7.5,-38.3;3.2,7.0,-38.9;3,4.8,-41.4;0.8,3.7,-39.1];


thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S14 Duaa --- %% U13

thisSubject.id  = 'S14';
thisSubject.num  = 14;

thisSubject.meg_fld  = 'meg22_064';
thisSubject.date_meg= '220310';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-4.6,9.0,-46.3;-5,8.7,-49.8;-1.4,11.3,-55;-0.9,9.1,-50.9;-0.6,9.5,-50.8];


thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S15 Max --- %% U14

thisSubject.id  = 'S15';
thisSubject.num  = 15;

thisSubject.meg_fld  = 'meg22_066';
thisSubject.date_meg= '220311';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw' 'run5_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4' 'Run5'};
thisSubject.runid = [1:3,NaN;4,NaN,NaN,NaN;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = [16,4,12]; % total, first file, second file
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = 'missed start of fourth block; runs 3-5 in a seperate file; lawnmowing during run4';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [0.2,15.1,-50.9;-0.5,15.0,-51.1;-1.0,15.1,-51.7;-2.1,11.5,-51.8;-1.1,10.8,-52.4;-0.6,10.9,-51.5];

% S15_2 %
% thisSubject.id  = 'S16';
% thisSubject.num  = 16;
% thisSubject.meg_runs = {'run3_raw' 'run4_raw' 'run5_raw'};
% thisSubject.meg_labs = {'Run3' 'Run4' 'Run5'};
% thisSubject.runid = [5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
% thisSubject.runblks = 12;
% thisSubject.note    = 'lawnmowing during run4';
% thisSubject.movement = [-2.1,11.5,-51.8;-1.1,10.8,-52.4;-0.6,10.9,-51.5];


thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;


%% --- S17 Jun Xun --- %% U15

thisSubject.id  = 'S17';
thisSubject.num  = 17;

thisSubject.meg_fld  = 'meg22_068';
thisSubject.date_meg= '220314';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = [16,4,12]; % total, first file, second file
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = 'run 2-4 in a seperate file';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-1.2,-5.5,-43.8;-1.2,-4.4,-61.5;-2.8,-3.5,-42.2;-4.8,1.4,-38.6;-3.6,-3.5,-40.5];

% S17_2 % 
% thisSubject.id  = 'S18';
% thisSubject.num  = 18;
% thisSubject.meg_fld  = 'meg22_068';
% thisSubject.date_meg= '220314';
% thisSubject.meg_runs = {'run2_raw' 'run3_raw' 'run4_raw'};
% thisSubject.meg_labs = {'Run2' 'Run3' 'Run4'};
% thisSubject.runid = [5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
% thisSubject.runblks = 12;
% thisSubject.movement = [-2.8,-3.5,-42.2;-4.8,1.4,-38.6;-3.6,-3.5,-40.5];


thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S19 Duncan--- %% U16

thisSubject.id  = 'S19';
thisSubject.num  = 19;

thisSubject.meg_fld  = 'meg22_071';
thisSubject.date_meg= '220315';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [1.5,15.6,-46.3;1.4,17.4,-46.2;0.7,15.9,-45.8;1.8,15.9,-48.4];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S20 Monami --- %% U17

thisSubject.id  = 'S20';
thisSubject.num  = 20;

thisSubject.meg_fld  = 'meg22_072';
thisSubject.date_meg= '220315';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-1.7,13.9,-45.1;1,13.6,-43.8;0.6,13.2,-44.8;1.1,10.7,-44.6];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S21 Arqum --- %% U18

thisSubject.id  = 'S21';
thisSubject.num  = 21;

thisSubject.meg_fld  = 'meg22_073';
thisSubject.date_meg= '220316';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-1.7,-0.2,-54.8;-3.6,-0.2,-55.2;-3.5,-0.4,-55.1;-3.6,0.2,-56.3];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S22 Jeff --- %% U19

thisSubject.id  = 'S22';
thisSubject.num  = 22;

thisSubject.meg_fld  = 'meg22_074';
thisSubject.date_meg= '220316';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-5.1,23.2,-50.6;-6.1,21,-51.8;-5.8,21.4,-51.9;-6.6,21.1,-53.5];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S23 Amelia --- %% U20

thisSubject.id  = 'S23';
thisSubject.num  = 23;

thisSubject.meg_fld  = 'meg22_077';
thisSubject.date_meg= '220318';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [0.1,17.5,-37.0;0.0,17.6,-37.1;2.4,15.9,-38.2;0.1,14.1,-38.2;-1.6,17.1,-38.4];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S24 Marco --- %%

thisSubject.id  = 'S24';
thisSubject.num  = 24;

thisSubject.meg_fld  = 'meg22_078';
thisSubject.date_meg= '220318';

thisSubject.meg_runs = {'run1_raw' 'run2_raw'};
thisSubject.meg_labs = {'Run1' 'Run2'};
thisSubject.runid = [1:4;5,NaN,NaN,NaN]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = 'abandoned at first block of run 2 because we got no good blocks';
thisSubject.usable = 0;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-5.7,6.2,-42.5;-6.0,6.8,-43.9];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S25 Meta --- %% U21

thisSubject.id  = 'S25';
thisSubject.num  = 25;

thisSubject.meg_fld  = 'meg22_079';
thisSubject.date_meg= '220321';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = 'banging from upstairs in run 4';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-0.6,14.0,-43.6;-0,12.2,-44.1;-0.8,13.6,-43.8;-1.5,14.2,-43.9];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S26 Camila --- %% U22

thisSubject.id  = 'S26';
thisSubject.num  = 26;

thisSubject.meg_fld  = 'meg22_080';
thisSubject.date_meg= '220322';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [1.3,13.0,-41.7;-1.1,10.2,-49.4;4.5,14.2,-38.6;1.7,9.7,-42.4];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

% --- S27 Lisa --- %%
% no trials - too magnetic

%% --- S28 Kate --- %% U23

thisSubject.id  = 'S28';
thisSubject.num  = 28;

thisSubject.meg_fld  = 'meg22_081';
thisSubject.date_meg= '220323';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-0.1,11.6,-48.7;0.6,12.2,-48];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S29 David --- %% U24

thisSubject.id  = 'S29';
thisSubject.num  = 29;

thisSubject.meg_fld  = 'meg22_083';
thisSubject.date_meg= '220324';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = 'no eyetracking';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-4.5,13.2,-51.3;-5.9,12.1,-51.3;-5.3,13.2,-57;-4.9,11.5,-57.2];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S30 Sukritta --- %% U25

thisSubject.id  = 'S30';
thisSubject.num  = 30;

thisSubject.meg_fld  = 'meg22_084';
thisSubject.date_meg= '220324';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw' 'run5_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4' 'Run5'};
thisSubject.runid = [1:4;5:8;9:12;13:16;17:20]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = [20,16,4]; % total, first file, second file
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = 'run 5 in a seperate file';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-3.3,8.0,-44.9;-3.9,5.4,-47.8;-4,4.7,-48.3;-2,4.0,-49.6;-2.4,6.9,-46.1];

% S30_2 %
% thisSubject.id  = 'S30';
% thisSubject.num  = 30;
% thisSubject.meg_runs = {'run5_raw'};
% thisSubject.meg_labs = {'Run5'};
% thisSubject.runid = [17:20]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
% thisSubject.runblks = 4;
% thisSubject.note    = '';
% thisSubject.movement = [-2.4,6.9,-46.1];


thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S31 Chris --- %% U26

thisSubject.id  = 'S31';
thisSubject.num  = 31;

thisSubject.meg_fld  = 'meg22_086';
thisSubject.date_meg= '220328';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [30,38];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-12.5,7.2,-49.6;-10.1,7.6,-51.5;-10.8,7.6,-53.0;-9.8,7.8,-53.4];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S34 Aaron --- %% U27

thisSubject.id  = 'S34';
thisSubject.num  = 34;

thisSubject.meg_fld  = 'meg22_088';
thisSubject.date_meg= '220329';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [1712];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-2.7,16.8,-45.8;-0.8,16.1,-47.5];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S35 Naomi --- %% U28

thisSubject.id  = 'S35';
thisSubject.num  = 35;

thisSubject.meg_fld  = 'meg22_089';
thisSubject.date_meg= '220330';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-2,7.6,-45.8;-2.8,7.8,-51.7;-2.3,12.2,-58;-2.5,10.6,-48.3;-1.9,13,-56.2;-1.4,13.1,-48.9];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S36 Alex --- %% U29

thisSubject.id  = 'S36';
thisSubject.num  = 36;

thisSubject.meg_fld  = 'meg22_090';
thisSubject.date_meg= '220331';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [0.4,6.1,-42.3;-0.3,4.3,-43.8;0.2,2.9,-44.3;0.4,2.2,-47.1];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S37 Madeleine --- %% U30

thisSubject.id  = 'S37';
thisSubject.num  = 37;

thisSubject.meg_fld  = 'meg22_091';
thisSubject.date_meg= '220401';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-9.0,10.2,-38.3;-8.5,12.4,-41.2;-8.5,10.9,-41.7;-8.6,11.8,-42.4];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S38 Joyce --- %%

thisSubject.id  = 'S38';
thisSubject.num  = 38;

thisSubject.meg_fld  = 'meg22_092';
thisSubject.date_meg= '220406';

thisSubject.meg_runs = {'run1_raw' 'run2_raw'};
thisSubject.meg_labs = {'Run1' 'Run2'};
thisSubject.runid = [1:4;5:8]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 8;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 0;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [0,0.3,-54.9;-0.5,-1.9,-56.6];

thisSubject = doExtras(thisSubject);
subjects = [subjects,thisSubject]; clear thisSubject;

%% --- S40 Alex --- %% U31

thisSubject.id  = 'S40';
thisSubject.num  = 40;

thisSubject.meg_fld  = 'meg22_093';
thisSubject.date_meg= '220406';

thisSubject.meg_runs = {'run1_raw' 'run2_raw' 'run3_raw' 'run4_raw'};
thisSubject.meg_labs = {'Run1' 'Run2' 'Run3' 'Run4'};
thisSubject.runid = [1:4;5:8;9:12;13:16]; % each row is a vector of the blocks in that run (e.g. [1:6;7:12])
thisSubject.runblks = 16;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [];

thisSubject.note    = '';
thisSubject.usable = 1;

thisSubject.mri      = '';
thisSubject.date_mri = '';
        
thisSubject.movement = [-9.7,4.5,-57.8;-10.1,6.3,-56.7;-10.6,5.2,-60.2;-11.3,3.7,-62];

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

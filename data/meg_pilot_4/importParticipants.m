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
thisSubject.checkTrigs = 0;
 
thisSubject.bad_eeg = [];
thisSubject.bad_meg = [1712];

thisSubject.note    = 'shifting in run3; noisy temporal channels in run5';
thisSubject.usable = 1;

thisSubject.mri      = 'CBU210277';
thisSubject.date_mri = '20210706';
        
thisSubject.movement = [-0.5 20.7 -48.0;-1.9 19.6 -46.1;-1.6 19.4 -46.2;-0.9 18.5 -47.4;-2.0 18.0 -47.2]

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

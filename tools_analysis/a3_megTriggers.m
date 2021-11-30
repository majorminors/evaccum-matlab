function a3_megTriggers(thisSubject)

addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/'));

addpath /hpc-software/matlab/cbu/

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab/';
addpath(genpath(fullfile(rootdir,'tools_analysis')))
droot = fullfile(rootdir,'data','meg_pilot_3',thisSubject.id);
behavioural_data = fullfile(rootdir,'data','meg_pilot_3','behavioural');

infld   = fullfile(droot, 'Preprocess');
behav   = fullfile(behavioural_data,[thisSubject.id '_EvAccum.mat']);
trfile  = fullfile(infld, 'run%s_%s_trl.mat');
megrt   = fullfile(behavioural_data,[thisSubject.id '_MEGRTs.mat']);
datafld = fullfile(droot, 'MaxfilterOutput');
tarfld  = infld;

filename    = 'Run%d_%s_trans.fif';
outfirst    = 'Run%d_%s_trans.mat';

addpath(droot);


[dname,IDnum] = getnames(droot,7);
bname         = getnames(behavioural_data,10);
bname = unique({bname{:}});

clc
clear settings
settings.overwrite = 1;




    
    
    settings.bEEG  = thisSubject.bad_eeg;
    settings.bMEG  = thisSubject.bad_meg;
    settings.behav = behav;
    settings.svbeh = megrt;
    settings.runid = thisSubject.runid;
    
    filenames = {};trfiles = {};
    
    for runi = 1:12 % different runs for each participant
        
        
        
        convertedfile= sprintf([tarfld outfirst], dname{subi}, runi, dname{subi});
        if ~exist(convertedfile,'file')
            
            fiffile = sprintf([datafld filename], dname{subi}, runi, dname{subi});
            if exist(fiffile,'file')
                %% convert fif files
                S=[];
                S.dataset       = fiffile;
                S.outfile       = convertedfile;
                S.save          = 0;
                S.reviewtrials  = 0;
                S.channels      = 'all';
                S.continuous    = 1;
                S.checkboundary = 0;
                
                % DO IT:
                D = spm_eeg_convert(S);
                
            end
            
        end
        if exist(convertedfile)
            filenames{end+1} = convertedfile;
            trfiles{end+1} = sprintf(trfile,dname{subi},num2str(runi),dname{subi});
        end
        
    end
    
    if ~exist(settings.svbeh,'file') || settings.overwrite
        fun_MEGtrgRT(filenames,settings.behav,settings.svbeh,settings,trfiles)
    end




return
end

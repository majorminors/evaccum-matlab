addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/'));

close all
clc
addpath /hpc-software/matlab/cbu/

addpath /imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Preprocessing/
addpath /imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Preprocessing/lib/

droot = '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/MEGData/';
infld= '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/MEGData/%s/MEEG/Preprocess/';
behav = '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Behav/MEGBehav/%s_MEGRT.mat';
trfile= '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/MEGData/%s/MEEG/Preprocess/run%s_%s_trl.mat';
megrt = '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Behav/%s_MEGtrg_RT.mat';
datafld = '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/MEGData/%s/MEEG/MaxfilterOutput/';
tarfld= '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/MEGData/%s/MEEG/Preprocess/';

filename    = 'Run%d_%s_trans.fif';
outfirst    = 'Run%d_%s_trans.mat';

addpath(droot);


[dname,IDnum] = getnames(droot,7);
bname         = getnames('/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Behav/MEGBehav/',10);
bname = unique({bname{:}});

clc
clear settings
settings.overwrite = 0;



for subi = 1:length(dname)
    
       ID = getMEGID(sprintf('AT_RDK2_%s',IDnum{subi}));
       
       settings.bEEG  = ID.bad_eeg;
       settings.bMEG  = ID.bad_meg;
       settings.behav = sprintf(behav,bname{subi});
       settings.svbeh = sprintf(megrt,dname{subi}); 
  
       filenames = {};trfiles = {};

        for runi = 1:12
                       
                
                
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
end


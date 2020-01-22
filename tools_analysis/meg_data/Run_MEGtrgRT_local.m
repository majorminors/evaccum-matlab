addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/'));

close all
clc
addpath /hpc-software/matlab/cbu/

addpath(genpath('/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/'))

droot   = '/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/megdata/';
infld   = fullfile(droot, '%s/MEEG/Preprocess/');
behav   = '/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/behavioural/%s/%s_EvAccum.mat';
trfile  = fullfile(infld, 'run%s_%s_trl.mat'); % going to have difficulty with subject naming conventions (Sxx vs subjx)
megrt   = '/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/behavioural/%s/%s_MEGRTs.mat';
datafld = fullfile(droot, '%s/MEEG/MaxfilterOutput/');
tarfld  = infld;

filename    = 'Run%d_%s_trans.fif';
outfirst    = 'Run%d_%s_trans.mat';

addpath(droot);


[dname,IDnum] = getnames(droot,7);
bname         = getnames('/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/behavioural/',10);
bname = unique({bname{:}});

clc
clear settings
settings.overwrite = 0;



for subi = 3 %:length(dname)
    
       ID = getMEGID(sprintf('DM_evaccumpilot_%s',IDnum{subi}));
       
       settings.bEEG  = ID.bad_eeg;
       settings.bMEG  = ID.bad_meg;
       settings.behav = sprintf(behav,bname{subi},bname{subi});
       settings.svbeh = sprintf(megrt,bname{subi},dname{subi}); 
  
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
end


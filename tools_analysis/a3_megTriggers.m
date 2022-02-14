function a3_megTriggers(thisSubject)

addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/'));

addpath /hpc-software/matlab/cbu/

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-for-alex/';
addpath(genpath(fullfile(rootdir,'tools_analysis')))
droot = fullfile(rootdir,'data','meg_pilot_3',thisSubject.id);
behavioural_data = fullfile(rootdir,'data','meg_pilot_3','behavioural');

tarfld  = fullfile(droot, 'Preprocess'); % output (target) folder
behav   = fullfile(behavioural_data,[thisSubject.id '_EvAccum.mat']);
trfile  = fullfile(tarfld, 'run%s_trl.mat');
megrt   = fullfile(behavioural_data,[thisSubject.id '_MEGRTs.mat']);
datafld = fullfile(droot, 'MaxfilterOutput');

filename    = 'Run%d_trans.fif';
outfirst    = 'Run%d_trans.mat';

addpath(droot);

clc
clear settings
settings.overwrite = 1;
    
settings.bEEG  = thisSubject.bad_eeg;
settings.bMEG  = thisSubject.bad_meg;
settings.behav = behav;
settings.svbeh = megrt;
settings.runid = thisSubject.runid;
settings.checkTrigs = thisSubject.checkTrigs;
settings.deleteMultiTrigs = thisSubject.deleteMultiTrigs;
settings.reduceTriggers = thisSubject.reduceTriggers;

filenames = {};trfiles = {};

for runi = 1:numel(thisSubject.meg_runs)
    

    
    convertedfile= sprintf([tarfld,filesep,outfirst], runi);
    if ~exist(convertedfile,'file')
        
        fiffile = sprintf([datafld,filesep,filename], runi);
        if exist(fiffile,'file')
            %% convert fif files
            S=[];
            S.dataset       = fiffile;
            S.outfile       = convertedfile;
            S.save          = 0;
            S.reviewtrials  = 0;
            S.channels      = 'all';
%             S.continuous    = 1;
            S.mode          = 'continuous';
            S.checkboundary = 1;
% spm_jobman('initcfg');
% matlabbatch{1}.spm.meeg.convert.dataset(1) = {fiffile};
% matlabbatch{1}.spm.meeg.convert.mode.continuous.readall = 1;
% matlabbatch{1}.spm.meeg.convert.channels{1}.all = 'all';
% matlabbatch{1}.spm.meeg.convert.outfile = '';
% matlabbatch{1}.spm.meeg.convert.eventpadding = 0;
% matlabbatch{1}.spm.meeg.convert.blocksize = 3276800;
% matlabbatch{1}.spm.meeg.convert.checkboundary = 1;
% matlabbatch{1}.spm.meeg.convert.saveorigheader = 0;
% matlabbatch{1}.spm.meeg.convert.inputformat = 'autodetect';
% 
% D = spm_jobman('run',matlabbatch);
            
            % DO IT:
            D = spm_eeg_convert(S);
            
        end
        
    end
    if exist(convertedfile)
        filenames{end+1} = convertedfile;
        trfiles{end+1} = sprintf(trfile,num2str(runi));
    end
    
end

if ~exist(settings.svbeh,'file') || settings.overwrite
    doMEGTriggers(filenames,settings.behav,settings.svbeh,settings,trfiles)
end



return
end

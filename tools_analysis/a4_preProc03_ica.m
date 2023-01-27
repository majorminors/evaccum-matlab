function a4_preProc03_ica(thisSubject,datadir,toolsdir,manual)
% here we will:
% 1. run an ica adapted from olafs mne script which will remove eye and heartrate artefacts


% set up paths
addpath /hpc-software/matlab/cbu/
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/');
addpath(fullfile(toolsdir,'fieldtrip-20190410')); ft_defaults;
addpath(genpath(fullfile(toolsdir,'lib')))
rootDir = fullfile(datadir,thisSubject.id);
preProcDir = fullfile(rootDir,'Preprocess');
addpath(rootDir);

% for this to work must set conda activate mne0.18
% in command window before opening matlab - later versions some of the functions deprecated (e.g. ica.plot_sources(exclude))%%
cd /imaging/local/software/mne_python/Utilities
addpath /imaging/local/software/mne_python/

% file definition
trlfile  = error('you havent defined this yet for all runs together');
if exist(fullfile(prexProcDir,'Cf_transmid.fif'),'file')
    preIcaFile = fullfile(preProcDir,'Cf_transmid.fif');
else
    error('no file I know to check for')
end
inputFile = preIcaFile; % needs a .fif file for the ica
[tmpPath, tmpName, tmpExt] = fileparts(preIcaFile);
outputFile = fullfile(tmpPath,[tmpName '_ica.fif']);
htmlOutput = fullfile(tmpPath,[tmpName '_ica.html']);

eog={'EOG001','EOG002'}; % eog channels
ecg = {'ECG003'}; % ecg channel

%% do ica

fprintf('identifying heart and blink artefacts with ICA\n');
icaCompCmd=sprintf(' python Fiff_Compute_ICA.py --FileRaw %s --FileICA %s --FileHTML %s --EOG %s %s --ECG %s --n_pca_comps %d',inputFile,outputFile,htmlOutput,eog{:},ecg{:},.99);
if exist('gradthresh','var') && exist('magthresh','var')
    icaCompCmd=sprintf(' python Fiff_Compute_ICA.py --FileRaw %s --FileICA %s --FileHTML %s --EOG %s %s --ECG %s --n_pca_comps %d --RejGrad %3.2f --RejMag %3.2f',inputFile,outputFile,htmlOutput,eog,ecg,.99,gradthresh,magthresh);
end

fprintf('\n\n%s\n\n',icaCompCmd);
system(icaCompCmd);

%apply ICA to rm artefacts
if rmArt
    fprintf('removing heart and blink artefacts\n');
    
    %   automatically remove components
    icaRejCmd=sprintf(' python Fiff_Apply_ICA.py --FileRawIn %s --FileICA %s',inputFile,outputFile);
    
    %	manually enter the components you want to remove
    if manual
        artefacts=input('Please enter the number of any components you wish to remove e.g. [1 4 5], then press return: ');
        icaRejCmd=sprintf(' python Fiff_Apply_ICA.py --FileRawIn %s --FileICA %s --ICAcomps %d %d %d',inputFile,outputFile,artefacts);
    end
    
    fprintf('\n\n%s\n\n',icaRejCmd);
    system(icaRejCmd);
end


end

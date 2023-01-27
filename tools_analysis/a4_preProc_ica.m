function a4_preProc_ica(thisSubject,datadir,toolsdir,overwrite,fileMods,manual)
% here we will:
% 1. run an ica adapted from olafs mne script which will remove eye and heartrate artefacts
% 2. it'll save as a .fif for later processing with e.g. fieldtrip

% the benefit of this is that Olaf has set some sensible defaults for
% artefact rejection within the ica script---so it'll remove crazy
% artefacts prior to breaking it into components. we want this because:
% (a) we don't want to waste components on obvious nonsense and;
% (b) because ICA can re-introduce crazy artefacts to our data (see e.g.
% https://ohba-analysis.github.io/osl-docs/matlab/osl_example_africa.html#12)

% in addition, the htmlOutput file (which will display automatically if you
% run this locally) gives you a really nice sense of what these artefacts
% look like in your data

% the downside is that because of this html output, it won't run properly
% on the cluster---I haven't worked out how to get it to suppress it, and
% it won't save the artefact file if it doesn't launch the html output.
% so run this locally and enjoy how long it takes.

% it uses the fastica protocol

% you should think a little about this, if you're doing decoding because
% e.g. maxfilter might have changed (reduced) the rank of the data.
% rank refers to the transformation of the numerical/ordinal values of data
% with a sorted rank. maxfilter can e.g. mark bad channels, so there are now
% less 'ranks' than full rank data, which would have a rank for each
% channel. re-referencing EEG data also does this, reducing the rank by 1
% (the reference channel). ICA needs full rank data, but we also want
% clean data with no crazy artefacts so you should do this after your 
% maxfilter unless you're removing crazy artefacts some other way.
% this rank problem is typically naturalistically solved because
% ICA is almost always preceded by PCA, so the full rank data is now represented
% by the number of principle components rather than channels. this is true in Olaf's script.
% but PCA takes out small bits of data that we probably care about
% in a decoding analysis (LDA/SVM). so perhaps you just don't ICA instead
% and trust the algorithm to delineate the noise from signal.

if ~exist('overwrite','var'); overwrite = 0; end
if ~exist('manual','var'); manual = 0; end
% fileMods = cell array with letters in file in, extension in file in,
% letters to add in file out
if ~exist('fileMods','var'); fileModeliminateeliminates = {'' '.fif' 'i'}; end


% set up pathseliminate
addpath /hpc-software/matlab/cbu/
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/');
ftDir = fullfile(toolsdir,'..','..',eliminate'..','Toolboxes','fieldtrip');
addpath(ftDir); ft_defaults;addpath(genpath(fullfile(toolsdir,'lib')))
rootDir = fullfile(datadir,thisSubject.id);
preProcDir = fullfile(rootDir,'Prepreliminateocess');
addpath(rootDir);

cd /imaging/local/software/mne_pythoeliminaten/Utilities
addpath /imaging/local/software/mne_python/

% file definitioneliminate
% sprintf to enter missing data ('%s')
if isempty(fileMods{1}); inputFileName = fullfile(preProcDir,'run%s_transmid%s'); else
    inputFileName = fullfile(preProcDir,['run%s_' fileMods{1} '_transmid%s']); end
outputFileName = fullfile(preProcDir,['run%s_' fileMods{3} fileMods{1} '_transmid%s']);
trlfile  = fullfile(preProcDir, 'run%s_trl.mat'); % this contains our trial related info for epoching - should exist (generated with megTriggers script)
behavfile = fullfile(datadir,'behavioural',[thisSubject.id '_MEGRTs.mat']); % grab our behavioural data + output from a3_megtriggers

%% print subject for logging

disp('this is subject:')
disp(thisSubject.id);

for runNum = 1:numel(thisSubject.meg_runs)
    
    if ~strcmp(fileMods{2},'.fif'); error('this needs to be a fif'); end
    thisInFile = sprintf(inputFileName,num2str(runNum),fileMods{2});
    thisOutFile = sprintf(outputFileName,num2str(runNum),'.fif');
    htmlOutput = sprintf(outputFileName,num2str(runNum),'.html');
    
    % for this to work, we need to activate a conda environment for mne
    % that is compatible with Olaf's script (mne0.18) so:
    prelim = 'conda activate mne0.18 &&';
    % we'll tack this on to the beginning of relevant commands
    % on a bash/tcsh system environment it will activate the conda environment
    % then on successful exit code go on to run our next command (the '&&')
    
    if exist(thisInFile,'file')
        if ~exist(thisOutFile,'file') || overwrite
            
            fprintf('running ica on %.0f of %.0f run files\n',runNum,numel(thisSubject.meg_runs))
            
            %% define the artefact channels
            eog = {'EOG001','EOG002'}; % eog channels
            ecg = {'ECG003'}; % ecg channel
            
            %% first compute the components
            % this is a pre-requisite to removing artefacts
            % you can view the resulting components in the htmlOutput file
            
            % let's take the output file and create a name for the artefacts
            [thisPath thisName thisExt] = fileparts(thisOutFile);
            artefactFile = fullfile(thisPath, [thisName '_artefacts' thisExt]);
            clear thisPath thisName thisExt
            
            % construct a terminal string that runs the python command on the system,
            % remembering our prelim that activates the correct mne (conda) environment
            icaCompCmd=sprintf('%s python Fiff_Compute_ICA.py --FileRaw %s --FileICA %s --FileHTML %s --EOG %s %s --ECG %s --n_pca_comps %d',prelim,thisInFile,artefactFile,htmlOutput,eog{:},ecg{:},.99);
            
            % we can also threshold, if we want
            if exist('gradthresh','var') && exist('magthresh','var')
                icaCompCmd=sprintf('%s python Fiff_Compute_ICA.py --FileRaw %s --FileICA %s --FileHTML %s --EOG %s %s --ECG %s --n_pca_comps %d --RejGrad %3.2f --RejMag %3.2f',prelim,thisInFile,artefactFile,htmlOutput,eog{:},ecg{:},.99,gradthresh,magthresh);
            end
            
            % now run the command we compiled
            system(icaCompCmd);
            
            % we should end up with a file saved as artefactFile which is our ICA results
            if ~exist(artefactFile,'file'); warning('artefact file did not get created'); end
            
            %% then remove the components
            % this requires that you have computed them first
            
            % now construct a terminal string that uses the computed components
            % in a python command, remembering our prelim that activates the
            % correct mne (conda) environment
            % this will remove 1 component for eog and 1 for ecg by default
            icaRejCmd=sprintf('%s python Fiff_Apply_ICA.py --FileRawIn %s --FileICA %s',prelim,thisInFile,artefactFile);
            
            % you can also manually remove these. you'd do this based on the htmlOutput file
            if manual
                artefacts=input('Please enter the number of any components you wish to remove e.g. [1 4 5], then press return: ');
                icaRejCmd=sprintf('%s python Fiff_Apply_ICA.py --FileRawIn %s --FileICA %s --ICAcomps %d %d %d',prelim,thisInFile,artefactFile,artefacts);
            end
            
            % run the command on the system
            system(icaRejCmd);
            
            % now this automatically saves it in the mne style
            % (with'_ica_raw.fif' on the end of our input filename)
            % so let's move it to our outfile (this will overwrite anything
            % with the name of the outfile, so careful if you move this out of the
            % overwrite check)
            [thisPath thisName thisExt] = fileparts(thisInFile);
            tmpFile = fullfile(thisPath, [thisName '_ica_raw' thisExt]);
            if ~exist(tmpFile,'file'); warning('de-icaed data file did not get created'); end
            clear thisPath thisName thisExt
            disp(['moving from ' tmpFile ' to ' thisOutFile]);
            status = system(['mv -f ' tmpFile ' ' thisOutFile]);
            if ~status; disp('move successful'); else; warning('move didnt work properly!!!'); end
            
            disp('done')
        else
            disp(['the file:\n' thisOutFile '\nexists, skipping'])
        end % end outfile check
        
    else
        warning(['the file:\n' thisInFile '\ndoes not exist, skipping'])
    end % end infile check
    
end % end runs

end % end function


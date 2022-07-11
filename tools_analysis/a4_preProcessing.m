function a4_preProcessing(thisSubject,datadir,toolsdir)

error('this script is an old attempt to convert the pipeline of A. Tomassini, but no longer works because Riks ICA expects artefact topographies and I never worked out how to get around that---as such this script should be considered unfinished')

addpath /hpc-software/matlab/cbu/
addpath(genpath('/neuro/meg_pd_1.2/'));       % FIFACCESS toolbox if use some meg_misc functions below
addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/'));
addpath(genpath(fullfile(toolsdir,'eeglab13_5_4b')))
addpath(genpath(fullfile(toolsdir,'meg_data')))


droot = fullfile(datadir,thisSubject.id);
infld   = fullfile(droot,'MaxfilterOutput');
tarfld  = fullfile(droot,'Preprocess');
behavDir   = fullfile(datadir,'behavioural');
behavfile = fullfile(behavDir,[thisSubject.id '_Evaccum.mat']);
megrtfile = fullfile(behavDir,[thisSubject.id '_MEGRTs.mat']);
icaSaveFile     = fullfile(droot,['ICAOutput/ICA_' thisSubject.id '.mat']);
outpt   = fullfile(droot, 'Preprocess/MSL_%s.mat'); % this I think it what is saved for overwriting
addpath(droot);


settings.fld = tarfld;

settings.freshstart = 1;%delete everything before starting again -not needed for 1st preprocessing run
settings.overwrite = 1;
settings.ICAoverwrite = 1;

settings.outfname = [tarfld,filesep,thisSubject.id '_trans.mat'];
%settings.ctr   = ID.ctr;
settings.bEEG  = thisSubject.bad_eeg;
settings.bMEG  = thisSubject.bad_meg;

filenames = {}; %number of files
outfirst= {};%converted files

for runNum = 1:numel(thisSubject.meg_runs)
    filenames{end+1} = sprintf([infld,filesep,thisSubject.meg_labs{runNum},'_trans.fif']);
    outfirst{end+1} = sprintf([infld,filesep,thisSubject.meg_labs{runNum},'_trans.mat']);
end

cd(tarfld);

%trial files needed for epoching
trls = getnames(tarfld,[],'run*_trl.mat');
trls = {trls{:,1}};

%% prepare data
%extract behavioural data
%load(behavfile)

if settings.freshstart
    delete *fRun*_trans.*;
    delete *SL_AT_RDK2_*; % change this to what you named it
end

for runi = 1:numel(filenames)
    
    if ~exist(outfirst{runi},'file')
        %% convert fif files
        S=[];
        S.dataset       = filenames{runi};
        S.outfile       = outfirst{runi};
        S.save          = 0;
        S.reviewtrials  = 0;
        S.channels      = 'all';
        S.continuous    = 1;
        S.checkboundary = 0;
        
        % DO IT:
        spm_eeg_convert(S);
        
    end
end

% behavioural data from triggers
if ~exist(megrtfile,'file')
    error('Prepare behavioural files first!')
end


%% Line and highpass filter
% electric line 50Hz in the uk needs to be removed
for runi = 1:length(outfirst)
    [pathstr,name,ext] = fileparts(outfirst{runi});
    outf = sprintf('%s/ff%s%s',pathstr,name,ext);
    
    if ~exist(outf,'file') || settings.overwrite
        
        
        %highpass
        S   = [];
        S.D = outfirst{runi};
        S.band = 'high';
        S.freq = 0.5; % removes the drift
        D = spm_eeg_filter(S);
        
        %line filter
        S   = [];
        S.D = D;
        S.band = 'stop';
        S.freq = [49 51];% This defines the notch filter frequency range, i.e. around 50Hz
        D = spm_eeg_filter(S);
        
        
        filenames{runi} = D.fullfile;
    else
        filenames{runi}=outf;
    end
end






%% Re-reference
% re-referencing to the average
for runi = 1:length(filenames)
    
    [pathstr,name,ext] = fileparts(filenames{runi});
    outf=sprintf('%s/M%s%s',pathstr,name,ext);
    
    
    
    if ~exist(outf,'file') || settings.overwrite
        setBadCh(settings,filenames{runi});
        
        D = reref(filenames{runi});
        filenames{runi} = D{2}.Dfname;
        
    else
        filenames{runi}=outf;
    end
end

for i = 1:numel(filenames);filenames{i}=char(filenames{i});end

%% Do ICA

clear D
D = spm_eeg_load(filenames{1});

%setup parameters
ICA.PCA_dim = 60;                       % Number of PCs for ICA
ICA.PCA_EEGdim = 60;
ICA.FiltPars = [0.1 40];                % filter bandpass [1 40] %% [0.05 20];
ICA.FiltPars = [ICA.FiltPars D.fsample];
ICA.TemAbsPval = .05/ICA.PCA_dim;       %0.05;  % Too liberal if not permuted?
ICA.SpaAbsPval = .05/ICA.PCA_dim;       %0.05;  % Too liberal if not permuted?
ICA.TemRelZval = 3;                     % Relative temporal threshold in Z-values
ICA.SpaRelZval = 2;                     % 3 is too strict? (since high topo correlation anyway)
ICA.Nperm = 1600;%1600;                       %1600 for 95%CI of p=0.01; 0 turns off permutation
ICA.VarThr = 0;                         % variance threshold (could be 100/PCA_dim, but suggest 0);
ICA.Rseed = 1;                          % to make reproducible
modalities = {'MEGMAG', 'MEGPLANAR', 'EEG'};


%path for templates
% compares ICA with my lines of interest (e.g. ECG HEOG VEOG etc)
% but won't work for EEG, because EEG artefact templates only exist for 70
% channel EEG cap
% arttopos = load('/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/meg_data/templates/MEGEEGArtifactTemplateTopographies');
s_ref_chans = {'EOG001','EOG002','ECG003'};%check D.chanlabels to determine what these are
%if ~settings.ctr
%    t_ref_chans = {'EOG061','EOG062','ECG063','EMG064'};%needed as we don't have tremor templates
%    ICA.FiltPars = repmat(ICA.FiltPars,[3,1]);
%else
t_ref_chans = s_ref_chans;
ICA.FiltPars = repmat(ICA.FiltPars,[3,1]);
%end

% takes spatial and temporal channels for comparison

ICA.s_ref_labs = s_ref_chans;
ICA.t_ref_labs = t_ref_chans;

% scrambles  (permutes) signal and compares lines to determine correlation
% against null distribution 
% (95th percentile) for definition as artifact


for iir = 1:numel(t_ref_chans); ICA.refs.tem{iir} = [];end % Reference signal for correlating with ICs
for iif = 1:numel(filenames)
    D = spm_eeg_load(filenames{iif});
    for a = 1:length(t_ref_chans)
        ICA.refs.tem{a} = [ICA.refs.tem{a},D(find(strcmp(D.chanlabels,t_ref_chans{a})),:)];%#ok %epoched data need to be reshaped into a 2-d matrix
    end
end

chans = {}; remove = {}; weights = {}; temcor = {}; spacor = {}; TraMat = {};sphere =[];%#ok

toDoW = 1:length(modalities);%
if exist(icaSaveFile,'file') && ~settings.overwrite
    
    load(icaSaveFile);
    
    toDoW(find(~cellfun(@isempty,all_ica.weights)))=[];%#ok
    %needed to decide whether to run ADJUST
    
    
end


if ~isempty(toDoW) || settings.ICAoverwrite
    
    for m = toDoW %startm:length(modalities)
        
        for refIdx = 1:size(ICA.refs.tem,2)
            ICA.refs.spa{refIdx} = cell(size(ICA.refs.tem{refIdx}));
        end
%         ICA.refs.spa = {arttopos.HEOG{m}', arttopos.VEOG{m}', arttopos.ECG{m}'}; % Assumes modalities ordered same way!!!
        
        ICA.d = [];
        chans = find(strcmp(D.chantype,modalities{m}));
        bad_ch = find(ismember(chans,D.badchannels));
        
        %remove badchannels from d and spat templates

        
        for iif = 1:numel(filenames)
            D = spm_eeg_load(filenames{iif});
            
            % D = spm_eeg_load( sprintf('_trdef_cRun%d_JS_trans',iif));
            ICA.d = [ICA.d,D(chans,:)];%All runs in one matrix
            
        end
        
        
        
        
        %run Rik's ICA with MEG templating
        [remove,weights,TraMat,temcor,spacor,varexpl,ICs,ICsEMG,sphere] = detect_ICA_artefacts_v2_PD(ICA);
        
        all_ica.temcor{m} = temcor;
        all_ica.spacor{m} = spacor;
        all_ica.remove{m} = remove;
        all_ica.TraMat{m} = TraMat;
        all_ica.varexpl{m} = varexpl;
        
        
        
        all_ica.weights{m} = weights;
        all_ica.ICs{m} = ICs(1:5000,:);
        %if ~settings.ctr
        %    all_ica.ICsEMG{m} = ICsEMG(1:5000,:);
        %end
        all_ica.sphere{m} = sphere;
        all_ica.chans{m} = chans;
        
    end
    
    saveicastuff = sprintf(' save %s  all_ica', icaSaveFile);
    
    eval(saveicastuff)
    
end




%% Do preprocessing: downsampling & epoching + label sorting
%Stimulus locked & Response locked

load(megrtfile,'megBehaviouralData'); % load behavioural file from MRGtrgRT script
numtrials = megBehaviouralData(:,end-1:end); % pull the last two columns of megBehaviouralData: what run the trial belonged to, and a unique number for the trial respectively
numtrials = str2double(numtrials); % convert to double
Xfiles=preprocPipeline(filenames,trls,numtrials);


%% Run ADJUST to detect EEG artefacts

% if length(all_ica.weights)==3 && (length(all_ica.remove)<3 || settings.ICAoverwrite)
%     
%     disp('running ADJUST');
%     clear D; D = spm_eeg_load(Xfiles);
%     
%     EEG = pop_fileio(Xfiles,'channels',all_ica.chans{3});%at
%     
%     chLab = D.sensors('EEG').label;
%     goodch = find(~ismember(1:70,D.badchannels));
%     chLab = D.chanlabels(goodch)';
%     chPos = D.sensors('EEG').chanpos;
%     chPos = D.sensors('EEG').chanpos;
%     chPos = chPos(goodch,:);
%     eeglocs = [num2cell(1:length(chLab))',num2cell(chPos),num2cell(chLab)];
%     dlmcell('eeglocations.xyz',eeglocs)
%     EEG = pop_chanedit(EEG, 'load',{'eeglocations.xyz' 'filetype' 'xyz'});
%     EEG.icasphere = all_ica.sphere{3};
%     EEG.icaweights = all_ica.weights{3};
%     EEG.icaact = []; % EEG.icaweights*EEG.icasphere*EEG.data;
%     EEG.setname = 'EEG';
%     remove = interface_ADJ_noplot(EEG,[subj,'_report']);
%     
%     
%     finalics  = setdiff(1:ICA.PCA_dim,remove);
%     iweights = pinv(EEG.icaweights);
%     TraMat    = iweights(:,finalics) * EEG.icaweights(finalics,:);
%     clear EEG;
%     if length(remove)>5
%         remove = remove(1:5);%remove only the 3 more influential components
%     end
%     all_ica.remove{3} = unique([all_ica.remove{3} remove]);
%     all_ica.TraMat{3} = TraMat;
%     
%     saveicastuff = sprintf(' save %s  all_ica TraMat', savename);
%     eval(saveicastuff)
% end

clear all_ica;

%% Apply ICA art-removal

disp('removing bad ICs for SL');
outputf{1} = sprintf('%s/SL_%s.mat',pathstr,thisSubject.id);
at_apply_ICA(Xfiles,outputf{1},icaSaveFile)



%% Re-reference
Xfiles = outputf{1};clear outf;


[pathstr,name,ext] = fileparts(Xfiles);
outf=sprintf('%s/M%s%s',pathstr,name,ext);

%if ~exist(outf,'file') || settings.overwrite

D = reref(Xfiles);



Xfiles = D{2}.Dfname{1};

%else
%   Xfiles={};
%end


%% clean up
delete *fRun*_trans.*;





%% ancillary functions
function  setBadCh(settings,fname)


D = spm_eeg_load(fname);

%remove old badchannels
D=D.badchannels(D.badchannels,0);
%identify eeg badchannels
badch = [];
if ~isempty(settings.bEEG)
    for ib = 1:numel(settings.bEEG)
        if settings.bEEG(ib)<10
            badch(ib) = D.indchannel(['EEG00',num2str(settings.bEEG(ib))]);%#ok
        else
            badch(ib) = D.indchannel(['EEG0',num2str(settings.bEEG(ib))]);%#ok
        end
    end
end

if ~isempty(settings.bMEG)
    for ib = 1:numel(settings.bMEG)
        if settings.bMEG(ib)<1000
            badch(end+1) = D.indchannel(['MEG0',num2str(settings.bMEG(ib))]);%#ok
        else
            badch(end+1) = D.indchannel(['MEG',num2str(settings.bMEG(ib))]);%#ok
        end
    end
end
%D = units(D, D.indchantype('EEG'), 'uV');
%update with new badchannels
D=D.badchannels(badch,1);
D.save;

function D = reref(fname)

matlabbatch{1}.spm.meeg.preproc.prepare.D = {fname};
matlabbatch{1}.spm.meeg.preproc.prepare.task{1}.avref.fname = 'avref_montage.mat';
matlabbatch{2}.spm.meeg.preproc.montage.D = {fname};
matlabbatch{2}.spm.meeg.preproc.montage.mode.write.montspec.montage.montagefile(1) = cfg_dep('Prepare: Average reference montage', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','avrefname'));
matlabbatch{2}.spm.meeg.preproc.montage.mode.write.montspec.montage.keepothers = 1;
matlabbatch{2}.spm.meeg.preproc.montage.mode.write.blocksize = 655360;
matlabbatch{2}.spm.meeg.preproc.montage.mode.write.prefix = 'M';
D = spm_jobman('run',matlabbatch);
spm_unlink('avref_montage.mat')



function Xfiles=preprocPipeline(infname,trls,numtrials)

spm_jobman('initcfg');


%% epoch + merge
lab = 'SL';


for iff = 1:numel(infname)
    
    clear matlabbatch;
    
    %Epoch+Rearrange labels
    matlabbatch{1}.spm.meeg.preproc.epoch.D(1) = cellstr(infname{iff});
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.trlfile = cellstr(trls{iff});
    
    matlabbatch{1}.spm.meeg.preproc.epoch.bc = 0;
    matlabbatch{1}.spm.meeg.preproc.epoch.eventpadding = 0;
    matlabbatch{1}.spm.meeg.preproc.epoch.prefix = lab;
    
    output = spm_jobman('run',matlabbatch);
    
    %% Add  info needed for MVPA
    D = spm_eeg_load(output{1}.Dfname{1});
    %info to store: [trial num, block num, here you can add anything you wish]
    thesetrials = numtrials(find(numtrials(:,1) == iff),2); % find the trial numbers related to this 'iff' which is the index for each file - files are the Maxfilter outputs ...trans.fif for each run - so what trials for this run
    info = num2cell([thesetrials';iff.*ones(1,numel(thesetrials))],1); % put that into 'info'
    D = trialtag(D, ':', info) ;save(D); % tack 'info' onto the trial data
    
    clear matlabbatch;
    matlabbatch{1}.spm.meeg.preproc.prepare.D(1) = cellstr(output{1}.Dfname);
    
    %matlabbatch{2}.spm.meeg.preproc.prepare.D(1) = cfg_dep('Epoching: Epoched Datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
    matlabbatch{1}.spm.meeg.preproc.prepare.task{1}.sortconditions.label = {
        'LcLr'
        'HcLr'
        'LcHr'
        'HcHr'
        }';
    
    matlabbatch{2}.spm.meeg.preproc.downsample.D(1) = cfg_dep('Prepare: Prepared Datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
    matlabbatch{2}.spm.meeg.preproc.downsample.fsample_new = 500;
    matlabbatch{2}.spm.meeg.preproc.downsample.method = 'resample';
    matlabbatch{2}.spm.meeg.preproc.downsample.prefix = 'd';
    
    
    output = spm_jobman('run',matlabbatch);
    
    
    
    
    %get ready for another round
    clear matlabbatch;
    outfiles{iff} = output{1}.Dfname; %#ok
    
    %check that number of trials corresponds
    D = spm_eeg_load(outfiles{iff}{1});
    trl = load(trls{iff},'trl');
    %if numel(D.conditions)~= size(trl.trl,1); error([D.fname,' discrepant number of trials!']);end
    clear trl D;
    
    
    
end


%% merge
S=[];

%S.D needs to be  a char array
S.D = char([outfiles{:}]);

S.recode = 'same';
S.prefix = 'X';
D = spm_eeg_merge(S);



Xfiles = D.fullfile;



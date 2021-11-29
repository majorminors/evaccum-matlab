function a1_importAndOrganiseScans(thisSubject)

% set data directory
rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';
datadir = fullfile(rootdir,'data','meg_pilot_3');
overwrite = 0; % turn this on for auto overwrite

% using ale's language to avoid errors in changing the names
base = fullfile(datadir,num2str(thisSubject.subj));

% if it exists choose whether to overwrite
fprintf('creating folder %s: \n',base);
if exist(base,'dir') && ~overwrite
    t.prompt = 'this folder already exists, overwrite? ([y]/n) ';
    t.ok = input(t.prompt,'s');
    if isempty(t.ok); t.ok = 'y';
    elseif ~strcmp(t.ok,'y'); error('aborted by user, not overwriting');
    elseif strcmp(t.ok,'y'); rmdir(base,'s');
    end
end

if ~exist(base,'dir'); mkdir(base); end
cd(base);

%create my folders
cellfun(@(x) mkdir(base,x),{'Preprocess','T1','MaxfilterOutput','ICAOutput','ScalpTimeStats'});

%% import MEG data
%locate raw files folder
source      = thisSubject.meg_folder;

%locate and copy raw files into my directories
sourcefiles = cellfun(@(x) sprintf('%s.fif',x),thisSubject.meg_runs, 'UniformOutput',false);
runlabs     = thisSubject.meg_labs;

runfiles    = cellfun(@(x) sprintf([base,'%s_raw.fif'],x), runlabs,'UniformOutput',false);

cd(source);
cellfun(@copyfile,sourcefiles,runfiles,'UniformOutput',false);

%% import MRI structurals
if ~isempty(thisSubject.mri)
    % locate the right destination folder
    T1Folder = [base,'T1/'];
    
    % locate MRI images
    mriFolder = dir(thisSubject.mri_folder); mriFolder = ['/mridata/cbu/',mriFolder.name];
    
    
    
    mriSubFld = dir(sprintf([mriFolder,'/%s_*'],thisSubject.date_mri));
    
    %double check we are picking the right MRI file
    if isempty(mriSubFld); error(sprintf('mismatch with MRI dates %s',base));end
    
    mriSubFld = mriSubFld.name;
    target = dir([mriFolder,'/',mriSubFld,'/*MPRAGE_*chn']);
    
    %in case MRI scan had to be repeated (e.g. if the volunteer moved)
    %choose the last one
    if size({target.name},2)>1; name = target(size({target.name},2)).name;clear target; target.name = name; end
    
    target = [mriFolder,'/',mriSubFld,'/',target.name];
    
    
    copyfile(target,T1Folder);
    
end

cd(base); save([thisSubject.subj,'_filesID.mat'],'runfiles','sourcefiles','thisSubject');

return
end
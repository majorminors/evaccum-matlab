function a1_importAndOrganiseScans(thisSubject,datadir,overwrite)

if ~exist('overwrite','var'); overwrite = 0; end

% using ale's language to avoid errors in changing the names
base = fullfile(datadir,thisSubject.id);

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
cellfun(@(x) mkdir(base,x),{'Preprocess','T1','MaxfilterOutput'});

%% import MEG data
disp('importing meg data')
%locate raw files folder
source      = thisSubject.meg_folder;

%locate and copy raw files into my directories
sourcefiles = cellfun(@(x) sprintf('%s.fif',x),thisSubject.meg_runs, 'UniformOutput',false);
runlabs     = thisSubject.meg_labs;

runfiles    = cellfun(@(x) sprintf([base,filesep,'%s_raw.fif'],x), runlabs,'UniformOutput',false);

cd(source);
cellfun(@copyfile,sourcefiles,runfiles,'UniformOutput',false);
disp('done importing meg data')

%% import MRI structurals
disp('importing mri struct(s)')
if isempty(thisSubject.mri)
    
    disp('we have not found mri structurals')
    
else
    
    disp('we have found mri structurals')
    
    % locate the right destination folder
    T1Folder = [base,filesep,'T1']; % we might need a filesep at the end
    
    % locate MRI images
    mriFolder = dir(thisSubject.mri_folder); mriFolder = ['/mridata/cbu/',mriFolder.name];
    
    
    
    mriSubFld = dir(sprintf([mriFolder,filesep,'%s_*'],thisSubject.date_mri));
    
    %double check we are picking the right MRI file
    if isempty(mriSubFld); error(sprintf('mismatch with MRI dates %s',base));end
    
    mriSubFld = mriSubFld.name;
    target = dir([mriFolder,filesep,mriSubFld,filesep,'*MPRAGE_*chn']);
    
    %in case MRI scan had to be repeated (e.g. if the volunteer moved)
    %choose the last one
    if size({target.name},2)>1; name = target(size({target.name},2)).name;clear target; target.name = name; end
    
    target = [mriFolder,filesep,mriSubFld,filesep,target.name];
    
    
    copyfile(target,T1Folder);
    
end
disp('done importing mri struct(s)')

cd(base); save([thisSubject.id,'_filesID.mat'],'runfiles','sourcefiles','thisSubject');

disp('done')

return
end

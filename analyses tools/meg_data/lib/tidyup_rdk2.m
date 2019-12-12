function [runfiles]=tidyup_rdk2(base,ID,overwrite)

    

    if ~exist(base,'dir') || overwrite
        display(sprintf('creating folder %s: ',base));
        if exist(base,'dir')
           rmdir(base,'s');
        end
        mkdir(base);
        
    end
    
    cd(base);
    
    %create my folders
    cellfun(@(x) mkdir([base,'MEEG/'],x),{'Preprocess','T1','Trials','EmptyRoom','MaxfilterOutput','Rest','ICAOutput','ScalpTimeStats'});
    
    %% import MEG data
    %locate raw files folder
    source      = ID.meg_folder;
    
    %locate and copy raw files into my directories
    sourcefiles = cellfun(@(x) sprintf('%s.fif',x),ID.meg_runs, 'UniformOutput',false);
    runlabs     = ID.meg_labs;

    runfiles    = cellfun(@(x) sprintf([base,'MEEG/%s_raw.fif'],x), runlabs,'UniformOutput',false);

    cd(source);
    cellfun(@copyfile,sourcefiles,runfiles,'UniformOutput',false);
    
%     if ~isempty(restfile)
%         copyfile(restfile{1},[base,'MEEG/Rest']);
%     end

%% import MRI structurals
if ~isempty(ID.mri)
% locate the right destination folder
   T1Folder = [base,'MEEG/T1/'];

% locate MRI images
   mriFolder = dir(ID.mri_folder); mriFolder = ['/mridata/cbu/',mriFolder.name];
   
   
   
   mriSubFld = dir(sprintf([mriFolder,'/%s_*'],ID.date_mri)); 
   
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
   
   cd([base,'MEEG']); save([ID.subj,'_filesID.mat'],'runfiles','sourcefiles','ID');  
    

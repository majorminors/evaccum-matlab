function a2_maxFilter(thisSubject)


overwrite = 0; % turn this on for auto overwrite

addpath /hpc-software/matlab/cbu/
rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab/';
addpath(genpath(fullfile(rootdir,'tools_analysis')))
droot = fullfile(rootdir,'data','meg_pilot_3',thisSubject.id);
maxfld= fullfile(droot,'MaxfilterOutput');
addpath(droot);



settings.overwrite = overwrite;

% get the names of the .fif files in the subjects' data folder
d = getnames(droot,[],'*raw.fif');

settings.maxfld = maxfld;


if~exist(fullfile(droot,['subjectInfo_',num2str(thisSubject.id),'.mat']),'file')
    save(fullfile(droot,['subjectInfo_',num2str(thisSubject.id),'.mat']),'thisSubject');% save subject's details
end

if ~isempty(thisSubject.bad_meg) %either use operator's bad channels or automatic bad channels detection
    if numel(thisSubject.bad_meg)>9; thisSubject.bad_meg = thisSubject.bad_meg(1:9);end %maxfilter doesn't work with more than 10 bad chans
    settings.badchans=['-autobad off -bad', sprintf(' %d',thisSubject.bad_meg)];
else
    settings.badchans = sprintf(' -autobad %d -badlimit %d',900,7);
end


for ifile = 1:length(d)
    settings.infname{ifile} = [droot,filesep,d{ifile}];
    %always match output's name to input's file to avoid mislabeling
    settings.outfname{ifile}= [droot,filesep,['Maxfiltered_' d{ifile}]];
end

[PATHSTR,NAME,EXT] = fileparts(settings.outfname{ifile});

if overwrite || ~exist([PATHSTR,filesep,NAME,'_trans',EXT],'file')
    doMaxFilter(settings)
end



return
end

function a2_maxFilter(thisSubject,datadir,toolsdir,overwrite)

warning('maxfilter will only work on f, g, and h nodes (i.e. compute nodes)')
% no idea how to force the script to use them though

if ~exist('overwrite','var'); overwrite = 0; end

addpath /hpc-software/matlab/cbu/
addpath(genpath(toolsdir))
droot = fullfile(datadir,thisSubject.id);
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
    infileNoRaw = strjoin(regexp(d{ifile},'_raw','split'),''); % split on the 'raw' filename bit that the MEG operators add and put it back together again
    settings.outfname{ifile}= [maxfld,filesep,infileNoRaw];
end

[PATHSTR,NAME,EXT] = fileparts(settings.outfname{ifile}); % just check one of these exists
if overwrite || ~exist([PATHSTR,filesep,NAME,'_trans',EXT],'file') % and if not (or overwrite)
    doMaxFilter(settings) % max it up!
end



return
end

%%
clear all
clear classes

close all
clc
addpath /hpc-software/matlab/cbu/

addpath(genpath('/group/woolgar-lab/projects/Dorian/evaccum/evaccum-matlab/tools_analysis/'))

droot = '/group/woolgar-lab/projects/Dorian/evaccum/evaccum-matlab/data/meg_pilot_1/megdata/';
maxfld= '/group/woolgar-lab/projects/Dorian/evaccum/evaccum-matlab/data/meg_pilot_1/megdata/%s/MEEG/MaxfilterOutput/'; % will this work
addpath(droot);


[dname,IDnum] = getnames(droot,7);




%% Prepare for parallel processing
clc





%%

overwrite = 0;

J = [];
ind = 0;

for subi = 1:length(dname)
    disp(subi)
    clear settings
    
    settings.overwrite = overwrite;
    
    d = getnames([droot,dname{subi},'/MEEG/'],[],'*raw.fif');
    
    settings.maxfld= sprintf(maxfld,dname{subi});
    
    
    ID = getMEGID(sprintf('DM_evaccumpilot_%s',IDnum{subi}));
    
    if~exist([droot,dname{subi},'ID_',num2str(dname{subi}),'.mat'],'file')
        save([droot,dname{subi},'ID_',num2str(dname{subi}),'.mat'],'ID');% save subject's details
    end
    
    if ~isempty(ID.bad_meg) %either use operator's bad channesl or automatic bad channels detection
        if numel(ID.bad_meg)>9; ID.bad_meg = ID.bad_meg(1:9);end %maxfilter doesn't work with more than 10 bad chans
        settings.badchans=['-autobad off -bad', sprintf(' %d',ID.bad_meg)];
    else
        settings.badchans = sprintf(' -autobad %d -badlimit %d',900,7);
    end
    
    
    for ifile = 1:length(d)
        settings.infname{ifile} = [droot,dname{subi},'/MEEG/', d{ifile}];
        %always match output's name to input's file to avoid mislabeling
        settings.outfname{ifile}= sprintf([maxfld,'Run%s%s.fif'],dname{subi},char(regexp(d{ifile},'(\d*_)*','match')),dname{subi});
    end
    
    [PATHSTR,NAME,EXT] = fileparts(settings.outfname{ifile});
    
    if overwrite || ~exist([PATHSTR,'/',NAME,'_trans',EXT],'file')
        fun_MaxFilter_RDKPD(settings)
    end
    
    
end


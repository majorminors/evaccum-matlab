function a5_visualise_data(thisSubject,datadir,toolsdir,manual)

error('not even really started')
% you need to be able to do trial definition
% and you can use ft_redefine trial to lock to response to visualise that
% too with ft_timelock or whatever


%% set up a fieldtrip config
cfg = [];
cfg.continuous = 'yes'; % this data is not epoched, it's continuous
cfg.dataset = outputFile;
cfg.trialfun = 'trlFromFile'; % this is my own trial function that simply reads in from the spm/osl trl file I made in a3_megTriggers
cfg.filename = sprintf(trlfile,num2str(runNum)); % I just arbitrarily add this to cfg struct so the function has something to read to find the file it needs
cfg = ft_definetrial(cfg);

end
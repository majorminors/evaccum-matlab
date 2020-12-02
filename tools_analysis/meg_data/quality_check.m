addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/'));
addpath(genpath('/group/woolgar-lab/projects/Dorian/EvAccum/tools_analysis/'))


droot = '/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/megdata';
dfld  = '/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/megdata/%s/MEEG/Preprocess/Run1_%s_trans.mat';
dfld2  = '/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/megdata/%s/MEEG/Preprocess/MSL_%s.mat';

savefl  = '/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/megdata/%s/MEEG/Preprocess/quality_%s.mat';


[dname, IDnum] = getnames(droot,7);

chantypes = {'EEG','MEGMAG','MEGPLANAR'};

%plot spectral power from signals in each modality
for subi = 3%1:numel(dname)
   ID = getMEGID(sprintf('DM_evaccumpilot_%s',IDnum{subi}));
   
   if ID.usable 
    fname = sprintf(dfld,dname{subi},dname{subi});
    D = spm_eeg_load(fname);
    clf
    for mods = 1:3
    idx = D.indchantype(chantypes{mods});
    %idx(ismember(idx,D.badchannels)) = []
    
    
    data = D(idx,:,:);
    
    
    data = data';
    %deglitch
    data = medfilt1(data,1,[],1);
    %detrend
    data = detrend(data);
    nwin = 500;
    nfft = 500;
    fs = D.fsample;
    pwelch(mean(data(1:100000,:),2), nwin, 0.5*nwin, nfft, fs)
    hold on;
    end
    axesHandlesToAllLines = findobj(gca,'Type','line');
    axesHandlesToAllLines(1).Color = 'b';
    axesHandlesToAllLines(2).Color = 'r';
    axesHandlesToAllLines(3).Color = 'k';
    legend({'EEG','MEGMAG','MEGPLANAR'})
    title(dname{subi})
    axis square
    set(gcf,'Color','white')
    
    export_fig( sprintf('/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/megdata/quality_check/%s_pow.jpeg',dname{subi}),'-transparent');
   
   end
   
end
close all
figure75= figure(75);
set(figure75,'Render','OpenGL','Units','pixels','Color',[.95 .95 .95],'Position',[78 391 1265 492],'PaperPosition',[-3.53  1.23 15.56  8.54])
for subi = 3%1:numel(dname)
    
   ID = getMEGID(sprintf('DM_evaccumpilot_%s',IDnum{subi}));
   
   if ID.usable 
   fname = sprintf(dfld2,dname{subi},dname{subi});
    D = spm_eeg_load(fname);
    clf
    for mods = 1:3
    idx = D.indchantype(chantypes{mods});
    idx(ismember(idx,D.badchannels)) = [];
    
    
    data = D(idx,:,:);
    data = medfilt1(data,1,[],3);
    subplot(1,3,mods)
    plot(D.time,median(data,3),'-k')
    grid on
    axis square
    xlabel('time (s)')
    ylabel('ERP/ERF')
    hold on;
    at_line(-1,1,'r')
    at_line(0,1,'g')
    title([dname{subi},'-',chantypes{mods}])
    set(gcf,'Color','white')
    end
export_fig( sprintf('/group/woolgar-lab/projects/Dorian/EvAccum/data/meg_pilot_1/megdata/quality_check/%s_ERPs.jpeg',dname{subi}),'-transparent');
   end
end

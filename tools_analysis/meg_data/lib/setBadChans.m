function setBadChans(fname,subj)

%bad channels
bad.TG = [3,67,68];
bad.LF = [7,8,3,22,25,27,70];
bad.SS = [3,7,8,66,70];
bad.IM = [3,4,37,46,45];
bad.EW = [29,61,63];
bad.JR = [12 61 63 67];
bad.SG = [3,8,29,66];
bad.EK = [16,40,66,67];
bad.IF = [3];
bad.JC = [2,38,39,56];
bad.BF = [3,46];
bad.MS = [];
bad.PK = [3,16,59];
bad.CM =[31, 65];
bad.JS = [51,69];
bad.KC = [3,45,68];
bad.LS = [56,69];
bad.LW = [2];


eval(sprintf(['badchan=bad.%s;'],subj));


%base = ['/imaging/at07/Matlab/Projects/CBU2015/RDKUnc/MEGData/',subj,'/MEEG/Preprocess/'];
%fname = ['_trdef_cRun',num2str(run),'_',subj,'_trans.mat'];

D = spm_eeg_load(fname);
ind = D.indchantype('EEG');

%remove old badchannels
D=D.badchannels(D.badchannels,0);

%update with new badchannels
D=D.badchannels(ind(badchan),1);

D.save;

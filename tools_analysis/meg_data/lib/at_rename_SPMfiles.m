droot = '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/MEGData/';
dfld = '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/MEGData/%s/MEEG/Preprocess/';

[dname,IDnum] = getnames(droot,7);
IDnum{1} ='04';IDnum{2} ='05';IDnum{3} ='08';
%%
for isub = 1:numel(dname)
    
    wrongfname = sprintf([dfld 'MSL_AT_RDK2_%s.mat'],dname{isub},num2str(IDnum{isub}));
    rightfname = sprintf([dfld 'MSL_%s.mat'],dname{isub},dname{isub});
    
    if exist(wrongfname)
        S =[];
        S.D = wrongfname;
        S.outfile = rightfname;
        D = spm_eeg_copy(S);        
        
        wrongfname = sprintf([dfld 'SL_AT_RDK2_%s.mat'],dname{isub},num2str(IDnum{isub}));
        rightfname = sprintf([dfld 'SL_%s.mat'],dname{isub},dname{isub});      
        S = [];
        S.D = wrongfname;
        S.outfile = rightfname;
        D = spm_eeg_copy(S);     
        
        wrongfname1 = sprintf([dfld 'MSL_AT_RDK2_%s.*'],dname{isub},num2str(IDnum{isub}));
        wrongfname2 = sprintf([dfld 'SL_AT_RDK2_%s.*'],dname{isub},num2str(IDnum{isub}));
        delete(wrongfname1) 
        delete(wrongfname2)

    end
end



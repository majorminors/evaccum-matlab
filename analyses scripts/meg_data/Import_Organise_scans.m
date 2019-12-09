
addpath /hpc-software/matlab/cbu/
addpath /imaging/at07/Matlab/Projects/CBU2016/RDK_PD/Preprocessing/

%numbers refer to the PID number
subjs = [51 52 53];
fld_tar = '/imaging/at07/Matlab/Projects/CBU2016/RDK_PD/MEGData/subj';

%decide whether to update (0 - e.g. after scanning new participants) or whether
%to start from scratch (1);
overwrite = 0;

% if not overwrite select only the new subjects
if ~overwrite
    exdir=arrayfun(@(x) exist([fld_tar,num2str(x),'/'],'dir'),subjs);
    subjs =subjs(exdir==0);
end


ParType = 2;   % Run on multiple Compute machines using parfar (best, but less feedback if crashes)

% open matlabpool if required
% matlabpool close force CBU_Cluster
if ParType
    if  isempty(gcp('nocreate'))%matlabpool('size')==0;
        P = cbupool(16);
        P.ResourceTemplate='-l nodes=^N^,mem=22GB,walltime=48:00:00';
        parpool(P);
    end
end

baseF ={};
parfor sb = 1:numel(subjs)

       baseF{sb} = [fld_tar,num2str(subjs(sb)),'/'];

       ID = getMEGID(sprintf('AT_RDK2_%s',num2str(subjs(sb))));
       tidyup_rdk2(baseF{sb},ID,overwrite);
       tidyrest_rdk2(baseF{sb},ID,overwrite);
end    

%parpool close force CBU_Cluster
delete(gcp)



addpath /hpc-software/matlab/cbu/
addpath /group/woolgar-lab/projects/EvAccum/data/meg_pilot_1/

%numbers refer to the PID number
subjs = [01 02 03];
fld_tar = '/group/woolgar-lab/projects/EvAccum/data/meg_pilot_1/S';

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

       ID = getMEGID(sprintf('DM_evaccumpilot_%s',num2str(subjs(sb))));
       tidyup_evaccum(baseF{sb},ID,overwrite);

end    

%parpool close force CBU_Cluster
delete(gcp)


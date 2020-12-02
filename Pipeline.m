%pipeline
%1 - Import_Organise_scans
%2 - Run_maxfilterRDK_PD -> somehow it fails when ran through the cluster
%    but it works fine on parpool
%3 - ConvertT1_PD + manual re-alignment
%4 - Run_MEGtrg_local
%5 - Run_RDK_preproc_PD
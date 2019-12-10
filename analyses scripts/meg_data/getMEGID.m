function ID = getMEGID(subject)
% retrieve scans IDs and details for the PD rdk study.
% meg .fif
% Ale's data (example): /megdata/cbu/rdk/megid/megdate/run
% windows: \\cbsu\data\Scandata\MEG\cbu\evaccum\megid\megdate\run
% unix: /megdata/cbu/evaccum/megid/megdate/run

switch subject
    
    case {'AT_RDK2_4','AT_RDK2_04'}
        
        ID.meg_fld  = 'meg16_0354';
        ID.meg_runs = {'rdk_1_2_raw','rdk_3_4_raw','rdk_5_6_raw','rdk_7_8_raw',...
            'rdk_9_10_raw','rdk_11_12_raw','rdk_13_14_raw','rdk_15_16_raw','rdk_17_18_raw',...
            'rdk_19_20_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.meg_rest = 'rdk_rest_raw';
        ID.meg_rlab = 'Run_rest';
        ID.mri      = 'CBU170666';
        ID.date_mri = '20171011';
        
        ID.bad_eeg = [52 73];
        ID.bad_meg = [0222 0813 1512];
        ID.date_meg= '170202';
        ID.emg     = 1;
        ID.tremor  ='left';
        ID.ctr     = 0;
        ID.tremor_s= 1+1+1+2+1+2+3;
        ID.age     = 66;
        ID.usable = 1;
        ID.motor  = 54;
        ID.LEDD   = 750;
        
        ID.subj  = 'AT_RDK2_04';
        
    case {'AT_RDK2_5','AT_RDK2_05'}
        
        ID.meg_fld  = 'meg17_0026';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU170329';
        ID.date_mri = '20170524';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
        
        ID.bad_eeg = [2 16 32 73];
        ID.bad_meg = [1412 0813];
        ID.date_meg= '170214';
        ID.emg     = 1;
        ID.note    = 'EMG recorded, however this was a control';
        ID.ctr     = 1;
        ID.age     = 51;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_05';
        
 case {'AT_RDK2_6','AT_RDK2_06'}
        
        ID.meg_fld  = 'meg17_0024';
        ID.meg_runs = {};
        ID.meg_labs = {};
        ID.mri      = '';
        ID.date_mri = '';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
        
        ID.bad_eeg = [2 16 32 73];
        ID.bad_meg = [0813 0640 1412];
        ID.date_meg= '170214';
        ID.emg     = 1;
        ID.note    = 'Aborted as volunteer could not keep their hand on the button box';
        ID.ctr     = 0;
        ID.tremor_s= 0+0+1+1+2+2+4; 
        ID.usable = 0;
        ID.subj  = 'AT_RDK2_06';  
        
    case {'AT_RDK2_8','AT_RDK2_08'}
        ID.meg_fld  = 'meg17_0090';
        ID.meg_runs = {'RDK2_2_raw','RDK2_2real_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU170263';
        ID.date_mri = '20170427';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [30 41 69 70 71 72 73];
        ID.bad_meg = [0222 0813];
        ID.date_meg= '170427';
        ID.emg     = 1;
        ID.ctr     = 0;
        ID.tremor_s= 0+0+1+1+0+1+1;      
        ID.age     = 57;
        ID.usable = 1;
        ID.motor  = 26;
        ID.LEDD   = 120;
        ID.subj  = 'AT_RDK2_08';
        
    case 'AT_RDK2_10'    
        ID.meg_fld  = 'meg17_0109';
        ID.meg_runs = {'RDK2_1_raw2','RDK2_2_raw','rdk2_3_raw','rdk2_4_rawreal',...
            'rdk2_5_raw','rdk2_6_raw','rdk2_7_raw','rdk2_8_raw','rdk2_9_raw',...
            'rdk2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU170326';
        ID.date_mri = '20170523';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        
        ID.bad_eeg = [12 14 16 17 23 38 39 40 51];
        ID.bad_meg = [0222 0813];
        ID.date_meg= '170523';
        ID.emg     = 1;
        ID.tremor  = 'left';
        ID.ctr     = 0;
        ID.tremor_s=0+1+1+1;
        ID.age     = 71;
        ID.motor  = 20;
        ID.LEDD   = 300;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_10';
        
    case 'AT_RDK2_11'
        ID.meg_fld  = 'meg17_0111';
        ID.meg_runs = {'RDK2_1_raw2','RDK2_2_raw2','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU170327';
        ID.date_mri = '20170523';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [4];
        ID.bad_meg = [1412 0813 2323 0222 1943 1712 1731];
        ID.date_meg= '170523';
        ID.emg     = 1;
        ID.note    = 'EMG recorded, however this was a control; suspiciously large ventricles but fine';
        ID.ctr     = 1;
        ID.age     = 72;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_11';
        
    case 'AT_RDK2_12'
        ID.meg_fld  = 'meg17_0125';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU170397';
        ID.date_mri = '20170703';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        
        ID.bad_eeg = [26 19 36 38 44 73 68 70];
        ID.bad_meg = [0222 1643 0813 1641];
        ID.date_meg= '170623';
        ID.emg     = 1;
        ID.ctr     = 0;
        ID.tremor_s= 0;
        ID.motor  = 42;
        ID.LEDD   = 450;
        ID.age     = 73;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_12';  
        
    case 'AT_RDK2_13'
        
        ID.meg_fld  = 'meg17_0127';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU170402';
        ID.date_mri = '20170704';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
             
        ID.bad_eeg = [16 17 38 39];
        ID.bad_meg = [0813 0933 1643];
        ID.date_meg= '170704';     
        ID.emg     = 1;
        ID.ctr     = 0;
        ID.tremor_s= 3;
        ID.tremor = 'both';
        ID.age     = 67;
        ID.motor  = 33;
        ID.LEDD   = 0;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_13';
        
    case 'AT_RDK2_15'
        
         
        ID.meg_fld  = 'meg17_0139';
        ID.meg_runs = {'RDK2_1b_raw','RDK2_2b_raw','RDK2_3b_raw','RDK2_4b_raw',...
            'RDK2_5b_raw','RDK2_6b_raw','RDK2_7b_raw','RDK2_8b_raw','RDK2_9b_raw',...
            'RDK2_10b_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU170447';
        ID.date_mri = '20170725';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        
        ID.bad_eeg = [24 71 72 73];
        ID.bad_meg = [0813 0222 1643 1731];
        ID.date_meg= '170725';     
        ID.emg     = 1;
        ID.notes   = 'HPIs moved during fitting the participant in the MEG and required extra checking; Maxfilter diagnostic did not flag any problem';
        ID.tremor  = 'right';
        ID.ctr     = 0;
        ID.tremor_s= 1+2+2;
        ID.age     = 60;
        ID.motor  = 17;
        ID.LEDD   = 2800;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_15';
        
    case 'AT_RDK2_16'
        
        ID.meg_fld  = 'meg17_0142';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU170455';
        ID.date_mri = '20170728';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [1 50];
        ID.bad_meg = [0813 0222 123 1643];
        ID.date_meg= '170728';     
        ID.emg     = 1;
        ID.note    = 'index button not working properly on first two blocks';
        ID.ctr     = 0;
        ID.tremor_s= 1;
        ID.tremor = 'right';        
        ID.age     = 55;
        ID.motor  = 15;
        ID.LEDD   = 75;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_16';
        
    case 'AT_RDK2_17'
        
           
        ID.meg_fld  = 'meg17_0164';
        ID.meg_runs = {'rdk2_1_2_raw_real','rdk2_3_4_raw','rdk2_5_6_raw','rdk2_7_8_raw',...
            'rdk2_9_10_raw','rdk2_11_12_raw','rdk2_13_14_raw','rdk2_15_16_raw','rdk2_17_18_raw',...
            'rdk2_19_20_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU170544';
        ID.date_mri = '20170831';
        ID.meg_rest = 'rdk2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [19 26 34 44 68 70 69 72 73 32];
        ID.bad_meg = [0813 0222 0121 2323];
        ID.date_meg= '170831';     
        ID.emg     = 1;
        ID.ctr     = 0;
        ID.tremor_s= 2;
        ID.tremor = 'right';    
        ID.age     = 61;
        ID.motor  = 22;
        ID.LEDD   = 90;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_17';
        
    case 'AT_RDK2_19'
        
            
        ID.meg_fld  = 'meg18_0060';
        ID.meg_runs = {'RDK2_1_raw_REAL','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU180345';
        ID.date_mri = '20180514';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [72 3];
        ID.bad_meg = [0813 1123 2122 0932 1412 1142];
        ID.date_meg= '180514';     
        ID.emg     = 1;
        ID.tremor  ='left';
        ID.ctr     = 0;
        ID.tremor_s= 1+1+1;
        ID.age     = 58;
        ID.motor  = 13;
        ID.LEDD   = 588;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_19';
        
    case 'AT_RDK2_20'
        
            
        ID.meg_fld  = 'meg19_0031';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190080';
        ID.date_mri = '20190205';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [58];
        ID.bad_meg = [0813];
        ID.date_meg= '190205';     
        ID.emg     = 0;
        ID.ctr     = 1;
        ID.age     = 64;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_20';
        
    case 'AT_RDK2_21'
            
        ID.meg_fld  = 'meg19_0060';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_9_raw','RDK2_9realb_raw',...
            'RDK2_4b_raw','RDK2_7b_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10','Run11'};
        ID.mri      = 'CBU190121';
        ID.date_mri = '20190219';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [58];
        ID.bad_meg = [0813 2543];
        ID.date_meg= '190219';     
        ID.emg     = 1;
        ID.notes   = '9 Runs, 10th run aborted';
        ID.ctr     = 0;
        ID.tremor_s= 3;
        ID.tremor = 'both';
        ID.age     = 53;
        ID.motor  = 19;
        ID.LEDD   = 295;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_21';
        
    case 'AT_RDK2_22'
            
        ID.meg_fld  = 'meg19_0055';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw','RDK2_9b_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10','Run11'};
        ID.mri      = 'CBU190134';
        ID.date_mri = '20190221';
        ID.meg_rest = 'RDK_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [55 58 68 71];
        ID.bad_meg = [0813 2141];
        ID.date_meg= '190218';     
        ID.emg     = 0;
        ID.ctr     = 1;
        ID.age     = 55;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_22';
        
    case 'AT_RDK2_24'    
        
           
        ID.meg_fld  = 'meg19_0073';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190150';
        ID.date_mri = '20190226';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [33 37 58];
        ID.bad_meg = [0813 2522];
        ID.date_meg= '190226';     
        ID.emg     = 1;
        ID.ctr     = 0;
        ID.tremor_s= 3;
        ID.tremor = 'both';
        ID.age     = 67;
        ID.motor  = 14;
        ID.LEDD   = 0;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_24';
        
    case 'AT_RDK2_25'
            
        ID.meg_fld  = 'meg19_0111';
        ID.meg_runs = {'rdk2_1_raw','rdk2_2_raw','rdk2_3_raw','rdk2_4_raw',...
            'rdk2_5_raw','rdk2_6_raw','rdk2_7_raw','rdk2_8_raw','rdk2_9_raw',...
            'rdk2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190215';
        ID.date_mri = '20190321';
        ID.meg_rest = '';%no resting state
        ID.meg_rlab = '';        
                
        ID.bad_eeg = [37 58];
        ID.bad_meg = [0813 2323];
        ID.date_meg= '190321';  
        ID.note_meg= 'rather strong right arm tremor - check emg';
        ID.emg     = 1;
        ID.ctr     = 0;
        ID.age     = 57;
        ID.motor  = 10;
        ID.LEDD   = 300;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_25';
        
    case 'AT_RDK2_28'
            
        ID.meg_fld  = 'meg19_0114';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190224';
        ID.date_mri = '20190325';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [58 71];
        ID.bad_meg = [0813 2323 2513];
        ID.date_meg= '190325';     
        ID.emg     = 1;
        ID.ctr     = 0;
        ID.tremor_s= 1+1+2+1+1+3;
        ID.tremor = 'left';
        ID.age     = 62;
        ID.motor  = 34;
        ID.LEDD   = 250;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_28';
        
    case 'AT_RDK2_29'
    
            
        ID.meg_fld  = 'meg18_0118';%note the mistake in the id
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190230';
        ID.date_mri = '20190326';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [58];
        ID.bad_meg = [0813 2323 1412];
        ID.date_meg= '190326';     
        ID.emg     = 1;
        ID.tremor = 'left';
        ID.ctr     = 0;
        ID.tremor_s= 1+2+3;
        ID.age     = 50;
        ID.motor  = 16;
        ID.LEDD   = 100;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_29';
        
    case 'AT_RDK2_30'
        
           
        ID.meg_fld  = 'meg19_0122';
        ID.meg_runs = {'RDK_2_1_raw','RDK2_2_raw','RD2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190238';
        ID.date_mri = '20190328';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [33 37 38 58 71];
        ID.bad_meg = [0813 2323 2032];
        ID.date_meg= '190328';     
        ID.emg     = 1;
        ID.tremor = 'left';
        ID.ctr     = 0;
        ID.tremor_s=1+1+2+2+3;
        ID.age     = 74;
        ID.motor  = 33;
        ID.LEDD   = 550;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_30';
        
    case 'AT_RDK2_31'
        
           
        ID.meg_fld  = 'meg19_0120';
        ID.meg_runs = {'rdk2_1_raw','rdk2_2_raw','rdk2_3_raw','rdk2_4_raw',...
            'rdk2_5_raw','rdk2_6_raw','rdk2_7_raw','rdk2_8_raw','rdk2_9_raw',...
            'rdk2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190239';
        ID.date_mri = '20190328';
        ID.meg_rest = 'rdk2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [37 49 58];
        ID.bad_meg = [0813 0122 0932];
        ID.date_meg= '190328';     
        ID.emg     = 0;
        ID.ctr     = 1;
        ID.age     = 70;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_31';
        
    case 'AT_RDK2_32'
            
        ID.meg_fld  = 'meg19_0127';
        ID.meg_runs = {'RDK2_1_raw_real','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDk2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190246';
        ID.date_mri = '20190401';
        ID.meg_rest = 'RDK2_1_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [37 52];
        ID.bad_meg = [0813 1412 2032];
        ID.date_meg= '190401';     
        ID.emg     = 0;
        ID.ctr     = 1;
        ID.age     = 72;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_32';
        
    case 'AT_RDK2_33'
            
        ID.meg_fld  = 'meg19_0132';
        ID.meg_runs = {'rdk2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190256';
        ID.date_mri = '20190404';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [58 39];
        ID.bad_meg = [0813 2323];
        ID.date_meg= '190404';     
        ID.emg     = 0;
        ID.ctr     = 1;
        ID.age     = 69;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_33';
        
     case 'AT_RDK2_34'
            
        ID.meg_fld  = 'meg19_148';
        ID.meg_runs = {'RDK2_1_RAW','RDK2_2_RAW','RDK2_3_RAW','RDK2_4_RAW',...
            'RDK2_5_RAW','RDK2_6_RAW','RDK2_7_RAW','RDK2_8_RAW','RDK2_9_RAW',...
            'RDK2_10_RAW'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190307';
        ID.date_mri = '20190425';
        ID.meg_rest = 'RDK2_REST_RAW';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [37 58];
        ID.bad_meg = [1422 1423 1433 0813];
        ID.date_meg= '190425';     
        ID.emg     = 1;
        ID.tremor  = 'left';
        ID.ctr     = 0;
        ID.tremor_s=1+1+1+2+3;
        ID.age     = 64;
        ID.motor  = 32;
        ID.LEDD   = 100;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_34';
        
    case 'AT_RDK2_35'
            
        ID.meg_fld  = 'meg19_0149';
        ID.meg_runs = {'rdk2_1_raw','rdk2_2_raw','rdk2_3_raw','rdk2_4_raw',...
            'rdk2_5_raw','rdk2_6_raw','rdk2_7_raw','rdk2_8_raw','rdk2_9_raw',...
            'rdk2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190319';
        ID.date_mri = '20190429';
        ID.meg_rest = 'rdk2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [37 58];
        ID.bad_meg = [2323 2022 0813];
        ID.date_meg= '190429';         
        ID.emg     = 0;
        ID.ctr     = 1;
        ID.age     = 69;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_35';
        
     case 'AT_RDK2_36'
            
        ID.meg_fld  = 'meg19_0152';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190328';
        ID.date_mri = '20190502';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [37 30];
        ID.bad_meg = [1412 2342 2112 0813];
        ID.date_meg= '190430';      
        ID.emg     = 0;
        ID.ctr     = 1;
        ID.age     = 68;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_36';
        
    case 'AT_RDK2_38'
            
        ID.meg_fld  = 'meg19_0156';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190294';
        ID.date_mri = '20190423';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = 37 ;
        ID.bad_meg = 0813;
        ID.date_meg= '190507';       
        ID.emg     = 0;
        ID.ctr     = 1;
        ID.age     = 75;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_38';
        
    case 'AT_RDK2_39'
            
        ID.meg_fld  = 'meg19_0158';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'rdk2_5_raw','rdk2_6_raw','rdk2_7_raw','rdk2_8_raw','rdk2_9_raw',...
            'rdk2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190344';
        ID.date_mri = '20190509';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [29 37 33 58];
        ID.bad_meg = [1412 0122 0123 0813];
        ID.date_meg= '190509';      
        ID.emg     = 0;
        ID.ctr     = 1;
        ID.age     = 64;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_39';
        
    case 'AT_RDK2_40'
            
        ID.meg_fld  = 'meg19_0164';
        ID.meg_runs = {'rdk2_1_raw','rdk2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190368';
        ID.date_mri = '20190516';
        ID.meg_rest = 'rdk2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [20 37 43 44 58] ;
        ID.bad_meg = 0813;
        ID.date_meg= '190516'; 
        ID.emg     = 0;
        ID.ctr     = 1;
        ID.age     = 75;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_40';
        
case 'AT_RDK2_41'
            
        ID.meg_fld  = 'meg19_0166';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190380';
        ID.date_mri = '20190520';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [37 58 29] ;
        ID.bad_meg = [0813 1731];
        ID.date_meg= '190520';   
        ID.emg     = 0;
        ID.ctr     = 1;
        ID.age     = 68;
        ID.note    = 'Volunteer lamented migraine during the scan';
        ID.usable = 0;
        ID.subj  = 'AT_RDK2_41';
        
case 'AT_RDK2_42'
            
        ID.meg_fld  = 'meg19_0168';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDk2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190391';
        ID.date_mri = '20190523';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [33 37] ;
        ID.bad_meg = [0813 2113 2032 2511];
        ID.date_meg= '190523';    
        ID.emg     = 0;
        ID.ctr     = 1;
        ID.age     = 71;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_42';
        
case 'AT_RDK2_43'
            
        ID.meg_fld  = 'meg19_0196';
        ID.meg_runs = {'rdk2_1_raw','rdk2_2_raw','rdk2_3_raw','rdk2_4_raw',...
            'rdk2_5_raw','rdk2_6_raw','rdk2_7_raw','rdk2_8_raw','rdk2_9_raw',...
            'rdk2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190466';
        ID.date_mri = '20190617';
        ID.meg_rest = 'rdk2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = 29 ;
        ID.bad_meg = 0813;
        ID.date_meg= '190617'; 
        ID.note    = {'first minute run #2 not usable: sensors re-tuned while recording','toilet break between 8th and 9th run'};
        ID.emg     = 1;
        ID.ctr     = 0;
        ID.tremor_s=1+1+1+1;
        ID.age     = 77;
        ID.motor  = 45;  
        ID.LEDD   = 450;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_43';
        
case 'AT_RDK2_44'
            
        ID.meg_fld  = 'meg19_0199';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190471';
        ID.date_mri = '20190618';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = 29 ;
        ID.bad_meg = [0813 1123 2323];
        ID.date_meg= '190618';         
        ID.emg     = 0;
        ID.ctr     = 1;
        ID.age     = 79;
        ID.usable = 0;
        ID.subj  = 'AT_RDK2_44';
        
case 'AT_RDK2_45'
            
        ID.meg_fld  = 'meg19_203';
        ID.meg_runs = {'rdk2_1_real_raw','rdk2_2_raw','rdk2_3_raw','rdk2_4_raw',...
            'rdk2_5_raw','rdk2_6_raw','rdk2_7_raw','rdk2_8_raw','rdk2_9_raw',...
            'rdk2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190477';
        ID.date_mri = '20190620';
        ID.meg_rest = 'rdk2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [] ;
        ID.bad_meg = [0813 1412 2323];
        ID.date_meg= '190620';             
        ID.emg     = 0;
        ID.ctr     = 1;
        ID.age     = 66;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_45';
        
case 'AT_RDK2_46'
            
        ID.meg_fld  = 'meg19_0226';
        ID.meg_runs = {'rdk2_1_raw','rdk2_2_raw_v2','rdk2_3_raw','rdk2_4_raw',...
            'rdk2_5_raw','rdk2_6_raw','rdk2_7_raw_v2','rdk2_8_raw','rdk2_9_raw',...
            'rdk2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190526';
        ID.date_mri = '20190711';
        ID.meg_rest = 'rdk2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [] ;
        ID.bad_meg = [0813];
        ID.date_meg= '190709';             
        ID.emg     = 0;
        ID.ctr     = 1;        
        ID.age     = 53;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_46';

case 'AT_RDK2_47'
            
        ID.meg_fld  = 'meg19_0230';
        ID.meg_runs = {'rdk2_1_raw','rdk2_2_raw','rdk2_3_raw','rdk2_4_raw',...
            'rdk2_5_raw','rdk2_6_raw','rdk2_7_raw','rdk2_8_raw','rdk2_9_raw',...
            'rdk2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190536';
        ID.date_mri = '20190712';
        ID.meg_rest = 'rdk2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [43:60];
        ID.bad_meg = [0813 1412 2323];
        ID.date_meg= '190712';             
        ID.emg     = 0;
        ID.ctr     = 1;        
        ID.note    = 'Damaged EEG cap cable 18 missing channels';
        ID.age     = 66;
        ID.usable = 0;
        ID.subj  = 'AT_RDK2_47';
        
case 'AT_RDK2_48'
            
        ID.meg_fld  = 'meg19_0241';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190551';
        ID.date_mri = '20190718';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [29 48] ;
        ID.bad_meg = [0813 0522 1913 2043];
        ID.date_meg= '190717';             
        ID.emg     = 1;
        ID.ctr     = 1;  
        ID.age     = 65;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_48';
        
case 'AT_RDK2_50'
            
        ID.meg_fld  = 'meg19_0243';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190549';
        ID.date_mri = '20190718';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [39 73 65];
        ID.bad_meg = [0813  2323];
        ID.date_meg= '190718';             
        ID.emg     = 1;
        ID.ctr     = 1;        
        ID.age     = 68;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_50';

case 'AT_RDK2_51'
            
        ID.meg_fld  = 'meg19_0381';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190850';
        ID.date_mri = '20191010';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [73];
        ID.bad_meg = [0813  2323 1123 1412];
        ID.date_meg= '191010';             
        ID.emg     = 1;
        ID.ctr     = 1;        
        ID.age     = 62;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_51';        

    case 'AT_RDK2_52'
            
        ID.meg_fld  = 'meg19_0387';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_5_raw','RDK2_6_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = '';
        ID.date_mri = '';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [];
        ID.bad_meg = [0813  0711];
        ID.date_meg= '191022';             
        ID.emg     = 1;
        ID.ctr     = 1;        
        ID.age     = 68;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_52';     
        
    case 'AT_RDK2_53'
            
        ID.meg_fld  = 'meg19_0393';
        ID.meg_runs = {'RDK2_1_raw','RDK2_2_raw','RDK2_3_raw','RDK2_4_raw',...
            'RDK2_6_raw','RDK2_6b_raw','RDK2_7_raw','RDK2_8_raw','RDK2_9_raw',...
            'RDK2_10_raw'};
        ID.meg_labs = {'Run1','Run2','Run3','Run4','Run5','Run6','Run7','Run8','Run9','Run10'};
        ID.mri      = 'CBU190076';
        ID.date_mri = '20190204';
        ID.meg_rest = 'RDK2_rest_raw';
        ID.meg_rlab = 'Run_rest';        
                
        ID.bad_eeg = [19 30];
        ID.bad_meg = [0813  2323];
        ID.date_meg= '191029';             
        ID.emg     = 1;
        ID.ctr     = 1;        
        ID.age     = 58;
        ID.usable = 1;
        ID.subj  = 'AT_RDK2_53';     
end






ID.meg_folder = sprintf('/megdata/cbu/rdk/%s/%s/',ID.meg_fld,ID.date_meg);
ID.mri_folder = sprintf('/mridata/cbu/%s_*',ID.mri);
if ~isfield(ID,'subj'); ID.subj = subject;end
%sanity check
if numel(ID.meg_runs) ~= numel(ID.meg_labs); error('Number of runs and labels does not match');end

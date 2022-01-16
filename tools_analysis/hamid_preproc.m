clc;
close all;
clear all;

ft_defaults;

toolboxdir = '/group/woolgar-lab/projects/Hamid/Projects/Vigilance/MEG/Analyses';
addpath([toolboxdir '/fieldtrip-lite-20190410']);
dataset_address='/megdata/cbu/mdconnectivity/';


%% start with an empty data struct and load MEG file
cfg = [];
photo_diode_based_triggers=0; % 1=photo_diode_based; 0=other_triggers
prestim=500;   % pre
posstim=1500;   % post
trigger_span=100;
baseline_span=1:100;
Subj_num=13; %
subject_code='211';
date='211214';

Data=[];
for run_num=1:1 % it can be 1:8
    clearvars triggers trigger_start_code ft ft_downsampled cfg data triggers triggers_code triggers_sent
    %% tackling incostincencies in saving the data
    run_str=num2str(run_num);
    trigger_shift=0;
    Run_char='run';
    %% Loading the data
    cfg.dataset = fullfile([dataset_address,'meg21_',subject_code,'/',date,'/',Run_char,run_str,'_raw.fif']);
    data = ft_preprocessing(cfg);
    
    %% Triggers
    % EOG  = 1 and 2
    % Data = 3 to 308
    % MISC = 309 to 320
    % Triggers = 321 to 328 (STI001 to STI008)
    % triggers 1:125: Stimulus appearance is 1:25*5 (125 codes) based on the phase of left stim
    % triggers 126:129: Cues in order about left SF, left Orient, right SF, right Orient
    % triggers 131:140: Probes 130+(1/2)*(1:5)
    % Button presses = 317 to 324 (STI009 to STI0016)
    % Photodiode trigger = 330; (STI010)
    % Combined triggers and button = 337 (STI101)
    % 338 = Sys201
    triggerchannels=[321:328]+trigger_shift;
    all_events_decimal_codes=bi2de([(data.trial{1}(triggerchannels,:)>4)'],'right-msb');
    triggers_code = all_events_decimal_codes>0 & all_events_decimal_codes<126;
    
    % Now we want to find the start of each trigger signal.
    trigger_start_code(1)=0;
    for count = 2:length(triggers_code)
        if triggers_code(count)>triggers_code(count-1) && triggers_code(count+trigger_span*0.8)>triggers_code(count-1)
            trigger_start_code(count)=1;
        else
            trigger_start_code(count)=0;
        end
    end
    
    % In triggers_sent we have all the positions a trigger started.
    triggers_sent = find(trigger_start_code);
    if length(triggers_sent)~=50
        err
    end
    %% Define/epoch trials
    cfg=[];
    cfg.trl = [triggers_sent-prestim;triggers_sent+posstim;repmat(-prestim,1,length(triggers_sent))]';
    
    % ft contains the cut data. One piece = one trial
    ft = ft_redefinetrial(cfg,data);
    
    %% removing unneeded channels
    for trl=1:length(ft.trial)
        ft.trial{1,trl}([1 2 309:end],:)=[];
    end
    c=0;
    for ch=3:308
        c=c+1;
        label2{c,1}=ft.label{ch,1};
    end
    ft.label=label2;
    
    %% Filtering an downsampling
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq   = [0.03 200];
    cfg.bpfilttype   = 'brickwall';
    cfg.bsfilter = 'yes';
    cfg.bsfreq   =[48 52];
    ft = ft_preprocessing(cfg,ft);
    
    cfg = [];
    cfg.resamplefs = 200;
    ft_downsampled = ft_resampledata(cfg,ft);
    
    
    %% Removing the baseline
    for tt=1:length(ft_downsampled.trial)
        for ch=1:size(ft_downsampled.trial{1,1},1)
            ft_downsampled.trial{1,tt}(ch,:)=squeeze(ft_downsampled.trial{1,tt}(ch,:))-repmat(nanmean(squeeze(ft_downsampled.trial{1,tt}(ch,baseline_span)),2),[1 size(ft_downsampled.trial{1,tt}(ch,:),2)]);
        end
    end
    
    Data=horzcat(Data,ft_downsampled.trial);
    [run_num]
end
%% Save data
% str_name='Stim';
%     save(['MEG_',sprintf('P%d',Subj_num) '_',str_name,'_aligned.mat'],'Data');
[Subj_num]



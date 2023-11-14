%%%%%%%%%%%%
%% set up %%
%%%%%%%%%%%%

clear all %#ok

% set up paths
addpath /hpc-software/matlab/cbu/
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219/');

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';
datadir = fullfile(rootdir,'data','meg_pilot_4'); addpath(datadir);
toolsdir = fullfile(rootdir,'tools_analysis');
ftDir = fullfile(toolsdir,'..','..','..','Toolboxes','fieldtrip');
jobdir = fullfile(rootdir,'job_logging');
addpath(ftDir); ft_defaults;
addpath(genpath(fullfile(toolsdir,'lib')))

toolbox = fullfile(rootdir,'..','..','Toolboxes','gramm'); addpath(toolbox); clear toolbox
% toolbox = fullfile(rootdir,'..','..','Toolboxes','bayesFactor'); addpath(toolbox); clear toolbox
toolbox = fullfile(rootdir,'..','..','Toolboxes','BFF_repo'); addpath(genpath(toolbox)); clear toolbox

erpFigDir = fullfile(datadir, 'erpFigs');
if ~exist(erpFigDir,'dir'); mkdir(erpFigDir); end

% we'll loop through subject data dirs, so just get the path from there
inputFileName = ['Preprocess' filesep 'timelocked_averages.mat'];

% some colours
teal = [0.2, 0.6, 0.7];
coral = [0.9, 0.4, 0.3];
lilac = [0.7, 0.5, 0.8];
maroon = [0.6579, 0.2821, 0.198];


%%%%%%%%%%%%%%%%%%%%%%%%
%% load and prep data %%
%%%%%%%%%%%%%%%%%%%%%%%%

% since our folders are all e.g. S01, S12 etc, let's load information about
% the directories that have that pattern
tmp = dir([datadir filesep 'S*']);
subjectFolders = tmp([tmp.isdir]); clear tmp
% and loop through them
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    disp('no parallel pool is currently initialised---initialising');
    workers = numel(subjectFolders);
    P = cbupool(workers);
    P.JobStorageLocation = jobdir;
    tempPath = fullfile(jobdir,'tmp');
    if ~exist(tempPath,'dir');mkdir(tempPath);end
    P.JobStorageLocation = tempPath;
    parpool(P,workers);
else
    disp('parallel pool has already been initialized---skipping');
end
parfor subjectNum = 1:numel(subjectFolders)
    % clear thisFile % can't do this in a parfor loop
    
    % full path to the inputfile
    thisFile = fullfile(subjectFolders(subjectNum).folder,subjectFolders(subjectNum).name,inputFileName);
    % skip this loop if we can't find one
    if ~exist(thisFile,'file'); continue; end
    % I seem to be doing a lot of checks in this loop and I can't remember why
    % so let's do another one to
    % account for the fact I can't clear the thisFile variable in a parfor
    % loop. let's check the subject directory part of thisFile matches the
    % subject number
    pathParts = strsplit(thisFile, '/');
    index = find(contains(pathParts, 'S'));
    subjectCode = pathParts{index}; %#ok
    if ~strcmp(subjectCode,subjectFolders(subjectNum).name); error('file doesnt match subject'); end
    
    fprintf('this is subject: %s...',subjectFolders(subjectNum).name)
    
     % grab info about all the meeg data files we care about
    theseFiles = dir([subjectFolders(subjectNum).folder filesep subjectFolders(subjectNum).name filesep inputFileName]);
    if isempty(theseFiles); continue; end
    
    % so here, load what you care about and play with them!
    fprintf('loading...')

    erpManips = {...
        'ec_responseLockedAverage' 'hc_responseLockedAverage'...
        'er_responseLockedAverage' 'hr_responseLockedAverage'...
        'ec_coherenceLockedAverage' 'hc_coherenceLockedAverage'...
        'er_coherenceLockedAverage' 'hr_coherenceLockedAverage'...
    };
    erpConds = {...
        'ecer_responseLockedAverage' 'ecer_coherenceLockedAverage'...
        'echr_responseLockedAverage' 'echr_coherenceLockedAverage'...
        'hcer_responseLockedAverage' 'hcer_coherenceLockedAverage'...
        'hchr_responseLockedAverage' 'hchr_coherenceLockedAverage'...
    };

    whichVars = {...
            erpManips{:} erpConds{:}...
        };
    
    data{subjectNum} = load(thisFile,whichVars{:});
    
    fprintf('loaded\n')
    
end; clear theseFiles thisFile subjectFolders subjectNum subjectCode index pathParts
clear erpManips erpConds
response = input('Do you want to delete the parallel pool? (y/n): ','s');
if strcmpi(response,'y') || strcmpi(response,'yes')
    disp('deleting...');
    delete(gcp('nocreate')); clear workers;
    system(['rm -rf ' tempPath]);
    if exist(tempPath,'file'); warning('I couldnt delete the job directory :('); end
else
    disp('not deleting...')
end

fprintf('loading complete with %.0f subjects\n',sum(~cellfun(@isempty,data)))

% get some layouts and calculate neighbours

disp('prep layouts and neighbours')

% meg is standard---can just use FTs version
megLayout = fullfile(ftDir,'template','layout','neuromag306all.lay');
tmp = load(fullfile(ftDir,'template','neighbours','neuromag306mag_neighb.mat'),'neighbours');
megNeighbours = tmp.neighbours; clear tmp
% or make our own from a subject:
% megLayout = fullfile(datadir,'S01','Preprocess','meg_layout.lay'); % note this layout will squeeze everything into an eeg circle in a topo plot---hard to read
% cfg = [];
% cfg.channel = 'MEG';
% cfg.method = 'distance';
% cfg.layout = megLayout;
% megNeighbours = ft_prepare_neighbours(cfg);

% for EEG, we have a modified 70channel cap to fit 64 channels, so although
% there is a template for 64 chan eeg e.g.:
% eegLayout = fullfile(ftDir,'template','layout','acticap-64ch-standard2.mat');
% eegNeighbours = fullfile(ftDir,'template','neighbours','easycap64ch-avg_neighb.mat');
% I think better to make this ourselves from a participant layout
eegLayout = fullfile(datadir,'S01','Preprocess','eeg_layout.lay');

% now neighbours

cfg = [];
cfg.channel = 'EEG';
cfg.method = 'distance';
cfg.layout = eegLayout;
eegNeighbours = ft_prepare_neighbours(cfg);


% and we'll collect useful sensors

% % we can plot these
% cfg = [];
% cfg.layout = eegLayout;
% % cfg.layout = megLayout;
% layout = ft_prepare_layout(cfg);
% % ft_plot_layout(layout);
% % clear layout


CPP = {'EEG040' 'EEG041' 'EEG042'};
frontal = {'EEG004' 'EEG002' 'EEG008'};

disp('done loading and prep')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% let's compile some averages and differences %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------%
%-- average coherence data --%
%----------------------------%

disp('averaging coherence data')

cfg = [];
ecRespAll = returnStructs(data, 'ec_responseLockedAverage');
ecRespAve = ft_timelockgrandaverage(cfg, ecRespAll{:});
cfg = [];
ecOnsAll = returnStructs(data, 'ec_coherenceLockedAverage');
ecOnsAve = ft_timelockgrandaverage(cfg, ecOnsAll{:});
cfg = [];
hcRespAll = returnStructs(data, 'hc_responseLockedAverage');
hcRespAve = ft_timelockgrandaverage(cfg, hcRespAll{:});
cfg = [];
hcOnsAll = returnStructs(data, 'hc_coherenceLockedAverage');
hcOnsAve = ft_timelockgrandaverage(cfg, hcOnsAll{:});

% get an overall average
cfg = [];
cOnsAve = ft_timelockgrandaverage(cfg, ecOnsAll{:}, hcOnsAll{:});
cfg = [];
cRespAve = ft_timelockgrandaverage(cfg, ecRespAll{:}, hcRespAll{:});

% get diffs
cOnsDiffAve = getDifference(ecOnsAve,hcOnsAve);
cRespDiffAve = getDifference(ecRespAve,hcRespAve);

disp('done')

%---------------------------------%
%-- average categorisation data --%
%---------------------------------%

disp('averaging categorisation data')

cfg = [];
erRespAll = returnStructs(data, 'er_responseLockedAverage');
erRespAve = ft_timelockgrandaverage(cfg, erRespAll{:});
cfg = [];
erOnsAll = returnStructs(data, 'er_coherenceLockedAverage');
erOnsAve = ft_timelockgrandaverage(cfg, erOnsAll{:});
cfg = [];
hrRespAll = returnStructs(data, 'hr_responseLockedAverage');
hrRespAve = ft_timelockgrandaverage(cfg, hrRespAll{:});
cfg = [];
hrOnsAll = returnStructs(data, 'hr_coherenceLockedAverage');
hrOnsAve = ft_timelockgrandaverage(cfg, hrOnsAll{:});

% get an overall average
cfg = [];
rOnsAve = ft_timelockgrandaverage(cfg, erOnsAll{:}, hrOnsAll{:});
cfg = [];
rRespAve = ft_timelockgrandaverage(cfg, erRespAll{:}, hrRespAll{:});

% get diffs
rOnsDiffAve = getDifference(erOnsAve,hrOnsAve);
rRespDiffAve = getDifference(ecRespAve,hcRespAve);

disp('done')

%--------------------------------%
%-- average conditionwise data --%
%--------------------------------%

disp('averaging conditionwise data')

cfg = [];
ecerRespAll = returnStructs(data, 'ecer_responseLockedAverage');
ecerRespAve = ft_timelockgrandaverage(cfg, ecerRespAll{:});
cfg = [];
ecerOnsAll = returnStructs(data, 'ecer_coherenceLockedAverage');
ecerOnsAve = ft_timelockgrandaverage(cfg, ecerOnsAll{:});
cfg = [];
echrRespAll = returnStructs(data, 'echr_responseLockedAverage');
echrRespAve = ft_timelockgrandaverage(cfg, echrRespAll{:});
cfg = [];
echrOnsAll = returnStructs(data, 'echr_coherenceLockedAverage');
echrOnsAve = ft_timelockgrandaverage(cfg, echrOnsAll{:});
cfg = [];
hcerRespAll = returnStructs(data, 'hcer_responseLockedAverage');
hcerRespAve = ft_timelockgrandaverage(cfg, hcerRespAll{:});
cfg = [];
hcerOnsAll = returnStructs(data, 'hcer_coherenceLockedAverage');
hcerOnsAve = ft_timelockgrandaverage(cfg, hcerOnsAll{:});
cfg = [];
hchrRespAll = returnStructs(data, 'hchr_responseLockedAverage');
hchrRespAve = ft_timelockgrandaverage(cfg, hchrRespAll{:});
cfg = [];
hchrOnsAll = returnStructs(data, 'hchr_coherenceLockedAverage');
hchrOnsAve = ft_timelockgrandaverage(cfg, hchrOnsAll{:});

disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% let's do an anova on the conditions %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('doing anova on conditions')

dataSets.onset = {ecerOnsAll echrOnsAll hcerOnsAll hchrOnsAll};
dataSets.response = {ecerRespAll echrRespAll hcerRespAll hchrRespAll};

% we can try erins, but this does a mixed model with subjects as a random
% effect: not quite an rm anova
% novas = fieldnames(dataSets);
% for i = 1:numel(novas)
%     thisNova = novas{i};
%     
%     fprintf('starting with %s...',thisNova)
%     
%     nSubjects = numel(dataSets.(thisNova){1});
%     nTimepoints = numel(dataSets.(thisNova){1}{1}.avg(1,:));
%     nConditions = numel({'ecer' 'echr' 'hcer' 'hchr'});
%     
%     fprintf('compiling...')
%     thisAnova = NaN(nSubjects,nTimepoints,nConditions); % initialise this
%     for subject = 1:nSubjects
%         for condition = 1:nConditions
%             % we have all timepoints, so we don't need to loop them
%             thisAnova(subject,:,condition) = getRawAmplitude(dataSets.(thisNova){condition}{subject}, CPP);
%         end; clear subject
%     end; clear condition
%     
%     fprintf('running anova...')
%     
%     [p, fdrh, stats] = doAnova(thisAnova);
%     results.(thisNova).p = p;
%     results.(thisNova).h = fdrh;
%     results.(thisNova).stats = stats;
%     
%     clear p fp t thisAnova nSubjects nTimepoints nConditions
%     
%     fprintf('done\n')
%     
% end; clear i thisNova
% 
% disp('done anova on conditions')
% 
% plot(1:size(results.onset.h,1),results.onset.h.difficulty_x_manipulation)
% plot(1:size(results.response.h,1),results.response.h.difficulty_x_manipulation)

% let's try an rm anova
novas = fieldnames(dataSets);
for i = 1:numel(novas)
    thisNova = novas{i};
    
    fprintf('starting with %s\n',thisNova)
    
    % some initialisations
    conditionNames = {'ecer' 'echr' 'hcer' 'hchr'};
    nSubjects = numel(dataSets.(thisNova){1});
    nTimepoints = numel(dataSets.(thisNova){1}{1}.avg(1,:));
    nConditions = numel(conditionNames);
    
    % let's grab the data per subject, per condition and compile it
    for subject = 1:nSubjects
        for condition = 1:nConditions
            % so this pulls out a row of mean amplitudes per timepoint
            % stack them row-wise for subjects and the 3rd dimension is
            % conditions
            tmp(subject,:,condition) = getRawAmplitude(dataSets.(thisNova){condition}{subject}, CPP);
        end; clear condition
    end; clear subject
    
    % let's reshape it now, so we have the conditions row-wise and
    % timepoints on the 3rd dimension
    % I don't actually remember why I did this---maybe it just made better
    % sense in my head?
    conditionwiseData = permute(tmp, [1, 3, 2]); clear tmp
    
    % couple more initialisations
    pVals = NaN(nTimepoints,3);
    meanSq = NaN(nTimepoints,3);
    fstat = NaN(nTimepoints,3);

    % now we loop through the timepoints and construct the anova
    for timepoint = 1:nTimepoints
        if timepoint>1; fprintf(repmat('\b', 1, numel(timepointStr))); end
        timepointStr = sprintf('timepoint %.0f of %.0f', timepoint, nTimepoints);
        fprintf(timepointStr);
        
        % arbitrary index numbers for the subjects so the model fitting can
        % identify the subject-specific effects
        subjects = 1:nSubjects;
        
        % now collect into a table the subject ids + their data for each
        % condition (second dimension) for this timepoint (3rd dimension)
        t = array2table([subjects',conditionwiseData(:,1,timepoint),conditionwiseData(:,2,timepoint),conditionwiseData(:,3,timepoint),conditionwiseData(:,4,timepoint)],...
            'VariableNames',{'Subject','ecer','echr','hcer','hchr'});
        % create a table that indicates how the levels of the two
        % manipulations vary across those conditions so we end up with:
        % 1,1: easy coherence, easy rule
        % 1,2: easy coherence, hard rule
        % 2,1: hard coherence, easy rule
        % 2,2: hard coherence, hard rule
        within = table(categorical([1;1;2;2]),categorical([1;2;1;2]),'VariableNames',{'Coherence','Rule'});
        % specify the design---each condition on the intercept (1)
        rm = fitrm(t,'ecer,echr,hcer,hchr ~ 1', 'WithinDesign', within);
        % and we want to look at the interaction, so specify that
        ranovatbl = ranova(rm, 'WithinModel','Coherence*Rule');
%         writetable(t)
        
        % the ranova table comes out with stats for the main effects and
        % the interaction now, but it's a bit weird, because it also
        % produces invisible values for the Error rows or
        % something? It looks like there are 4 p values in pVal, but
        % actually the fourth is always 0.5 and this is where there is no
        % value showing that corresponds to the Error row. Very odd.
        % Anyway, make sure to select the actual row and not the *apparent*
        % position. Here I loop though the actual row numbers.
        for testIdx = 1:3
            rowNumbers = [3,5,7]; % the correct row numbers
            pVals(timepoint,testIdx) = ranovatbl.pValue(rowNumbers(testIdx));
            meanSq(timepoint,testIdx) = ranovatbl.MeanSq(rowNumbers(testIdx));
            fstat(timepoint,testIdx) = ranovatbl.F(rowNumbers(testIdx));
        end; clear testIdx rowNumbers
        clear ranovatbl conds rm t subjects within
    end; clear timepoint
    fprintf('\n')
    
    % now Niko Kriegeskorte's FDR correction for that set of p-values
    % if we want to correct each test seperately
%     fdrp = NaN(1,3);
%     for pidx = 1:3
%         if ~isempty(FDRthreshold(pVals(:,pidx)))
%             fdrp(:,pidx) = FDRthreshold(pVals(:,pidx));
%         end
%     end; clear pidx
    % if we want to correct across tests instead
    fdrp = FDRthreshold(pVals);
    % now let's get where we can reject the null
    if isempty(fdrp); fdrp = 0; end
    fdrh = pVals < fdrp;
    varNames = {'Coherence', 'Rule', 'Interaction'};
    
    % and compile them into tables for ease of reference
    results.(thisNova).p = array2table(pVals, 'VariableNames', varNames);
    results.(thisNova).h = array2table(fdrh, 'VariableNames', varNames);
    results.(thisNova).meanSq = array2table(meanSq, 'VariableNames', varNames);
    results.(thisNova).fstat = array2table(fstat, 'VariableNames', varNames);
    clear pVals fdrh meanSq fstat
    
end; clear i thisNova
% 
% scatter(1:numel(results.onset.h.Interaction),results.onset.h.Interaction)
% scatter(1:numel(results.response.h.Interaction),results.response.h.Interaction)
% scatter(1:numel(results.onset.p.Interaction),results.onset.p.Interaction<0.05)
% scatter(1:numel(results.response.p.Interaction),results.response.p.Interaction<0.05)
    

clear dataSets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get subjectwise differences %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('getting differences')

for subject = 1:numel(ecRespAll)
    fprintf('\ngetting subjectwise differences for subject %.0f of %.0f\n\n',subject,numel(ecRespAll))
    
    cOnsDiffAll{subject} = getDifference(ecOnsAll{subject}, hcOnsAll{subject});
    cRespDiffAll{subject}= getDifference(ecRespAll{subject}, hcRespAll{subject});
    rOnsDiffAll{subject} = getDifference(erOnsAll{subject}, hrOnsAll{subject});
    rRespDiffAll{subject}= getDifference(erRespAll{subject}, hrRespAll{subject});
    
    ecOnsDiffAll{subject} = getDifference(ecerOnsAll{subject}, echrOnsAll{subject});
    hcOnsDiffAll{subject}= getDifference(hcerOnsAll{subject}, hchrOnsAll{subject});
    ttestOnsDiffAll{subject} = getDifference(ecOnsDiffAll{subject}, hcOnsDiffAll{subject});
    
    ecRespDiffAll{subject} = getDifference(ecerRespAll{subject}, echrRespAll{subject});
    hcRespDiffAll{subject}= getDifference(hcerRespAll{subject}, hchrRespAll{subject});
    ttestRespDiffAll{subject}= getDifference(ecRespDiffAll{subject}, hcRespDiffAll{subject});

end; clear subject

ttestOnsDiffAve = ft_timelockgrandaverage(cfg, ttestOnsDiffAll{:});
ttestRespDiffAve = ft_timelockgrandaverage(cfg, ttestRespDiffAll{:});
ttestOnsDiffAve.bfs = testDiffsAcrossTime(ttestOnsDiffAll,CPP);
ttestRespDiffAve.bfs = testDiffsAcrossTime(ttestRespDiffAll,CPP);


% disp('getting averages')
% 
% cfg = [];
% consdiffave = ft_timelockgrandaverage(cfg, consdiffall{:});
% crespdiffave = ft_timelockgrandaverage(cfg, crespdiffall{:});
% ronsdiffave = ft_timelockgrandaverage(cfg, ronsdiffall{:});
% rrespdiffave = ft_timelockgrandaverage(cfg, rrespdiffall{:});

disp('done')

%% save those differences, if we want

disp('saving')

save(fullfile(datadir,'erp_differences.mat'),...
    'cOnsDiffAll', 'cRespDiffAll', 'rOnsDiffAll', 'rRespDiffAll', 'cOnsDiffAve', 'cRespDiffAve', 'rOnsDiffAve', 'rRespDiffAve'...
    )

disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the timecourse of univariate activity %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%% coh onset in eeg
cOnsDiffAve.bfs = testDiffsAcrossTime(cOnsDiffAll,CPP);
plotTimecourse(eegLayout,[teal;coral],... % layout and colours
    {ecOnsAve,hcOnsAve},cOnsDiffAve,... % {average,data},averageDifference (with bayes factors)
    {ecOnsAll,hcOnsAll},... % {subjectwise,data}
    [-0.2 1.5],[],CPP,... % xlims,ylims,channels
    {'EEG','uV'},... % sensor type and units for ylablel
    {'easy coh';'hard coh'},'northwest','CPP (CP1 CPz CP2) in EEG: Coherence Onset',... % legend, legend location, plot title
    [erpFigDir filesep 'eeg_coh_ons_ERP.png']) % save loc
%% coh response in eeg
cRespDiffAve.bfs = testDiffsAcrossTime(cRespDiffAll,CPP);
plotTimecourse(eegLayout,[teal;coral],...
    {ecRespAve,hcRespAve},cRespDiffAve,...
    {ecRespAll,hcRespAll},...
    [],[],CPP,...
    {'EEG','uV'},...
    {'easy coh';'hard coh'},'northwest','CPP (CP1 CPz CP2) in EEG: Coherence Response',...
    [erpFigDir filesep 'eeg_coh_resp_ERP.png'])
%% rule onset in eeg
rOnsDiffAve.bfs = testDiffsAcrossTime(rOnsDiffAll,CPP);
plotTimecourse(eegLayout,[teal;coral],...
    {erOnsAve,hrOnsAve},rOnsDiffAve,...
    {erOnsAll,hrOnsAll},...
    [-0.2 1.5],[-3.5e-06 7.5e-06],CPP,...
    {'EEG','uV'},...
    {'easy cat';'hard cat'},'southeast','CPP (CP1 CPz CP2) in EEG: Categorisation Onset',...
    [erpFigDir filesep 'eeg_cat_ons_ERP.png'])
%% rule response in eeg
rRespDiffAve.bfs = testDiffsAcrossTime(rRespDiffAll,CPP);
plotTimecourse(eegLayout,[teal;coral],...
    {erRespAve,hrRespAve},rRespDiffAve,...
    {erRespAll,hrRespAll},...
    [],[-3.5e-06 7.5e-06],CPP,...
    {'EEG','uV'},...
    {'easy cat';'hard cat'},'northwest','CPP (CP1 CPz CP2) in EEG: Categorisation Response',...
    [erpFigDir filesep 'eeg_cat_resp_ERP.png'])
%% onset of all conditions in eeg
plotTimecourse(eegLayout,[lilac;teal;coral;maroon],... % layout and colours
    {ecerOnsAve, echrOnsAve, hcerOnsAve, hchrOnsAve},[],... % {average,data},averageDifference (with bayes factors)
    {ecerOnsAll,echrOnsAll,hcerOnsAll,hchrOnsAll},... % {subjectwise,data}
    [-0.2 1.5],[],CPP,... % xlims,ylims,channels
    {'EEG','uV'},... % sensor type and units for ylablel
    {'EasyCoh EasyCat';'EasyCoh HardCat';'HardCoh EasyCat';'HardCoh HardCat'},'northwest','CPP (CP1 CPz CP2) in EEG: Onset in all Conditions',... % legend, legend location, plot title
    [erpFigDir filesep 'eeg_allConds_ons_ERP.png'],... % save loc
    results.onset.h) % rsa significance
plotTimecourse(eegLayout,[lilac;teal;coral;maroon],... % layout and colours
    {ecerOnsAve, echrOnsAve, hcerOnsAve, hchrOnsAve},ttestOnsDiffAve,... % {average,data},averageDifference (with bayes factors)
    {ecerOnsAll,echrOnsAll,hcerOnsAll,hchrOnsAll},... % {subjectwise,data}
    [-0.2 1.5],[],CPP,... % xlims,ylims,channels
    {'EEG','uV'},... % sensor type and units for ylablel
    {'EasyCoh EasyCat';'EasyCoh HardCat';'HardCoh EasyCat';'HardCoh HardCat'},'northwest','CPP (CP1 CPz CP2) in EEG: Onset in all Conditions',... % legend, legend location, plot title
    [erpFigDir filesep 'eeg_allConds_ons_ERP_bayes.png']) % save loc

%% response of all conditions in eeg
plotTimecourse(eegLayout,[lilac;teal;coral;maroon],... % layout and colours
    {ecerRespAve, echrRespAve, hcerRespAve, hchrRespAve},[],... % {average,data},averageDifference (with bayes factors)
    {ecerRespAll,echrRespAll,hcerRespAll,hchrRespAll},... % {subjectwise,data}
    [],[],CPP,... % xlims,ylims,channels
    {'EEG','uV'},... % sensor type and units for ylablel
    {'EasyCoh EasyCat';'EasyCoh HardCat';'HardCoh EasyCat';'HardCoh HardCat'},'northwest','CPP (CP1 CPz CP2) in EEG: Response in all Conditions',... % legend, legend location, plot title
    [erpFigDir filesep 'eeg_allConds_resp_ERP.png'],... % save loc
    results.response.h) % rsa significance
plotTimecourse(eegLayout,[lilac;teal;coral;maroon],... % layout and colours
    {ecerRespAve, echrRespAve, hcerRespAve, hchrRespAve},ttestRespDiffAve,... % {average,data},averageDifference (with bayes factors)
    {ecerRespAll,echrRespAll,hcerRespAll,hchrRespAll},... % {subjectwise,data}
    [],[],CPP,... % xlims,ylims,channels
    {'EEG','uV'},... % sensor type and units for ylablel
    {'EasyCoh EasyCat';'EasyCoh HardCat';'HardCoh EasyCat';'HardCoh HardCat'},'northwest','CPP (CP1 CPz CP2) in EEG: Response in all Conditions',... % legend, legend location, plot title
    [erpFigDir filesep 'eeg_allConds_resp_ERP_bayes.png']) % save loc
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the topology of univariate activity %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% where is easy coh different from hard coh

% let's get some zlimits
        
minValue = @(x, pref) min(min(x.avg(startsWith(x.label, pref), :)));
maxValue = @(x, pref) max(max(x.avg(startsWith(x.label, pref), :)));

zlimEegOns = [...
    min([minValue(cOnsDiffAve,'EEG'),minValue(rOnsDiffAve,'EEG'),minValue(cOnsAve,'EEG'),minValue(rOnsAve,'EEG')]),...
    max([maxValue(cOnsDiffAve,'EEG'),maxValue(rOnsDiffAve,'EEG'),maxValue(cOnsAve,'EEG'),maxValue(rOnsAve,'EEG')])...
    ];
zlimEegResp = [...
    min([minValue(cRespDiffAve,'EEG'),minValue(rRespDiffAve,'EEG'),minValue(cRespAve,'EEG'),minValue(rRespAve,'EEG')]),...
    max([maxValue(cRespDiffAve,'EEG'),maxValue(rRespDiffAve,'EEG'),maxValue(cRespAve,'EEG'),maxValue(rRespAve,'EEG')])...
    ];
zlimMegOns = [...
    min([minValue(cOnsDiffAve,'MEG'),minValue(rOnsDiffAve,'MEG'),minValue(cOnsAve,'MEG'),minValue(rOnsAve,'MEG')]),...
    max([maxValue(cOnsDiffAve,'MEG'),maxValue(rOnsDiffAve,'MEG'),maxValue(cOnsAve,'MEG'),maxValue(rOnsAve,'MEG')])...
    ];
zlimMegResp = [...
    min([minValue(cRespDiffAve,'MEG'),minValue(rRespDiffAve,'MEG'),minValue(cRespAve,'MEG'),minValue(rRespAve,'MEG')]),...
    max([maxValue(cRespDiffAve,'MEG'),maxValue(rRespDiffAve,'MEG'),maxValue(cRespAve,'MEG'),maxValue(rRespAve,'MEG')])...
    ];

%% EEG topo of onset across coherence
doTopos({'EEG'},eegNeighbours,ecOnsAll,hcOnsAll,...
    zlimEegOns,0.05,[-0.3,1.5],[6,6],36,...
    cOnsDiffAve,eegLayout,'eeg_coh_ons_on_diff_sig_ERP_clusters','eeg_coh_ons_sig_ERP_clusters',datadir,erpFigDir)
%% EEG topo of response across coherence
doTopos({'EEG'},eegNeighbours,ecRespAll,hcRespAll,...
    zlimEegResp,0.05,[-0.6,0.2],[4,4],16,...
    cRespDiffAve,eegLayout,'eeg_coh_resp_on_diff_sig_ERP_clusters','eeg_coh_resp_sig_ERP_clusters',datadir,erpFigDir)
%% EEG topo of onset across categorisation
doTopos({'EEG'},eegNeighbours,erOnsAll,hrOnsAll,...
    zlimEegOns,0.05,[-0.3,1.5],[6,6],36,...
    rOnsDiffAve,eegLayout,'eeg_cat_ons_on_diff_sig_ERP_clusters','eeg_cat_ons_sig_ERP_clusters',datadir,erpFigDir)
%% EEG topo of response across categorisation
doTopos({'EEG'},eegNeighbours,erRespAll,hrRespAll,...
    zlimEegResp,0.05,[-0.6,0.2],[4,4],16,...
    rRespDiffAve,eegLayout,'eeg_cat_resp_on_diff_sig_ERP_clusters','eeg_cat_resp_sig_ERP_clusters',datadir,erpFigDir)
%% EEG topo of onset across coherence ON AVERAGE, NOT DIFF
doTopos({'EEG'},eegNeighbours,ecOnsAll,hcOnsAll,...
    zlimEegOns,0.05,[-0.3,1.5],[6,6],36,...
    cOnsAve,eegLayout,'eeg_coh_ons_on_ave_sig_ERP_clusters','eeg_coh_ons_sig_ERP_clusters',datadir,erpFigDir)
%% EEG topo of response across coherence ON AVERAGE, NOT DIFF
doTopos({'EEG'},eegNeighbours,ecRespAll,hcRespAll,...
    zlimEegResp,0.05,[-0.6,0.2],[4,4],16,...
    cRespAve,eegLayout,'eeg_coh_resp_on_ave_sig_ERP_clusters','eeg_resp_ons_sig_ERP_clusters',datadir,erpFigDir)
%% let's plot that frontal stuff: coherence
cOnsDiffAve.bfs = testDiffsAcrossTime(cOnsDiffAll,frontal);
plotTimecourse(eegLayout,[teal;coral],... % layout and colours
    {ecOnsAve,hcOnsAve},cOnsDiffAve,... % {average,data},averageDifference (with bayes factors)
    {ecOnsAll,hcOnsAll},... % {subjectwise,data}
    [-0.2 1.5],[],frontal,... % xlims,ylims,channels
    {'EEG','uV'},... % sensor type and units for ylablel
    {'easy coh';'hard coh'},'northwest','Frontal sensors in EEG: Coherence Onset',... % legend, legend location, plot title
    [erpFigDir filesep 'eeg_coh_ons_frontal_ERP.png']) % save loc
cRespDiffAve.bfs = testDiffsAcrossTime(cRespDiffAll,frontal);
plotTimecourse(eegLayout,[teal;coral],... % layout and colours
    {ecRespAve,hcRespAve},cRespDiffAve,... % {average,data},averageDifference (with bayes factors)
    {ecRespAll,hcRespAll},... % {subjectwise,data}
    [],[],frontal,... % xlims,ylims,channels
    {'EEG','uV'},... % sensor type and units for ylablel
    {'easy coh';'hard coh'},'northwest','Frontal sensors in EEG: Coherence Response',... % legend, legend location, plot title
    [erpFigDir filesep 'eeg_coh_resp_frontal_ERP.png']) % save loc
%% let's plot that frontal stuff: categorisation
rOnsDiffAve.bfs = testDiffsAcrossTime(rOnsDiffAll,frontal);
plotTimecourse(eegLayout,[teal;coral],... % layout and colours
    {erOnsAve,hrOnsAve},rOnsDiffAve,... % {average,data},averageDifference (with bayes factors)
    {erOnsAll,hrOnsAll},... % {subjectwise,data}
    [-0.2 1.5],[],frontal,... % xlims,ylims,channels
    {'EEG','uV'},... % sensor type and units for ylablel
    {'easy cat';'hard cat'},'northwest','Frontal sensors in EEG: Categorisation Onset',... % legend, legend location, plot title
    [erpFigDir filesep 'eeg_cat_ons_frontal_ERP.png']) % save loc
rRespDiffAve.bfs = testDiffsAcrossTime(rRespDiffAll,frontal);
plotTimecourse(eegLayout,[teal;coral],... % layout and colours
    {erRespAve,hrRespAve},rRespDiffAve,... % {average,data},averageDifference (with bayes factors)
    {erRespAll,hrRespAll},... % {subjectwise,data}
    [],[],frontal,... % xlims,ylims,channels
    {'EEG','uV'},... % sensor type and units for ylablel
    {'easy cat';'hard cat'},'northwest','Frontal sensors in EEG: Categorisation Response',... % legend, legend location, plot title
    [erpFigDir filesep 'eeg_cat_resp_frontal_ERP.png']) % save loc
%% let's plot a suss sensor: coherence
cOnsDiffAve.bfs = testDiffsAcrossTime(cOnsDiffAll,'EEG061');
plotTimecourse(eegLayout,[teal;coral],... % layout and colours
    {ecOnsAve,hcOnsAve},cOnsDiffAve,... % {average,data},averageDifference (with bayes factors)
    {ecOnsAll,hcOnsAll},... % {subjectwise,data}
    [-0.2 1.5],[],'EEG061',... % xlims,ylims,channels
    {'EEG','uV'},... % sensor type and units for ylablel
    {'easy coh';'hard coh'},'northwest','Iz in EEG: Coherence Onset',... % legend, legend location, plot title
    [erpFigDir filesep 'eeg_coh_ons_Suss.png']) % save loc

%% MEG topo of onset across coherence
doTopos({'MEG'},megNeighbours,ecOnsAll,hcOnsAll,...
    zlimMegOns,0.1,[-0.1,1.5],[4,4],16,...
    cOnsDiffAve,megLayout,'meg_coh_ons_on_diff_sig_ERP_clusters','meg_coh_ons_sig_ERP_clusters',datadir,erpFigDir)
%% MEG topo of response across coherence
doTopos({'MEG'},megNeighbours,ecRespAll,hcRespAll,...
    zlimMegResp,0.05,[-0.6,0.2],[4,4],16,...
    cRespDiffAve,megLayout,'meg_coh_resp_on_diff_sig_ERP_clusters','meg_coh_resp_sig_ERP_clusters',datadir,erpFigDir)
%% MEG topo of onset across categorisation
doTopos({'MEG'},megNeighbours,erOnsAll,hrOnsAll,...
    zlimMegOns,0.1,[-0.1,1.5],[4,4],16,...
    rOnsDiffAve,megLayout,'meg_cat_ons_on_diff_sig_ERP_clusters','meg_cat_ons_sig_ERP_clusters',datadir,erpFigDir)
%% MEG topo of response across categorisation
doTopos({'MEG'},megNeighbours,erRespAll,hrRespAll,...
    zlimMegResp,0.05,[-0.6,0.2],[4,4],16,...
    rRespDiffAve,megLayout,'meg_cat_resp_on_diff_sig_ERP_clusters','meg_cat_resp_sig_ERP_clusters',datadir,erpFigDir)


%% subfunctions

function difference = getDifference(data1,data2)

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
difference = ft_math(cfg, data1, data2);

return
end

function figHandle = plotSensors(structure,channelRow,time)
% structure should be a fieldtrip data structure
% channelRow should be the row number of the channel you want
% time should be the range of timepoints you want

selected_data = structure.avg(channelRow,time);
selected_time = structure.time(time);
figure;
figHandle = plot(selected_time, selected_data)

return
end

function bfs = testDiffsAcrossTime(dataStruct,channels,complementary,nullInterval)

if ~exist('nullInterval','var'); nullInterval = '0.5,Inf'; end

for i = 1:numel(dataStruct)
    % if multiple channels, get the average, else get the channel
    if size(dataStruct{i}.avg(find(ismember(dataStruct{i}.label,channels)),:),1) > 1
        diffs(:,i) = mean(dataStruct{i}.avg(find(ismember(dataStruct{i}.label,channels)),:));
    else
        diffs(:,i) = dataStruct{i}.avg(find(ismember(dataStruct{i}.label,channels)),:);
    end
end; clear i

% add the R module and get the path to Rscript
[status, result] = system('module add R && which Rscript');
if status == 0
    disp('R module added successfully');
    RscriptPath = strtrim(result);
    disp(['Rscript path: ' RscriptPath]);
else
    error('Failed to add R and/or locate Rscript');
end

% now run the rscript version of the bayes analysis
%   we can also get the bf for the complementary interval
%   by specifying complementary = 2. Let's set a default:
if ~exist('complementary','var'); complementary = 1; end
bfs = bayesfactor_R_wrapper(diffs,'Rpath',RscriptPath,'returnindex',complementary,...
    'args',['mu=0,rscale="medium",nullInterval=c(' nullInterval ')']);

% alternatively they have implemented it in matlab
% [bfs, bfs_complementary_interval] = bayesfactor(diffs, 'interval',[-Inf Inf]);

return
end

function [pvals, fdrh, stats] = doAnova(thisAnova)
% so we'll get six results from this:
% main effect of subjects
% main effect of difficulty
% main effect of manipulations
% interaction between subjects and diff
% interaction between subjects and manip
% interaction of diff and manip
variableNames = {'subjects' 'difficulty' 'manipulation'};
initResults = @(x) NaN(size(x,2), 6); % function to init a matrix per timepoint for six results
pvals = initResults(thisAnova); 
fstat = initResults(thisAnova); 
meanSq = initResults(thisAnova); 
subjs = repmat((1:size(thisAnova,1))',1,4); % create subject ids: a row with a unique number for each subject and then replicated horizontally four times for the four conditions
cohLevel = repmat([1 1 2 2],size(thisAnova,1),1); % repmat the levels of coherence for the number of subjects: [easy easy hard hard]
ruleLevel = repmat([1 2 1 2],size(thisAnova,1),1); % repmat the levels of rule for the number of subjects: [easy hard easy hard]
% anovan then combines the difficulty and manipulation vectors to create the four distinct groups:
% 1,1: easy coherence, easy rule
% 1,2: easy coherence, hard rule
% 2,1: hard coherence, easy rule
% 2,2: hard coherence, hard rule
for timepoint = 1:size(thisAnova,2) % loop through timepoints
    if timepoint>1; fprintf(repmat('\b', 1, numel(timepointStr))); end
    timepointStr = sprintf('timepoint %.0f of %.0f', timepoint, size(thisAnova,2));
    fprintf(timepointStr);
    tthisAnova = squeeze(thisAnova(:,timepoint,:)); % now we extract a 3D slice: all subjects, all conditions for a single timepoint
    % we squeeze() it so it makes it 2D (because the 3rd dimension is not
    % necessary as timepoint is only 1
    [pvals(timepoint,:), t] = anovan(tthisAnova(:),{subjs(:),cohLevel(:),ruleLevel(:)},'random',1,'display','off','model','interaction','varnames',variableNames);
    for tidx = 1:6 % loop through our six results to compile our test statistics for each timepoint
        % t is basically a table: statistic names across the top, and test equations down the left side
        fstat(timepoint,tidx) = t{tidx+1,6}; % grab the fstat
        meanSq(timepoint,tidx) = t{tidx+1,5}; % grab the mean square
    end; clear tidx

end; clear timepoint

% now Niko Kriegeskorte's FDR correction for that set of p-values
% if we want to correct each test seperately
fdrp = NaN(1,6);
for pidx = 1:6
    if ~isempty(FDRthreshold(pvals(:,pidx)))
        fdrp(:,pidx) = FDRthreshold(pvals(:,pidx));
    end
end; clear pidx
% if we want to correct across tests instead
% fdrp = FDRthreshold(pvals);
% now let's get where the null is rejected
fdrh = pvals < fdrp;

% let's make a table of these, so we can see what the cols correspond to
testNames = t(2:7,1)'; testNames = cellfun(@(x) strrep(x,'*','_x_'),testNames,'UniformOutput',false);
pvals = array2table(round(pvals,3), 'VariableNames', testNames);
fdrh = array2table(fdrh, 'VariableNames', testNames);

% and compile any other useful stats
stats.fstat = array2table(fstat, 'VariableNames', testNames);
stats.meanSq = array2table(meanSq, 'VariableNames', testNames);
stats.fdrp = fdrp;

return
end



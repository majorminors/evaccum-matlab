function fit_MEG_LBA_LR_lag_sourcePW_phaseperm(settings)

%estimates statistical dependency (i.e. correlation) between already smoothed
%envelopes and simulated signals
%Also estimates MI between simulated signals and baseline useful for
%cluster statistics

rng(17,'twister') % for reproducibility

load(settings.ModData);
load(settings.ModRT);
load(settings.dataRT,'trialinfo');

bands = settings.bandn;
%[datatmp,dInfo] =prepareVertices(settings.data,settings.band,settings.invnum,settings.orth);
[datatmp,dInfo] =prepareVertices(settings.data,settings.band,settings.invnum,settings.orth);

sprintf('Fitting %s at %d-%d Hz',settings.data,settings.band)
%check trial consistency
dInfo.behav(isnan(dInfo.behav(:,3)),3)=0;
if sum((trialinfo(:,3)-dInfo.behav(:,2)));error(['wrong trial mapping:', settings.data]);end%#ok



%Bw = settings.Bw;
datatmp = single(datatmp);%reduce memory requirements
datatmp = permute(datatmp,[3 1 2]); %tr x ch x time

subj   = settings.subNum;
%Lockid = settings.Lock;
design = settings.design;
lBand  = settings.lBand;
fprintf('processing subj # %d \n',subj)



if ~exist(settings.datafd,'dir');mkdir(settings.datafd);end

%Get parameters
[~,pidx]=min(bestval{subj});%#ok %best fit 
[~,params{1},params{2},params{3},params{4}]=getModelParam_cell_RDK(design,4,bestpar{subj}(pidx,:));%#ok


%Get good RTs; NOTE: in RT.data bad trials have already been removed!
bad = trialinfo(:,2)==0;

data.time = dInfo.time;
rts=[];for i = 1:length(RT{subj}.data);rts = cat(1,rts,RT{subj}.data{i});end;%#ok
rts = sortrows(rts,1);

%remove RTs shorter than T0
tooshort = find(trialinfo(:,1)<=params{1}.T0); 
bad(tooshort) = 1;

if ~isempty(tooshort); rts(ismember(rts(:,1), tooshort),:)=[];end

rts(:,5) = trialinfo(~bad,3);
data.rts = rts;



data.locki   = settings.locki;
data.fsample = 1/dInfo.fsample;


%construct regressors- for each trial only the ramping activity is
%simulated. The regressor is built by concatenating all trials. This allows
%for fitting all trials at once since we expect lags to differ only across
%ROIs

[~,signals,times,~]=buildReg(data,params);

% if settings.surrogate%build surrogate distribution
%    load(settings.lagsf)
%    idx = randperm(length(signals));
%    for is = 1:length(signals)
%        tmp1{is} = signals{idx(is)};
%        tmp2{is} = times{idx(is)};
%    end
%    signals = tmp1;
%    times   = tmp2;
%    clear tmp1 tmp2;
% end

npoints = cellfun(@(x)(length(x)),signals);

%% fit each band separately 

%check if previous sessions need to be recovered
%bands = 1:length(lBand);  

if exist(settings.outf,'file') && ~settings.overwrite;load(settings.outf,'logband'); 
    %bands = bands(~ismember(bands,logband)) ;%#ok
    sprintf('Recovering %s band %s', settings.outf,lBand{bands(1)})
end      



    data.meg = datatmp(~bad,:,:);%tr ch time
    clear datatmp;

    %% fit each ROI separately 
%MI = zeros(size(data.meg,2),numel(lBand));
clear tmp;
for chan_n = 1:size(data.meg,2)
    
    
    data.chan_n = chan_n;
    data.T0     = params{1}.T0;
    data.T0s    = floor(data.T0/data.fsample);%T0 in samples
    data.lagms  = 10/1000;%lag resolution in ms
    startlags   = 0.01/data.fsample;%(samples) perceptual processes take longer than 10ms
    data.lag    = data.lagms/data.fsample;%lag resolution in samples
    data.lags   = startlags+1:data.lag:data.T0s;%lags in samples
%     if settings.surrogate
%         data.lags = lags(chan_n,band);
%     end
    %[b(band).c(chan_n).R, b(band).c(chan_n).LR] = LR_trials_s(data,signals,npoints,times,settings.smoothtype);%#ok
     tmp.c(chan_n).R = LR_trials_s(data,signals,npoints,times,settings.smoothtype);% 235.18 seconds on the CBU's cluster
      
    
end



if exist(settings.outf,'file') ; load(settings.outf,'b','logband');
   b(bands) = tmp;%#ok
   logband(bands) = bands;
   save(settings.outf,'b','logband','-append');
else
    Readme.modeltype = {'1-only ramping component;','2-full spline model','3-full spike model'};
    Readme.fittype   = {'Spearman Rho ','Linear Regression','analysis trial-by-trial(loc) and pooled'};
    Readme.smoothtype= {'1 Median filter win=51','2  win = 10','3 no smooth'};
    Readme.timestamp = date;% use datenum to transform date into the number of days from 1/1/0000
    b(bands) = tmp;%#ok
    conds = [data.rts(:,1),data.rts(:,5)];%#ok
    data = rmfield(data,'meg');
    logband(bands) = bands;%#o
    save(settings.outf,'b','settings','Readme','conds','logband','data'); 
end
clear b;


function R = LR_trials_s(data,regressor1,npoints,timesT,smoothtype)
% Spearman's RHO on single trial
totTrials = size(data.rts,1);
time = data.time;
data2fit = permute(data.meg(:,data.chan_n,:),[1,3,2]);%trial x time




%% Prepare data 
%more elegant but 4 times slower! consider using for github version.
t=zeros(size(timesT,2),1);%uses times just to double check number of trials
% tic
% t = cellfun('length',times); %RT in samples
% Nt = numel(t);
% 
% if data.locki == 2
%     onoff{2} = dsearchn(time',0).*ones(1,Nt);
%     onoff{1} = onoff{2}-t+1;    
%     
% else
%     onoff{1} = dsearchn(time',0).*ones(1,Nt);
%     onoff{2} = onoff{1}+t-1;%from 0 to RT
%        
% end
% 
% signal = arrayfun(@(x,y,z) (data2fit(x,y:z)'),1:Nt,onoff{1},onoff{2},'UniformOutput',0);
% vecLen = [zeros(1,Nt); npoints-1; (onoff{2}-onoff{1})+1]';
% toc




for nTrial = 1:totTrials
    
    t(nTrial) = size(timesT{nTrial},2);
    
if data.locki == 2
    onoff(2) = dsearchn(time',0);
    onoff(1) = onoff(2)-t(nTrial)+1;    
    
else
    onoff(1) = dsearchn(time',0);
    onoff(2) = onoff(1)+t(nTrial)-1;%from 0 to RT
       
end

   
   if smoothtype == 1
    datatmp = smoothdata(data2fit(nTrial,:),'movmedian',51);
   elseif smoothtype  == 2
    datatmp = smoothdata(data2fit(nTrial,:),'movmedian',11);        
       
   else
       datatmp = data2fit(nTrial,:);%default option is no smoothing        

   end


signal{nTrial}   = datatmp(onoff(1):onoff(2))';%#ok data
vecLen(nTrial,:) = [0 npoints(nTrial)-1 numel(signal{nTrial})];%#ok
end
clear datatmp data2fit timesT;


%% Generate lagged regressors
TotLength = sum(vecLen(:,3),1);
TotLags   = length(data.lags);
regressor1 = cellfun(@transpose,regressor1,'UniformOutput',0);% 1 x time
regressor2 = zeros(TotLength,TotLags);

for nLag = 1:TotLags
        
        %shifts the signal according to lag: data points are taken from the time
        %window corresponding to the ramping part at different lags.        
        lags = num2cell([vecLen(:,1:2)+data.lags(nLag).*ones(size(vecLen,1),2) vecLen(:,3)],2)';%[lag veclen+lag length epoch]                
        %Full shape model activation (e.g. spline or spike) is assembled here
        %by adding 0 and 1 at the start/end tails
        tmpreg=cellfun(@(x,y) ([zeros(1,y(1)-1),x',ones(1,y(3)-y(2))*x(end)]'),regressor1,lags,'UniformOutput',0);%create full model curve
        regressor2(:,nLag) = cat(1,tmpreg{:});clear tmpreg;
        
end


%% Generate the null distribution of correlations

NPerm = 10000;%per lag
Msignal = cat(1,signal{:});
%shuffle trials
rsign = zeros(TotLength,NPerm);
nSamples = TotLength;
nVars = 1;
mean_term = mean(Msignal);
std_term = std(Msignal,1); % Use N-1 standard deviation to match original code


% Generate surrogate data though method 3 in Hurtado et al 2004. Statistical
% method for detection of phase locking episodes in neural oscillations. J
% Neurophysiol. 10.1152/jn.00853.2003.
%
% This scrambles the phase spectrum whilst preserving the amplitude spectrum
% The surrogate data is rescaled to the same mean and standard deviation

for i = 1:NPerm    
    
% get amplitude spectrum
amp_spec = fft(Msignal, nSamples,1);      % no need to take abs as we can just add random phases
% Generate phase noise
% first element of spectrum is DC component, subsequent elements mirror
% each other about the Nyquist frequency, if present
n_components = floor((nSamples-1)/2); % Number of components that *aren't* DC or Nyquist
noise = rand(n_components, nVars) .* (2 * pi);
if rem(nSamples,2) % If odd number of samples, then Nyquist frequency is NOT present
    newPhase = [zeros(1,nVars); 1i .* noise; -1i .* flipud(noise)]; % second half uses conjugate phases, mirrored in order
else % Otherwise, include zero phase shift for the Nyquist frequency
    newPhase = [zeros(1,nVars); 1i .* noise; zeros(1,nVars); -1i .* flipud(noise)]; % second half uses conjugate phases, mirrored in order
end
% Make new phase scrabled spectrum
rand_spec = exp(newPhase).*amp_spec;

% Create new time_course
surrogate = ifft(rand_spec, nSamples);
% demean
surrogate = bsxfun(@minus, surrogate, mean(surrogate));

% Normalise time_series
surrogate = bsxfun(@times, surrogate, std_term./std(surrogate,1));  
% Reset the mean
surrogate = bsxfun(@plus, surrogate, mean_term);

   rsign(:,i) = surrogate;
end
clear Msignal surrogate noise;


%rank data separately to speed up
regressor2 = tiedrank(regressor2);
rsign      = tiedrank(rsign);
   
%calculate Spearman's correlation
tmpr = zeros(nLag,NPerm);
for i = 1:nLag
    tmpr(i,:) = corr(regressor2(:,i),rsign);%Pearson's R on ranked data == Spearman's Rho but faster.
end    
    
R.surrogate = single(tmpr);


%% Calculate statistics for the actual time series
R.Rho = corr(regressor2,cat(1,signal{:}),'type','Spearman');




function [regressor,signals,times,conds]=buildReg(data,params)

regressor = [];
totTrials = size(data.rts,1);

for trial_num = 1:totTrials
    rt_i = data.rts(trial_num,2);  finger = data.rts(trial_num,3); avoid = data.rts(trial_num,4);cond = data.rts(trial_num,5);
    
    n_units = 4;
    if avoid == 0; n_units = 1;end    
    %theoretical signal
    [signals{trial_num},times{trial_num}]=EA_LBA(round(params{1}.T0/2,2),params{cond},rt_i,n_units,finger,avoid,data.fsample,3);%#ok 1spike 2spline   3only ramping
    regressor = [regressor signals{trial_num}];%#ok
    
end

conds = [data.rts(:,1) data.rts(:,5)];
function CondEntropy


clear all
addpath('/imaging/at07/Matlab/CommonScripts/exportFig');

%% --------------------
sname = {'BF','CM','DS','EK','EW','IF','IM','JC','JS','LF','PK','SG','SS','TG'};%JR GA ms

stemf = '_MEGBehav.mat';

%initialise cells
Ava_opt = {}; dataRaw = {}; CohLevs = {};
 figure(3)
 set(gcf,'Position',[3 -21 1596 1059]);

plotii = 0; 
for idxsubj = 1:length(sname)
    idxsubj
    %[foo1{subj,1},foo2{subj,1},dataRaw{subj,1}]=xlsread(sprintf('../data/%d.xls', sname(subj)));
    %% -----------------
    clear ResponseArray LocStim 
    
    load([sname{idxsubj},stemf])
   
    %% for this 1st VERSION let's ignore the errors ~15% however CHECK HOW TO USE THEM
idxCorr =ResponseArray(:,end)==1;
%  if ~isempty(idxErr)   
%      ResponseArray(idxErr,end)=[];
%      LocStim(idxErr,:)=[];
%  end

    %extract unavailable options;
    Ava_opt{idxsubj} = LocStim;
    dataRaw{idxsubj} = ResponseArray(:,[3 2 14 13]);
    CohLevs{idxsubj} = conf.sesInfo.cohSet_feature;
    
    idxE = 0;
for icL = [2 4]   

    %% Calculate conditional entropy for repetition choice trials 
    idxE = idxE+1;
    %select coh level
    tmpData = dataRaw{idxsubj}(dataRaw{idxsubj}(:,2)==CohLevs{idxsubj}(icL),:);
    
     %get rid of too late (inf) or wrong (nan) responses
    tmpData((isnan(tmpData(:,3)) |isinf(tmpData(:,3))),:)=[];
    
    
    %find choice trials    
    t1 = find(tmpData(:,1)==3);
    %t1 = t;
    t0 = t1-1;
    
   
    
    
    
    if ~isempty(t0<=0); t1(t0<=0)=[]; t0(t0<=0)=[];end
        
    %consider only  repetition allowed trials
    idxt=[];
    for aa = 1:numel(t1); idxt(aa)=LocStim(t1(aa),tmpData(t0(aa),3))==1;end %loc stim at t1 needs to include choice at t0
    
    C_Entropy{idxsubj}(idxE)=ConditionalEntropy(tmpData(t1(logical(idxt)),3),tmpData(t0(logical(idxt)),3));

    plotii = plotii+1; 
    %% Calculate transition probabilities : may want to add behav training session/collapse across subjs
    for finger = 1:4
        
    t1 = find(tmpData(:,1)==3 & tmpData(:,3)==finger);
    t0 = t1-1;
    t2 = t1(1:end-1)+1;
    if ~isempty(t0<=0); t1(t0<=0)=[]; t0(t0<=0)=[];end
    
    for finger2 = 1:4; Pchoice(finger,finger2,idxE)=mean(tmpData(t0,3) == finger2);end
    
    
    subplot(5, 4,plotii)
    imagesc(Pchoice(:,:,idxE));
    caxis([0,0.8]);
    axis square
    colorbar
    set(gca, 'XTick',[1:4],'XTickLabel',{'index','middle','ring','little'},'FontSize',9)
    set(gca, 'YTick',[1:4],'YTickLabel',{'index','middle','ring','little'},'FontSize',9)
    title(['CondEntropy:',num2str(C_Entropy{idxsubj}(idxE))],'FontWeight','bold','FontSize',12)
    set(gca,'LineWidth',2)
    set(gcf,'Color','white')
    xlabel('CHOICE T_1','FontWeight','bold','FontSize',10);
    ylabel('CHOICE T_0','FontWeight','bold','FontSize',10);
    text(-1.5,0,['S',num2str(idxsubj)],'FontSize',15,'Color','white','BackgroundColor','black')
    
    
    end
PIndChoice{idxsubj} = Pchoice;

end   
 
if plotii == 20 || idxsubj == numel(sname)
    plotii = 0;
    export_fig(['/imaging/at07/Matlab/Projects/CBU2015/RDKUnc/Model/Cond_Entropy.pdf'],'-pdf','-append')
    clf
end



    
end


% ConditionalEntropy: Calculates conditional entropy (in bits) of Y, given X
% by Will Dwinnell
%
% H = ConditionalEntropy(Y,X)
%
% H  = calculated entropy of Y, given X (in bits)
% Y  = dependent variable (column vector)
% X  = independent variable(s)
% 
% Note: requires 'Entropy' and 'MutualInformation' functions
%
% Example: 'X' (1 bit of entropy), 'Y' (2 bits of entropy)
%   Given 'X', 'Y's conditional entropy is 1 bit.
%
% Note: Estimated entropy values are slightly less than true, due to
% finite sample size.
%
% Y = ceil(4 * rand(1e3,1));  X = double(Y <= 2);
% ConditionalEntropy(Y,X)
%
% Last modified: Nov-12-2006

function H = ConditionalEntropy(Y,X)

% Axiom of information theory
H = Entropy(Y) - MutualInformation(X,Y);


% God bless Claude Shannon.

% EOF

% MutualInformation: returns mutual information (in bits) of the 'X' and 'Y'
% by Will Dwinnell
%
% I = MutualInformation(X,Y);
%
% I  = calculated mutual information (in bits)
% X  = variable(s) to be analyzed (column vector)
% Y  = variable to be analyzed (column vector)
%
% Note: Multiple variables may be handled jointly as columns in matrix 'X'.
% Note: Requires the 'Entropy' and 'JointEntropy' functions.
%
% Last modified: Nov-12-2006

function I = MutualInformation(X,Y)

if (size(X,2) > 1)  % More than one predictor?
    % Axiom of information theory
    I = JointEntropy(X) + Entropy(Y) - JointEntropy([X Y]);
else
    % Axiom of information theory
    I = Entropy(X) + Entropy(Y) - JointEntropy([X Y]);
end

% Entropy: Returns entropy (in bits) of each column of 'X'
% by Will Dwinnell
%
% H = Entropy(X)
%
% H = row vector of calculated entropies (in bits)
% X = data to be analyzed
%
% Example: Measure sample entropy of observations of variables with
%   1, 2, 3 and 4 bits of entropy.
%
% Note: Estimated entropy values are slightly less than true, due to
% finite sample size.
%
% X = ceil(repmat([2 4 8 16],[1e3,1]) .* rand(1e3,4));
% Entropy(X)
%
% Last modified: Nov-12-2006

function H = Entropy(X)

% Establish size of data
[n m] = size(X);

% Housekeeping
H = zeros(1,m);

for Column = 1:m,
    % Assemble observed alphabet
    Alphabet = unique(X(:,Column));
	
    % Housekeeping
    Frequency = zeros(size(Alphabet));
	
    % Calculate sample frequencies
    for symbol = 1:length(Alphabet)
        Frequency(symbol) = sum(X(:,Column) == Alphabet(symbol));
    end
	
    % Calculate sample class probabilities
    P = Frequency / sum(Frequency);
	
    % Calculate entropy in bits
    % Note: floating point underflow is never an issue since we are
    %   dealing only with the observed alphabet
    H(Column) = -sum(P .* log2(P));
end

% JointEntropy: Returns joint entropy (in bits) of each column of 'X'
% by Will Dwinnell
%
% H = JointEntropy(X)
%
% H = calculated joint entropy (in bits)
% X = data to be analyzed
%
% Last modified: Aug-29-2006

function H = JointEntropy(X)

% Sort to get identical records together
X = sortrows(X);

% Find elemental differences from predecessors
DeltaRow = (X(2:end,:) ~= X(1:end-1,:));

% Summarize by record
Delta = [1; any(DeltaRow')'];

% Generate vector symbol indices
VectorX = cumsum(Delta);

% Calculate entropy the usual way on the vector symbols
H = Entropy(VectorX);


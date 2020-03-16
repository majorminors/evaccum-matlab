function bestpar=fitparams_template(errfunc,startpar,mod_feature,data2fit,paramconstrains)
% Fitting parameters of a computational model
% Required inputs:
%  model - string containing a name of Matlab function implementing the connectionist model,
%          it should get as input the vector of parameters and return the values of statistics, i.e.
%          function statistics = model (parameters)
%  rt_data - vector for experiment response time
%  startpar - starting values of parameters
%  randiter - when = 0 -> the provided parameters are used as starting point for optimization
%             when > 0 -> number of itearions of random search for starting point (range of search space
%             is [0, 2*startpar], i.e. around the starting point; default = 0
%  paramconstrains - allowed range of parameters, matrix containing 2 columns with length = no of parameters
%	     the first column contains minimal, the second maximal values of params; default = no contraints
%
% Outputs:
%  bestpar - best sets of parameters found during each session
%  bestval - error for the corresponding set of parameters
%  bestat - values of the statistics for the corresponding set of parameters
%  p - significance of difference between model and experimental statistics

%  goalstat - the values of the statistics you want to achieve
%
%  Written by Jiaxiang Zhang 10/08/2010

if length(data2fit)==1
    if isstruct(data2fit)
        goalstat=data2fit;
    else
        goalstat=data_stats(data2fit);
    end
else
    if isstruct(data2fit{1})
        goalstat=data2fit;
    else
        for j=1:length(data2fit)
            goalstat{j} = data_stats(data2fit{j});
        end
    end
end

%       allObs: observations across responses
%       indvObs: observations for individual response
%       priorProb: probability for individual response
%  for allObs and indvObs:
%        column 1 - quantils
%        column 2 - quantils boundaries
%        column 3 - the probability mass in each bin
%        column 4 - the observed frequencies in each bin


% prepare output variables
nopara = length (startpar); % num of parameters
goalen = length (goalstat); % num of goal statistics
bestpar = zeros (1, nopara);
bestval = ones (1, 1) * 10000;
bestat = zeros (1, goalen);

% handling starting parameters equal to 0
startscaled = ones (1, nopara);
for i = 1:nopara
    if startpar(i) == 0
        startpar(i) = mean(parange(i,:));
        if startpar(i) == 0
            startpar(i) = (parange(i,1) + 3*parange(i,2)) / 4;
        end
        startscaled(i) = 0;
    end
end

% startpar(startpar==0)=eps;
bestpar=startpar;

% do the optimization
% options = zeros(1, 18);
% options (5) = 1;    % use Nelder-Mead Simplex
% options (14) = 100;
% options (18) = 0.3;
% disp ('Optimizing parameters');
% bestpar (:) = subplex (errfunc, bestpar(:)', options, [], ...
%     model,goalstat);
options=optimset('fminsearch');
options.Display='iter';
options.MaxFunEvals=10000;
% bestpar=fminsearchbnd(errfunc,bestpar(:)',zeros(1,length(startpar)),zeros(1,length(startpar))+10,options,mod_feature,goalstat);
bestpar=fminsearchbnd(errfunc,bestpar(:)',paramconstrains(:,1)',paramconstrains(:,2)',options,mod_feature,goalstat);

end






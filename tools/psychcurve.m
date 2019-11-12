function [result,hline,hdata] = psychcurve(data)
% data [intensity trialscorrect tottrials]


if nargin == 0 %run demo
    data =    [...
    0.10,   45.0000,   90.0000;...
    0.15,   50.0000,   90.0000;...
    0.20,   44.0000,   90.0000;...
    0.25,   44.0000,   90.0000;...
    0.30,   52.0000,   90.0000;...
    0.35,   53.0000,   90.0000;...
    0.40,   62.0000,   90.0000;...
    0.45,   64.0000,   90.0000;...
    0.50,   76.0000,   90.0000;...
    0.60,   79.0000,   90.0000;...
    0.70,   88.0000,   90.0000;...
    0.80,   90.0000,   90.0000;...
    1.00,   90.0000,   90.0000];

end


% To start psignifit you need to pass a struct, which specifies, what kind
% of experiment you did and any other parameters of the fit you might want
% to set:

% You can create a struct by simply calling [name]=struct

options             = struct;   % initialize as an empty struct

%Now you can set the different options with lines of the form
%[name].[field] as in the following lines:

options.sigmoidName = 'norm';   % choose a cumulative Gaussian as the sigmoid
options.expType     = '2AFC';   % choose 2-AFC as the paradigm of the experiment

%% Now run psignifit
% Now we are ready to run the main function, which fits the function to the
% data. You obtain a struct, which contains all the information about the
% fitted function and can be passed to the many other functions in this
% toolbox, to further process the results.

result = psignifit(data,options);

%result is a struct which contains all information obtained from fitting your data. 
%Perhaps of primary interest are the fit and the confidence intervals:

result.Fit
result.conf_Intervals

% This gives you the basic result of your fit. The five values reported are:
%    the threshold
%    the width (difference between the 95 and the 5 percent point of the unscaled sigmoid)
%    lambda, the upper asymptote/lapse rate
%    gamma, the lower asymptote/guess rate
%    eta, scaling the extra variance introduced (a value near zero indicates 
%         your data to be basically binomially distributed, whereas values 
%         near one indicate severely overdispersed data)
% The field conf_Intervals returns credible intervals for the values provided 
% in options.confP. By default these are 68%, 90% and 95%. With default settings 
% you should thus receive a 5x2x3 array, which contains 3 sets of credible intervals 
% (lower and upper end = 2 values) for each of the 5 parameters.

%% visualize the results
% For example you can use the result struct res to plot your psychometric
% function with the data:

[hline,hdata]=plotPsych(result);


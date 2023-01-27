function sendToCluster(functionToRun, inputArgs, jobpaths)
% functionToRun must be function handle, e.g.
%   functionToRun = @yourScript;
%   functionToRun = str2func(yourScript);
% jobpaths must be a cell array of path strings: whatever paths it might
% need to do it's job (although you can also just 'addpath()' in the script 
% inputArgs must be cell array of whatever input arguments your script needs 
% some notes below, but details can be found at http://intranet.mrc-cbu.cam.ac.uk/home/cluster-matlab/

disp('assembling job')
J.task=functionToRun; % this is a function handle to the script that will run the job
J.n_return_values=0; % this is an integer that represents the number of values returned by the job script
J.input_args={inputArgs{:}}; % this is a cell array with one element for each input argument to be passed to the job script
% this is vector of numbers representing the jobs that must complete before running the current job
% jobs are represented by their index number in the job array, so you could
% in theory, I guess, run a series of consecutive jobs? but the example the cbu cluster intranet page has
% doesn't use this, so I don't know what it's for really
J.depends_on=0; 

workers = 1; % how many workers
gbPerWorker = 4; % how many GB of RAM
secondsOfWalltime = 3600; % how many seconds to allow it to run (14400 secs = 4 hours)

disp('creating scheduler object')
clear S
S = cbu_scheduler('custom',{'compute',workers,gbPerWorker,secondsOfWalltime});

% S.jobStorageLocation = '/path/to/specific/directory'
% default will be '/imaging/<user name>/.cbu-cluster/matlab-jobs'

% this stops it traversing the directory structure of all the functions in 
% scripts, which slows things down if you add to the path a library of
% functions for example
jobopts.AutoAttachFiles=false;

disp('submitting job')
cbu_qsub(J,S,jobpaths,jobopts);

disp('done')

end
function fun_fitbehav_LBA_PD(settings)
% Fit LBA model to RDK-tapping task
% Use for  3-choice task: CONDITIONS segregated
% JZ 15/10/2014 - 4-choice tapping task
% AT 17/02/2016 adapted to 3-choice RDK task
% DM 02/2020 adapted to 2-choice RDK task

rseed = settings.rseed;
rng(rseed,'twister') % for reproducibility

%clear all
%%
savename = settings.savename;
data = settings.data;
randiter  = settings.randiter;
nosession = settings.nosession;
Model_Feature = settings.modfeat; % model variant for this job
numParam   = getModelParam_cell_RDK(Model_Feature,2); % the 2 here reflects response options
parange=[zeros(1,numParam);zeros(1,numParam)+10]';

% optional but doesn't save time (thought it would be faster than normal
% gradient descent, but isn't)
if ~isfield(settings,'bayesOptim')
    bayesian = 0;
else
    bayesian = settings.bayesOptim;
    if bayesian
        addpath(genpath('/imaging/at07/Matlab/Toolboxes/bads'));
    end
end

num_subjs = sum(~cellfun('isempty',{data.id})); % get the number of not empty arrays from the first field of the 'data' structure

%initialise cells
data2fit={};

for idxsubj = 1:num_subjs
    %%
    clear rt conds resp acc % clear vars from behavioural data to avoid reading mixed subject data
    fprintf('loading subject %1.0f \n',data(idxsubj).id);
    
    % pull the data
    subjdata = data(idxsubj).data; % put the data in a readable variable
    dataValid = []; % now we'll strip invalid responses out
    validLabs = {};
    for i = 1:size(subjdata,1)
        if subjdata{i,5} >= 0 % if there's a valid response
            dataValid(end+1,:) = [subjdata{i,2} subjdata{i,3} subjdata{i,4} subjdata{i,5}]; % add the following rows to this new variable in order: condition code, response, rt, accuracy
            temp = subjdata{i,1}; % extract valid condition label (can't do in one step with multi-level structures and non-scalar indexing)
            validLabs{end+1} = temp; % add the valid condition label
        end
    end; clear i temp;
    
    % converting again to readable variables
    conds = dataValid(:,1);
    easyOrHard = [];
    for i = 1:length(conds)
        if conds(i) == 2
            easyOrHard = 2;
            conds(i) = 1;
        elseif conds(i) == 3
            easyOrHard = 1;
            conds(i) = 2;
        elseif conds(i) == 4
            conds(i) = 2;
        end
    end
    resp = dataValid(:,2);
    rt = dataValid(:,3);
    acc = dataValid(:,4);
    
    % do descriptive stats    
    minRT(idxsubj) = min(rt(rt>=0.1)); % not used
    

    %%
    %Pool conditions according to design matrix & calculate stats
    unique_conds = unique(conds);
    for level = 1:length(unique_conds) % for all conditions
        
        trialn = conds == level;
        dataRaw{idxsubj,level}  = [resp(trialn) rt(trialn)];
        
        data2fit{idxsubj,level} = data_stats(dataRaw{idxsubj,level});

        % add some info to data2fit (the condition labels as a string, and
        % the minRT - in the case that the min rt is smaller than the
        % non-decision time - let's leave this for now, just in case)
        data2fit{idxsubj,level}.cond = validLabs{trialn};
        data2fit{idxsubj,level}.minRT = round(minRT(idxsubj),2);
        %parange(end,1) = round(minRT(idxsubj),2);%constrain the lower bound of T0 to the shortest RT
    end
    %% fit the model
    if easyOrHard == 2
        [bestpar{idxsubj,1},bestval{idxsubj,1},BIC{idxsubj,1}]=fitparams_refine_template_RDK('fiterror_cell_HARD',Model_Feature,{data2fit{idxsubj,1:2}},randiter,nosession,[],parange,bayesian);%#ok
    elseif easyOrHard == 1
        [bestpar{idxsubj,1},bestval{idxsubj,1},BIC{idxsubj,1}]=fitparams_refine_template_RDK('fiterror_cell_EASY',Model_Feature,{data2fit{idxsubj,1:2}},randiter,nosession,[],parange,bayesian);%#ok
    end
end

save(savename,'bestpar','bestval','BIC','rseed','settings');%,'bestpar_PA','bestval_PA','BIC_PA')

end
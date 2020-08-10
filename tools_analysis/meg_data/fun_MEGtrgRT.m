function fun_MEGtrgRT(filenames,behavfile,savebehav,settings,trfiles)%#ok


MEG_RT=[];alltrl = []; allconditionlabels = {};tmplabs ={};%#ok
load(behavfile) % matlab behavioural file

%runid = repmat((1:numel(filenames))',[1,2]).*2-repmat([1 0],[numel(filenames),1]); % tells the block ids within one run

% runid looks like:
%     runid =
%
%      1     2   - is the run, and which blocks
%      3     4   - is run 2 with those blocks
%      5     6   - etc
%      7     8
%      9    10

subjectswitcher = input('enter pilot subject (e.g. S01): ','s');
%subjectswitcher = string(subjectswitcher);

switch subjectswitcher
    case 'S01'
        % 'S01' - 6 runs, with 6 blocks per run
        runid = [1:6;7:12;13:18;19:24;25:30;31:36];
    case 'S02'
        % S02 - 3 runs, run 1 no data, run 2 is two blocks, run 3 is 4 blocks
        % runid = [[];1:2;3:6];
        runid = NaN(3,4); runid(2,1:2) = 1:2; runid(3,:) = 3:6;
    case 'S03'
        % S03 - 5 runs, 6 blocks per run ex last run is 1 block, and first
        % 8 trials of run 4 belong to the end of the last block in run 3 (we miss 13 trials = 21 - the 8 in run 4
        % runid = [1:6;7:12;13:18;19:24;25]; % so block 19 has 8 trials from end of block 18, and block 18 is missing 21 trials (8 of which are in the next block)
        runid = [1:6;7:12;13:18;19:24;25,NaN,NaN,NaN,NaN,NaN];
end

conditions = {'LcLr', 'LcHr','HcLr','HcHr'}; % 2x2 coherence and rule

for runi = 1:numel(filenames)
    %%
    clear trl conditionlabels
    D = spm_eeg_load(filenames{runi});
    
    
    % double check that trigger responses correspond to behavioral results
    trig = D(D.indchannel('STI101'),:)/1e6;  % 1e6 new scaling introduced by SPM12;
    ttime = D.time;
    
    
    
    
    % ale's code - block = ResponseArray(ismember(ResponseArray(:,1),runid(runi,:)),:); % (d.stim_mat_all - this selects the elements from the behavioural file from the block/s associated with the run with runid + the reaction time data - so runi on the 3rd element
    % block = d.stim_mat_all(:,:,runid(runi,:));
    % so can do the above and work with 3 dimensions, but requires us to change block refs to suit, eg :
    % if numel(trigger_onset)~= size(block,1)*size(block,3)
    % if sum(trialID ~= block(:,10,:))
    % or can reshape d.stim_mat_all as below, but concerned I'll lose block related information? maybe can add a row...
    % block = permute(d.stim_mat_all,[1,3,2]); % swaps 3rd and 2nd dimensions
    % block = reshape(block,[],size(d.stim_mat_all,2),1); % reshapes so rows are autocalculated, columns defined by stimulus matrix columns, and forces 3rd dimension into the 2d structure
    % block = block(ismember(block(:,1),runid(runi,:))); - but this doesn't
    % work because you don't have a column saying what's what. so just take 3rd
    % dimensions needed?
    
    vecrunid = runid(runi,:); % select vector of blocks in the run
    vecrunid(isnan(vecrunid)) = []; % get rid of any NaNs if they exist
    block = d.stim_mat_all(:,:,vecrunid); % pull the stimulus matrix info for those blocks
    block = permute(block,[1,3,2]); % swaps 3rd and 2nd dimensions
    block = reshape(block,[],size(d.stim_mat_all,2),1); % reshapes so rows are autocalculated, columns defined by stimulus matrix columns, and forces 3rd dimension into the 2d structure
    
    if strcmp(subjectswitcher,'S03')
        if runi == 3
            % need to remove the last 21 trials and save 8 for the next run
            S03savetheseblocks = block(end-7:end,:);
            block(end-20:end,:) = [];
        elseif runi == 4
            % need to get the last 8 trials of run 3 at the beginning
            block = [S03savetheseblocks;block];
        end
    end
    
    switch subjectswitcher
        case 'S01'
            % for S01, who has no distinct coh v trial onset triggers
            triggers = d.stim_mat_all(:,9)'; % trial type triggers
            keys     = [4096 8192 16384 -32768]; % index medium ring little - on right button box
            resp_onset = (find(ismember(round(diff(trig)),keys))+1)';
            trigger_onset  = (find(ismember(round(diff(trig)),triggers))+1)'; % onset of all triggers
            trialID  = [0; round(diff(trig)')]; trialID = trialID(trigger_onset); % for S01 - replaces other two
            
            if numel(trialID) ~= numel(block(:,9)) || sum(trialID ~= block(:,9))
                    tmp_tronset = [];
                    tmp_trialID =[];
                    iresp = 0;
                    twindow = 2000;
%                     nilresp = 0;
                    for itr = 1:numel(block(:,9))
                        iresp = iresp+1; % counter for button presses
                        
%                         if nilresp
%                             
%                             searchtwin=[0 round(diff(trig(resp_onset(iresp)-1700:resp_onset(iresp))))]; % cut out the time window in the trial interval preceding the button press
%                             ttrigger = find(searchtwin>0 & searchtwin <= 128); % find any triggers in there - might be able to make this smaller? This is in case is codes two 64 triggers as one, but I don't think that can happen
%                             ttrigger = ttrigger(end); % save only the most recent trigger (in case we took any from a previous trial)
%                             
%                         else
                          searchtwin=[0 round(diff(trig(resp_onset(iresp)-twindow:resp_onset(iresp))))]; % cut out the time window in the trial interval preceding the button press
        
                            ttrigger = find(searchtwin>0 & searchtwin <= 128); % find any triggers in there - might be able to make this smaller? This is in case is codes two 64 triggers as one, but I don't think that can happen
                            % is this an appropriate button press response?
                            if isempty(ttrigger) % this is a cue accept response, not a trial response
                                iresp = iresp+1;% go to the next button press
                                if runi == 3 && iresp == 2 && itr == 1; iresp = iresp+1;end % resp 1 accepts end of training for keyswap on this run, resp 2 is the accept cue resp, resp 3 is the start of the trials
                                if runi == 5 && iresp == 2 && itr == 1; iresp = iresp+1;end % not sure why there are two responses here, but there are
                                if runi == 6 && iresp == 2 && itr == 1; iresp = iresp+1;end % same - unsure - wondering if it's something in the script, but it didn't happen on run 4...
                                % do this again
                                searchtwin=[0 round(diff(trig(resp_onset(iresp)-twindow:resp_onset(iresp))))]; % cut out the time window in the trial interval preceding the button press
                                ttrigger = find(searchtwin>0 & searchtwin <= 128); % find any triggers in there - might be able to make this smaller? This is in case is codes two 64 triggers as one, but I don't think that can happen
                            end
                            ttrigger = ttrigger(end); % save only the most recent trigger (in case we took any from a previous trial)
                            
%                         end

                        % check that there is a response between upcoming triggers
                        % however, currently the code only goes one trigger
                        % ahead, which will find the response you just used
                        % to code ttrigger - need to go two ahead, but then
                        % if you have two fake triggers (due to the slow
                        % analogue trigger system) it won't work. needs
                        % thought.
%                         searchtwin = [0 round(diff(trig(ttrigger+100:ttrigger+1700)))]; % search from 100ms to 1800ms (i.e. 100ms plus trial of 1700ms) - must be less than 2.7s
%                         trigs_exist = find(searchtwin>0 & searchtwin <= 128);
%                         if trigs_exist
%                             searchtwin = [0 round(diff(trig(ttrigger:trigs_exist)))];
%                             resp_exists = find(searchtwin>2000);
%                             if ~resp_exists
%                                 nilresp = 1;
%                             end
%                         end

                        
                        tmp_trialID(itr) = block(itr,9);
                        
                        tmp_tronset(itr) = resp_onset(iresp)-twindow + ttrigger-1;
                    end                            % alternatively, can code by responses, idx for
                            % trial, then a counter as an idx for responses
                            % - where no tri
                    trigger_onset   = tmp_tronset';
                    trialIDtrg  =  tmp_trialID';
                end                            % alternatively, can code by responses, idx for
                            % trial, then a counter as an idx for responses
                            % - where no tri
        
        otherwise
            % find onset trials, coherence, trial type, and rt
            triggers = d.stim_mat_all(:,10,1)'; % trial type triggers
            coh_trigger = 1; % value of coherence onset trigger
            keys     = [4096 8192 16384 -32768]; % index medium ring little - on right button box
            trigger_onset  = (find(round(diff(trig)) == 1)+1)'; % onset of all triggers - SUPERFLOUS
            resp_onset_neg = (find(ismember(round(diff(trig)),-keys))+1)';
            coh_onset = (find(ismember(round(diff(trig)),coh_trigger))+1)'; % coherence trigger onset
            trial_type_onset = (find(ismember(round(diff(trig)),triggers))+1)'; % trial type trigger onset
            trial_type_onset_neg = (find(ismember(round(diff(trig)),-triggers))+1)'; % trial type trigger onset
            
            cohtrigsID  = [0; round(diff(trig)')]; cohtrigsID = cohtrigsID(coh_onset); % gets the amplitude of the trigger (which should be the thing you specified as the trigger) and then strips anything that is not a coherence onset trigger
            trialIDtrg  = [0; round(diff(trig)')]; trialID = trialIDtrg(trial_type_onset); % same for trial type
            
            if numel(trialID) ~= numel(block(:,10)) || sum(trialID ~= block(:,10))
                
                if numel(coh_onset) == numel(block(:,10)) % if coh  onset is reliable
                    tmp_tronset = [];
                    tmp_trialID =[];
                    for ic = 1:numel(coh_onset)
                        
                        searchtwin=[0 round(diff(trig(coh_onset(ic)-400:coh_onset(ic))))];
                        ttrigger = find(searchtwin>1 & searchtwin <= max(triggers));
                        ttrigger = ttrigger(1);
                        
                        tmp_trialID(ic) = block(ic,10);
                        
                        
                        tmp_tronset(ic) = coh_onset(ic)-400 + ttrigger-1;
                        
                    end
                    trial_type_onset   = tmp_tronset';
                    trialIDtrg  =  tmp_trialID';
                else
                    if strcmp(subjectswitcher,'S03')
                        % so to recover the triggers using the rt from
                        % matlab, we can go to response 2, then search back
                        % for the rt time window (probably plus 100ms or
                        % so) and if there's nothing, check the next
                        % response - but keep track of these, maybe even
                        % add another trigger, so you can see it in the
                        % plot.
                        if runi == 4 || runi == 5
                            iresp = 0;
                            for itr = 1:numel(block(:,10))
                                iresp = iresp+1; % counter for button presses
                                if iresp == 1
                                    ttrigger = resp_onset_neg(itr)+300;
                                elseif iresp == 9
                                    iresp = 0;
                                    ttrigger = ttrigger+1800;
                                else
                                    ttrigger = ttrigger+1800;
                                end
                                tmp_tronset(itr) = ttrigger;
                                
                                tmp_trialID(itr) = block(itr,10);

                            end
                            trigger_onset   = tmp_tronset';
                            trialIDtrg  =  tmp_trialID';
                            
                            if runi == 4
                                trial_type_onset = [trial_type_onset;trigger_onset(length(trial_type_onset)+1:end)];
                            elseif runi == 5
                                trial_type_onset = trigger_onset;
                            end
                            
                            %157797 is cue acceptor trigger end
                            
                            %bad_trigger = 255;
                            %bad_triggers = (find(trig == bad_trigger)+1)';% (find(ismember(round(diff(trig)),bad_trigger))+1)';
                            %first_bad_trigger = bad_triggers(1);
                            
                            % so now from here, we can go back to the big
                            % difference between the last response of the
                            % last subblock of 8, and the accept cue response
                            % of the subblock of 8 that's broken - then go 8
                            % by 8 thereafter (1.8s between onsets commencing after response trigger off).
                            
                            %searchtwin=[0 round(diff(trig(first_bad_trigger-14400:first_bad_trigger)))];
                            %ttrigger = find(abs(searchtwin>500)-1800);
                            
                            
                        end
                    end
                end
                
                
                
            end
            
    end
    
    switch subjectswitcher
        case 'S01'
            figure;plot(trig);hold on;plot(trigger_onset,trig(trigger_onset),'-o');
        otherwise
            if runi < 5
                figure;plot(trig);hold on;plot(trial_type_onset,trig(trial_type_onset),'-o');
            end
        end

    %check
    if numel(trigger_onset) ~= size(block,1)
        if runi < 5
        error('Inconsistent number of onset triggers vs trials')
        end
    end
    switch subjectswitcher
        case 'S01'
            if numel(trialID) ~= numel(block(:,9)) || sum(trialID ~= block(:,9)) % for S01
                if sum(trialIDtrg ~= block(:,9))
                    error('Discrepancy between triggers coded conditions and trials')
                end
            end
        otherwise
            if numel(trialID) ~= numel(block(:,10)) || sum(trialID ~= block(:,10))
              if runi < 5
                if sum(trialIDtrg ~= block(:,10)) % compare with the d.stim_mat_all column with trigger amplitudes
                    error('Discrepancy between triggers coded conditions and trials')
                end
              end
            end
    end
    
    % find responses and empty trials
    clear tmp_RT tKeyp responset
    
    switch subjectswitcher
        case 'S01'
            % for S01
            for ions = 1:numel(trigger_onset)
                
                if ions < numel(trigger_onset)
                    tmp    =trig(trigger_onset(ions):(trigger_onset(ions+1)-1));
                else
                    tmp    =trig(trigger_onset(ions):end);
                end
                
                tmpres = (find(ismember(diff(tmp),keys),1,'first')+1)';%first to avoid multiple button presses % time of button press within trial period (iti to iti)
                
                if ~isempty(tmpres)
                    responset(ions)  = trigger_onset(ions)+tmpres-1;%#ok % time of button press within block
                    tmp_RT(ions)       = ttime(responset(ions))-(ttime(trigger_onset(ions))+.034);%#ok 34ms projector's delay % this is reaction time (i.e. time of button press from coherence onset)
                else
                    responset(ions) = 0;%#ok
                    tmp_RT(ions)      = 0;%#ok
                end
            end
        otherwise
            for ions = 1:numel(trial_type_onset)
                if runi < 5
                if ions < numel(trial_type_onset)
                    tmp    =trig(trial_type_onset(ions):(trial_type_onset(ions+1)-1));
                else
                    tmp    =trig(trial_type_onset(ions):end);
                end
                end
                
                if runi < 5
                tmpres = (find(ismember(diff(tmp),keys),1,'first')+1)';%first to avoid multiple button presses % time of button press within trial period (iti to iti)
                else
                    tmpres = [];
                end
                if ~isempty(tmpres)
                    responset(ions)  = trial_type_onset(ions)+tmpres-1;%#ok % time of button press within block
                    if length(coh_onset) < ions
                        responset(ions) = 0;
                        tmp_RT(ions) = 0;
                    else
                        tmp_RT(ions)       = ttime(responset(ions))-(ttime(coh_onset(ions))+.034);%#ok 34ms projector's delay % this is reaction time (i.e. time of button press from coherence onset)
                    end
                else
                    responset(ions) = 0;%#ok
                    tmp_RT(ions)      = 0;%#ok
                end
                
            end
    end
    
    %check
    % rts = block(:,13); % will need to get rts from d.rts (this is a column in this script)
    
    rts = d.rt(vecrunid,:)';rts = rts(:); % transpose relevant rows of d.rt and then place each block under one another in the same row
    
    if strcmp(subjectswitcher,'S03')
        if runi == 3
            % need to remove the last 21 trials and save 8 for the next run
            S03savetheserts = rts(end-7:end,:);
            rts(end-20:end,:) = [];
        elseif runi == 4
            % need to get the last 8 trials of run 3 at the beginning
            rts = [S03savetheserts;rts];
        end
    end
    
    idx = rts~=0 & rts>0; % where there is a rt
    
    if strcmp(subjectswitcher,'S03')
        if runi == 3
            tmp_RT(end) = rts(end);
        end
    end
          

    tmp = [rts(idx) transpose(tmp_RT(idx))]; % this compares rts here with rts there
    tmp = corr(tmp);
    
    
    if tmp(2,1) <0.99 || isnan(tmp(2,1))   % takes the correlation from the output of corr
        if runi < 4            
            error('large discrepancy between stim RT and MEG RT!')
        end
    end
    
  
    
    % re-compose matrix of behavioural data
    % so block looks like:
    %  1)  cue direction (1-4)
    %  2)  cue direction in degrees
    %  3)  dot motion direction condition (1-8)
    %  4)  dot motion direction in degrees
    %  5)  coherence difficulty (1 = easy, 2 = hard)
    %  6)  matching distance from cue direction (absolute value) - used to calc match and match difficulty
    %  7)  match blue (1) or match orange (2)
    %  8)  matching difficulty (1 = easy, 2 = difficult)
    %  9)  gives you a unique number for each trial condition
    % 10)  gives a number based on 9 for meg triggers
    for i = 1:length(block) % create a row that gives you a number for each condition in your 2x2
        if block(i,5)==1 && block(i,8)==1
            conditionlabels(i,1) = string(conditions{1});
        elseif block(i,5)==1 && block(i,8)==2
            conditionlabels(i,1) = string(conditions{2});
        elseif block(i,5)==2 && block(i,8)==1
            conditionlabels(i,1) = string(conditions{3});
        elseif block(i,5)==2 && block(i,8)==2
            conditionlabels(i,1) = string(conditions{4});
        end
    end
    accuracyrow = d.correct(vecrunid,:)'; accuracyrow = accuracyrow(:); % may as well put this in the new matrix - so transpose relevant rows and then place them one after another in a column
    if strcmp(subjectswitcher,'S03')
        if runi == 3
            % need to remove the last 21 trials and save 8 for the next run
            S03savetheseaccs = accuracyrow(end-7:end,:);
            accuracyrow(end-20:end,:) = [];
        elseif runi == 4
            % need to get the last 8 trials of run 3 at the beginning
            accuracyrow = [S03savetheseaccs;accuracyrow];
        end
    end
    
    %MEG_RT = cat(1,MEG_RT,[block(:,[1:4 12 14 15 ]), tmp_RT']);% 1block# 2coh 3act 4cond 5trigger 6resKey 7acc 8rt % compiles matlab trial information with the new MEG generated rts
    if runi < 5
    MEG_RT = cat(1,MEG_RT,[block(:,[1 3 5 8 10]), conditionlabels, accuracyrow, rts, tmp_RT']);
    end
    % re-composed matrix looks like:
    %  1)  cue direction (1-4)
    %  2)  dot motion direction condition (1-8)
    %  3)  coherence difficulty (1 = easy, 2 = hard)
    %  4)  matching difficulty (1 = easy, 2 = difficult)
    %  5)  gives a number based on 9 for meg triggers
    %  6)  conditions (HcHr, HcLr, LcHr, LcLr)
    %  7)  accuracy (0, 1, or -1 for missed trials)
    %  8)  rts from behavioural data
    %  9)  rts from MEG
    
    %trial definition for epoching
    
    % set the trial window, including buffer for edge effects during filtering
    % - in the textbook - lowest frequency will determine
    switch subjectswitcher
        case 'S01'
            twind = [-1500 2500];%in samples at 1000Hz sampling % what time window 1500 = 1.5s around the coherence onset
            trl      = (trigger_onset+34)+twind;%[trial-beginning trial-end] accounts for 34ms projector's delay.
            trl(:,3) = -1500;% epoched trials begin 1500ms before coherence onset
        otherwise
            twind = [-1500 2500];%in samples at 1000Hz sampling % what time window 1500 = 1.5s around the coherence onset
            trl      = (coh_onset+34)+twind;%[trial-beginning trial-end] accounts for 34ms projector's delay.
            trl(:,3) = -1500;% epoched trials begin 1500ms before coherence onset
    end
    
    %add to session trial def
    alltrl = cat(1,alltrl,trl);
    
    tmplabs = cat(1,tmplabs,conditionlabels);
    conditionlabels = cellstr(conditionlabels); % need to convert for the preprocessing
    save(trfiles{runi},'trl','conditionlabels')
end

% don't think we need this
% idbad = MEG_RT(:,7)<=0 | MEG_RT(:,7)>2.1; % defines bad trials based on reaction time will need to check long RT and column number
% MEG_RT(idbad,8) = 0;
% MEG_RT(idbad,[6 7]) = -1;

%%

% add trial numbers to MEG_RT
MEG_RT(:,end+1) = 1:length(MEG_RT(:,1));
%allconditionlabels = {conditions{MEG_RT(:,8)}}';%#ok
allconditionlabels = MEG_RT(:,6);

%check for consistency
if runi < 5
if sum(strcmp(allconditionlabels,tmplabs)) ~= numel(tmplabs); error(['incongruent condition labels:' filenames{runi}]);end
end

save(savebehav,'MEG_RT','settings','alltrl','allconditionlabels')
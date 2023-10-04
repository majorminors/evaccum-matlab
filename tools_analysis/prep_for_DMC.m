%% setup

clear all


rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';
scriptdir = fullfile(rootdir,'tools_analysis'); cd(scriptdir)
datadir = fullfile(rootdir,'data','meg_pilot_4','behavioural'); addpath(datadir);
table_save_file = fullfile(datadir,'DMC_table_%s.csv');

allSubjects = importParticipants();

idx=0;
for subjectidx = 1:numel(allSubjects)
    clearvars d block trial
    disp('starting with subject: ');
    
    thisSubject = allSubjects{subjectidx}; % let's make this easier to reference
    disp(thisSubject.id)
    % let's skip participants who are not useable
    if ~thisSubject.usable
        disp('marked not usable - skipping')
        continue
    end
    
    behavData = fullfile(datadir,[thisSubject.id '_EvAccum.mat']);
    load(behavData,'d');
    
    
    %% create output table
    
    block = size(d.rt,1);
    
    while block
        
        trial = size(d.rt,2);
        
        while trial
            idx=idx+1;
            
            Var1(idx)=idx;
            s(idx)=thisSubject.num;
            % we don't want this kind of granularity yet
            %         if thisSubjectData.stim_array{trial}.cue_dir == 1 || thisSubjectData.stim_array{trial}.cue_dir == 3
            %             S{idx}='s1';
            %         elseif thisSubjectData.stim_array{trial}.cue_dir == 2 || thisSubjectData.stim_array{trial}.cue_dir == 4
            %             S{idx}='s2';
            %         end
            if d.stim_mat_all(trial,5,block) == 1
                Coh{idx}='cohE';
            elseif d.stim_mat_all(trial,5,block) == 2
                Coh{idx}='cohH';
            end
            if d.stim_mat_all(trial,8,block) == 1
                Angle{idx}='angleE';
            elseif d.stim_mat_all(trial,8,block) == 2
                Angle{idx}='angleH';
            end
            if strcmp(d.correct_resp{block,trial},'RR')
                R{idx}='r1';
                % let's just code for correctness for our first modelling attempt
                if d.correct(block,trial) == 1
                    S{idx} = 's1'; else; S{idx} = 's2';
                end
            elseif strcmp(d.correct_resp{block,trial},'RB')
                R{idx}='r2';
                % let's just code for correctness for our first modelling attempt
                if d.correct(block,trial) == 1
                    S{idx} = 's2'; else; S{idx} = 's1';
                end
            end
            RT(idx)=d.rt(block,trial);
            
            trial = trial-1;
        end
        block = block-1;
    end
    
end

Var1=Var1';s=s';S=S';Coh=Coh';Angle=Angle';R=R';RT=RT';

outputTable = table(Var1,s,S,Coh,Angle,R,RT);
thisSaveFile = sprintf(table_save_file,'all');
fprintf('saving table as %s\n', thisSaveFile);
writetable(outputTable,thisSaveFile); % save to csv

outputTable = table(Var1,s,S,Coh,R,RT);
thisSaveFile = sprintf(table_save_file,'coh');
fprintf('saving table as %s\n', thisSaveFile);
writetable(outputTable,thisSaveFile); % save to csv

outputTable = table(Var1,s,S,Angle,R,RT);
thisSaveFile = sprintf(table_save_file,'ang');
fprintf('saving table as %s\n', thisSaveFile);
writetable(outputTable,thisSaveFile); % save to csv


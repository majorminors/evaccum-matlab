clear all

theseSubjects = [5,15,17,30];

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';
scriptdir = fullfile(rootdir,'tools_analysis');
datadir = fullfile(rootdir,'data','meg_pilot_4'); addpath(datadir);
behavDataDir = fullfile(datadir,'behavioural');

allSubjects = importParticipants();

for subjectidx = 1:numel(allSubjects)
    
    thisSubject = allSubjects{subjectidx}; % let's make this easier to reference
    
    if ~ismember(thisSubject.num,theseSubjects)
        continue
    end
    
    disp('doing subject: ')
    disp(thisSubject.id)
    
    d1 = load(fullfile(behavDataDir,[thisSubject.id '_1_EvAccum.mat']),'d');
    d2 = load(fullfile(behavDataDir,[thisSubject.id '_2_EvAccum.mat']),'d');
    
    outdname = fullfile(behavDataDir,[thisSubject.id '_EvAccum.mat']);
    
    totalRuns = thisSubject.runblks(1);
    r1 = thisSubject.runblks(2);
    r2 = thisSubject.runblks(3);
    
    d.rt = [d1.d.rt(1:r1,:);d2.d.rt(1:r2,:)]; if size(d.rt,1) ~= totalRuns; error('rts dont match blocks'); end
    d.correct = [d1.d.correct(1:r1,:);d2.d.correct(1:r2,:)]; if size(d.correct,1) ~= totalRuns; error('correct dont match blocks'); end

    sm1 = d1.d.stim_mat_all(:,:,1:r1);
    sm2 = d2.d.stim_mat_all(:,:,1:r2);
    if d1.d.easy_rule ~= d2.d.easy_rule
        
        warning('rule values dont match---matching second file to first')
        
        % so now we need to swap the easy rule and hard rule of the second
        % stim mat to match the first, and we'll use that as the base
        
        d.easy_rule = d1.d.easy_rule;
        d.hard_rule = d1.d.hard_rule;
        
        for i1 = 1:size(sm2,3)
            
            for i2 = 1:size(sm2,1)
                
                %  8)  matching difficulty (1 = easy, 2 = difficult)
                
                if sm2(i2,8,i1) == 1
                    sm2(i2,8,i1) = 2;
                elseif sm2(i2,8,i1) == 2
                    sm2(i2,8,i1) = 1;
                end
                
            end
            
        end; clear i1 i2
        
    end
    d.stim_mat_all = cat(3,sm1,sm2);
    if size(d.stim_mat_all,3) ~= totalRuns; error('stim mat not correct size'); end
    
    if round(d1.d.easy_coherence,2) ~= round(d2.d.easy_coherence,2); error('coherence doesnt match'); end
    
    % save d d1 d2
    save(outdname,'d','d1','d2');
    
end

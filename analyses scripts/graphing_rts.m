%%
conds = a.d.stim_mat_all(:,[5 8]);
rts   = a.d.rt(1,:)';
acc   = a.d.correct(1,:)';
%%
count = 0;
for i = 1:2
    for ii = 1:2
        count = count + 1;
        idx = conds(:,1)==i & conds(:,2)==ii;
        idx = idx & acc;
        subplot(2,2,count)
        hist(rts(idx))
        
        
    end
end
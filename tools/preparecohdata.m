function data = preparedata(d)
%%
tmp=[d.stim_mat_all(:,6) d.correct'];
tmp = sortrows(tmp,1);
ints = unique(tmp(:,1));
data =[];
for i = 1:numel(ints)
    
    idx = tmp(:,1)==ints(i);
    data(i,1)=ints(i);
    data(i,2)=sum(tmp(idx,2));
    data(i,3)=sum(idx);
    
end
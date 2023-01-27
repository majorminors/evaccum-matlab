%parameter recovery for LBA script

rng(17,'twister') % for reproducibility

%fix parameters
model_var = [1 4];
simpars = [2.14228548713214, 1.10671026363223, 0.149709201018512, 0.36448066743315, 3.9544457786746, 9.98344328857036,0.912320163409241,0.891749025721009];
[num_param,parLL,parHL,parLH,parHH]=getModelParam_cell_RDK(model_var,2,simpar);
samplesize = 1000;%per condition

%generate simulated RT distributions

parLabs = {'parLL','parHL','parLH','parHH'};
model = 'LBA_spec_FT3_c_even_template';   

clear rt_iter act_iter data
tmpdata=[];
for ivar = 1:4
    
eval(sprintf('model_param = %s;',parLabs{ivar}));
numAct = 2;
[rt_iter,act_iter]= feval (model, model_param.N,model_param.B,model_param.C0,model_param.Ame,model_param.Astd,model_param.T0,ones(1,numAct)./numAct);

%only positive RTs
if sum(rt_iter<=0 | act_iter==0)~=0
        indx=find(rt_iter<=0 | act_iter==0 );
        act_iter(indx)=[];
        rt_iter(indx)=[];
end

ntrials = max(size(rt_iter'));
idxtr   = randsample(ntrials,samplesize,'false');
ntrials = size(idxtr);
tmpdata = [tmpdata; [nan(ntrials),ivar.*ones(ntrials),act_iter(idxtr)',rt_iter(idxtr)',ones(ntrials)]];

end

%scramble trials in random order
ntrials = max(size(tmpdata));
tmpdata = tmpdata(randsample(ntrials,ntrials,'false'),:);
data(1).data = num2cell(tmpdata);
data(1).id = 1;
%fit simulated data
settings.randiter = 100;
settings.nosession= 100;
settings.save = 0;
settings.modfeat = model_var;
settings.rseed   = 17;
[bestpar,bestval,BIC]=fun_fitbehav_LBA_PD(settings,data);
%compare estimated parameters with simulated ones
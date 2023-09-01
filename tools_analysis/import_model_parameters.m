clear all

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';

subs = readtable(fullfile(rootdir,'dmc_lba','subjectParams.csv'));
mns = readtable(fullfile(rootdir,'dmc_lba','meanParams.csv'));

plot([mns.b_ecer mns.b_echr mns.b_hcer mns.b_hchr],'linestyle','none','marker','o')
plot([mns.v1_ecer mns.v1_echr mns.v1_hcer mns.v1_hchr],'linestyle','none','marker','o')

plot([subs.b_ecer subs.b_echr subs.b_hcer subs.b_hchr]','linestyle','-','marker','o')
plot([subs.v1_ecer subs.v1_echr subs.v1_hcer subs.v1_hchr]','linestyle','-','marker','o')
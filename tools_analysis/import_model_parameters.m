clear all

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';

dmcdir = fullfile(rootdir,'dmc_lba');

subs = readtable(fullfile(dmcdir,'subjectParams.csv'));
mns = readtable(fullfile(dmcdir,'meanParams.csv'));

% plot([mns.b_ecer mns.b_echr mns.b_hcer mns.b_hchr],'linestyle','none','marker','o')
% plot([mns.v1_ecer mns.v1_echr mns.v1_hcer mns.v1_hchr],'linestyle','none','marker','o')
% 
% plot([subs.b_ecer subs.b_echr subs.b_hcer subs.b_hchr]','linestyle','-','marker','o')
% plot([subs.v1_ecer subs.v1_echr subs.v1_hcer subs.v1_hchr]','linestyle','-','marker','o')

% Load data from the subs table
subs.b_ec = mean(subs{:, {'b_ecer', 'b_echr'}}, 2);
subs.b_hc = mean(subs{:, {'b_hcer', 'b_hchr'}}, 2);
subs.b_er = mean(subs{:, {'b_ecer', 'b_hcer'}}, 2);
subs.b_hr = mean(subs{:, {'b_echr', 'b_hchr'}}, 2);
b_columns = subs(:, {'Var1' 'b_ecer', 'b_echr', 'b_hcer', 'b_hchr', 'b_ec', 'b_hc', 'b_er', 'b_hr'});
b_vals = stack(b_columns,b_columns.Properties.VariableNames(2:end));
b_vals.Properties.VariableNames{1} = 'Subject';
b_vals.Properties.VariableNames{2} = 'Param';
b_vals.Properties.VariableNames{3} = 'Value';

writetable(b_columns,fullfile(dmcdir,'b_for_jasp.csv'))

subs.v1_ec = mean(subs{:, {'v1_ecer', 'v1_echr'}}, 2);
subs.v1_hc = mean(subs{:, {'v1_hcer', 'v1_hchr'}}, 2);
subs.v1_er = mean(subs{:, {'v1_ecer', 'v1_hcer'}}, 2);
subs.v1_hr = mean(subs{:, {'v1_echr', 'v1_hchr'}}, 2);
v1_columns = subs(:, {'Var1' 'v1_ecer', 'v1_echr', 'v1_hcer', 'v1_hchr', 'v1_ec', 'v1_hc', 'v1_er', 'v1_hr'});
v1_vals = stack(v1_columns,v1_columns.Properties.VariableNames(2:end));
v1_vals.Properties.VariableNames{1} = 'Subject';
v1_vals.Properties.VariableNames{2} = 'Param';
v1_vals.Properties.VariableNames{3} = 'Value';

writetable(v1_columns,fullfile(dmcdir,'v1_for_jasp.csv'))

% so v1: f(3)=40.043 MS=8.262 p<0.001
% all post-hoc t-tests significantly difference <0.001, except HcEr and
% HcHr (0.412)!

% in contrast, the b looks actually like they aren't different!


writetable(join(b_columns,v1_columns,'Keys','Var1'),fullfile(dmcdir,'params_for_jasp.csv'))

forPlot = vertcat(b_vals,v1_vals);
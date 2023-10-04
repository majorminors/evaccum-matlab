clear all

rootdir = '/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab';

subs = readtable(fullfile(rootdir,'dmc_lba','subjectParams.csv'));
mns = readtable(fullfile(rootdir,'dmc_lba','meanParams.csv'));

% plot([mns.b_ecer mns.b_echr mns.b_hcer mns.b_hchr],'linestyle','none','marker','o')
% plot([mns.v1_ecer mns.v1_echr mns.v1_hcer mns.v1_hchr],'linestyle','none','marker','o')
% 
% plot([subs.b_ecer subs.b_echr subs.b_hcer subs.b_hchr]','linestyle','-','marker','o')
% plot([subs.v1_ecer subs.v1_echr subs.v1_hcer subs.v1_hchr]','linestyle','-','marker','o')

% Load data from the subs table
b_columns = subs(:, {'Var1' 'b_ecer', 'b_echr', 'b_hcer', 'b_hchr'});

Subjects = b_columns.Var1;
Values = table2array(b_columns(:, 2:end)); % If values start from column 2
Conditions = b_columns.Properties.VariableNames(2:end); % Assuming b_columns is a table

% If you need to repeat Subjects and Conditions according to the number of conditions
Subjects_repeated = repelem(Subjects, numel(Conditions));
Conditions_repeated = repmat(Conditions, numel(Subjects), 1);
Values_flattened = reshape(Values.',[],1);

v1_columns = subs(:, {'Var1' 'v1_ecer', 'v1_echr', 'v1_hcer', 'v1_hchr'});

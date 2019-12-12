%% set up

close all;
clearvars;
clc;

% enter filename
datafilenames = ["S01_EvAccum_stitched", "S02_EvAccum", "S03_EvAccum", "S04_EvAccum", "S05_EvAccum",...
    "S06_EvAccum", "S07_EvAccum", "S08_EvAccum", "S09_EvAccum", "S10_EvAccum", "S11_EvAccum",...
    "S12_EvAccum", "S13_EvAccum", "S14_EvAccum", "S15_EvAccum"];
datasuffix = '_datasummary';



% directory mapping
rootdir = 'C:\Users\doria\Google Drive\04 Research\05 Evidence Accumulation\01 EvAccum Code';
datadir = fullfile(rootdir, 'data');
analysisdir = fullfile(datadir, 'analysis');
%subdir = fullfile(analysisdir, datafilenames);


for i = 1:length(datafilenames)
    
% load file
data = load(fullfile(analysisdir, datafilenames(i), strcat(datafilenames(i),  datasuffix)));

data_summaries.(datafilenames(i)) = data.data_summary;
clear data;

mean_summary.(datafilenames(i)) = mean(data_summaries.(datafilenames(i)),3);

if i == 1
    
    [rows, columns] = size(mean_summary.(datafilenames(i)));
    data_mat = zeros(rows,columns,length(datafilenames));
    clear rows columns
    
end

data_mat(:,:,i) = mean_summary.(datafilenames(i));   

end

mean_summary.total = mean(data_mat,3);

% you'll end up with a table like this (without the labels):
%                 | easy coh |  hard coh | easy rule | hard rule |
% average rt      |          |           |           |           |
% percent correct |          |           |           |           |
% average correct |          |           |           |           |
%                      easy coherence        hard coherence
%                 |easy rule | hard rule | easy rule | hard rule |
% average rt      |          |           |           |           |
% percent correct |          |           |           |           |
% average correct |          |           |           |           |

% repeated 3x horizontally: | all results | results for blue | results for orange |
% for each participant included in data file names, as well as a mean
% summary of all

totalsem.rt = [nansem(data_mat(3,1,:)),nansem(data_mat(3,2,:)),nansem(data_mat(3,3,:)),nansem(data_mat(3,4,:));...
    nansem(data_mat(6,1,:)),nansem(data_mat(6,2,:)),nansem(data_mat(6,3,:)),nansem(data_mat(6,4,:))];

totalsem.pc = [nansem(data_mat(2,1,:)),nansem(data_mat(2,2,:)),nansem(data_mat(2,3,:)),nansem(data_mat(2,4,:));...
    nansem(data_mat(5,1,:)),nansem(data_mat(5,2,:)),nansem(data_mat(5,3,:)),nansem(data_mat(5,4,:))];

% you'll end up with tables like this (without the labels) for rts and pc:
%          | easy coh |  hard coh | easy rule | hard rule |
% sem      |          |           |           |           |
%                easy coherence   |    hard coherence
%          |easy rule | hard rule | easy rule | hard rule |
% sem      |          |           |           |           |


% lab = categorical({'easy coh', 'hard coh', 'easy rule', 'hard rule'});
% lab = reordercats(lab,{'easy coh', 'hard coh', 'easy rule', 'hard rule'});
% bar = bar(lab,data_summary_all(3,:));
% bar.FaceColor = 'flat';
% bar.CData(1,:) = [.5 0 .5];
% bar.CData(2,:) = [.5 0 .10];
% bar.CData(3,:) = [0 0 255];
% bar.CData(4,:) = [.10 0 .5];
% hold on
% er = errorbar(lab,data_summary_all(3,:),sem_summary(1,:));    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% hold off
% % 
% labtwo = categorical({'easy coh easy rule', 'easy coh hard rule', 'hard coh easy rule', 'hard coh hard rule'});
% labtwo = reordercats(labtwo,{'easy coh easy rule', 'easy coh hard rule', 'hard coh easy rule', 'hard coh hard rule'});
% bar2 = bar(labtwo,data_summary_all(6,:));
% bar2.FaceColor = 'flat';
% bar2.CData(1,:) = [.5 0 .5];
% bar2.CData(2,:) = [.5 0 .10];
% bar2.CData(3,:) = [0 0 255];
% bar2.CData(4,:) = [.10 0 .5];
% hold on
% er2 = errorbar(labtwo,data_summary_all(6,:),sem_summary(2,:));    
% er2.Color = [0 0 0];                            
% er2.LineStyle = 'none';  
% hold off

% save all that
save(fullfile(analysisdir, 'mean_summary_total'), 'mean_summary');

%% get sem without nans
function semval = nansem(vector_data)
% Recall that s.e.m. = std(x)/sqrt(length(x));
nonan_std = nanstd(vector_data);
nonan_len = length(vector_data(~isnan(vector_data)));
% Plug in values
semval = nonan_std / sqrt(nonan_len);
end
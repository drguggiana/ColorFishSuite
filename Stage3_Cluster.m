%% Cluster groups of fish together
%USE ONLY 1 GROUP AT A TIME
%% clean up
clearvars
close all force
addpath(genpath('E:\Behavioral data\Matlab'))

%define whether to do regressor analysis
reg_var = 0;
%define whether to analyze the fourier results (assignment will happen
%anyway)
four_var = 0;
%define whether to do the gain analysis, even if it's the right program
%(i.e. cluster only)
gain_on = 0;
%% Load the files and define paths

%get the folder where the image files are
tar_path_all = uipickfiles('FilterSpec',...
    'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\Stage2_threshold');

%get the number of experiments selected
num_exp = length(tar_path_all);

%define the list of labels to sort the files
% label_list = {'_thres.mat','_thresBEH.mat'};
label_list = {'_thres.mat'};

%get the program used
[~,name_whole,~] = fileparts(tar_path_all{1});
name_parts = strsplit(name_whole,'_');
%start all the control variables off
gains_var = 0;
coneiso_var = 0;
exc_var = 0;
%if it's a gain program
if strcmp(name_parts{1},'p17b')==1
    gains_var = 1;
elseif strcmp(name_parts{1},'p15b')==1
    coneiso_var = 1;
elseif strcmp(name_parts{1},'p6')==1
    exc_var = 1;
end

%get the number of each type of data (the round is to avoid the nonscalar
%warning for the for loop)
num_data = round(num_exp./length(label_list));

%allocate memory for the different types of files
name_cell = cell(num_data,length(label_list));

%for the types of files
for f_type = 1:length(label_list)
    %get the coordinates of the file names
    name_map = ~cellfun(@isempty,strfind(tar_path_all,label_list{f_type}));
    %store them in the corresponding layer of the name cell
    name_cell(:,f_type) = tar_path_all(name_map);
end
%% Define/load constants

%get the number of time points
time_num = load(name_cell{1,1},'time_num');
time_num = time_num.time_num;
%get the number of stimuli
stim_num2 = load(name_cell{1,1},'stim_num2');
stim_num2 = stim_num2.stim_num2;
%and the color info
col_out = load(name_cell{1,1},'col_out');
col_out = col_out.col_out;

%define the pre and post periods
pre_time = false(time_num,1);
pre_time(6:0.25*time_num) = 1;
stim_time = false(time_num,1);
stim_time(0.25*time_num+1:0.75*time_num) = 1;
post_time = false(time_num,1);
post_time(0.75*time_num+1:end) = 1;

%path for the clustering BIC results
bic_name = strcat('E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\bic_files\',...
    'GroupFile_',num2str(num_data),'_fish');
%% Load the traces (prep for clustering)

%define a list of the stimuli to keep
% keep_stim = [1 3 13 14 10 12];
keep_stim = [];
%load and concatenate all the traces

%load the conc_trace variable into the trace_all cell
conc_trace = load(name_cell{1},'conc_trace');
%extract the matrix of values from the struct
conc_trace = conc_trace.conc_trace;

%load the fish origin information
fish_ori = load(name_cell{1},'fish_ori');
fish_ori = fish_ori.fish_ori;

%get the number of fish
fish_num = length(unique(fish_ori(:,1)));

if ~isempty(keep_stim)
    %load a temp matrix for calculations
    temp_mat = trace_all{fish};
    %allocate memory for the stimuli to keep
    keep_mat = zeros(size(temp_mat,1),time_num*length(keep_stim));
    
    %initialize a counter for indexing the matrix
    keepc = 1;
    %for all the stimuli to keep
    for kstim = keep_stim
        %store the desired stim in the matrix
        keep_mat(:,keepc:keepc+time_num-1) = ...
            temp_mat(:,(kstim-1)*time_num+1:kstim*time_num);
        %update the counter
        keepc = keepc + time_num;
    end
    %load the modified matrix into the original storage
    trace_all{fish} = keep_mat;
    
    %if it's the first fish
    if fish == 1
        %redefine the number of stimuli
        stim_num2 = length(keep_stim);
        %also modify the color matrix
        col_out = col_out(keep_stim,:,:);
    end
end

%get the number of traces
trace_num = size(conc_trace,1);
 
% %concatenate the responses
% conc_trace = vertcat(trace_all{:});
% %and create a vector with the fish of origin for each seed and the original
% %seed number
% fish_ori = zeros(size(conc_trace,1),2);
% %initialize a counter for the index within the seeds
% fish_count = 1;
% %for all the experiments
% for fish = 1:num_data
%     %get the number of seeds in the fish
%     seed_num = size(trace_all{fish},1);
%     
%     %insert it in the final vector
%     fish_ori(fish_count:fish_count+seed_num-1,1) = fish;
%     %also insert the number of each seed
%     fish_ori(fish_count:fish_count+seed_num-1,2) = 1:seed_num;
%     %update the counter
%     fish_count = seed_num + fish_count;
% end

% %flip the middle stimuli in conc_trace to match wavelength ONLY FOR P17b
% conc_trace = conc_trace(:,[1:80,161:240,81:160,241:320]);
%% OFF OPTIONAL Make the first 5 frames zero
%this is done for all files before 7/7 since the saving of the files at the
%end of each rep would trigger fast flashing of stimuli due to
%randomization. I added a delay period after a randomization of 10 seconds
%to avoid the issue

% %make a vector to kill 1 trial
% trial_vec = zeros(time_num,1);
% %kill the first 5 frames
% trial_vec(1:8) = 1;
% %rep the vector to cover the entire experiment
% trial_full = logical(repmat(trial_vec,stim_num2,1));
% %kill the frames in the actual data
% conc_trace(:,trial_full) = 0;
%% Exclude the outsides of the trial (pre and post stim)

%reshape the matrix for conversion
conc_trace2 = reshape(conc_trace,trace_num,time_num,[]);
%eliminate the undesired time fragments
conc_trace2 = conc_trace2(:,21:60,:);
%rewrite the time_num and conc_trace variables
time_num2 = size(conc_trace2,2);
conc_trace2 = reshape(conc_trace2,trace_num,[]);
%% Normalize by stimulus

% %reshape the matrix to make the stimuli more accessible
% conc_trace3 = reshape(conc_trace2,trace_num,time_num2,[]);
% 
% %for all the stimuli
% for stim = 1:stim_num2
%     %overwrite the values in conc_trace3
%     conc_trace3(:,:,stim) = normr_1(conc_trace3(:,:,stim),0);
% end
% 
% %overwrite conc_trace2
% conc_trace2 = reshape(conc_trace3,trace_num,[]);
%% Cluster the traces
tic
%define the sPCA parameters to use
bounds_top = 1:time_num2:size(conc_trace2,2);
bounds_bottom = [bounds_top(2:end)-1,size(conc_trace2,2)];
bounds = [bounds_top;bounds_bottom];
K = ones(1,stim_num2).*4;
t_bins = ones(stim_num2,1).*10;
pca_vec = ones(stim_num2,1).*1;

%define the vector of cluster numbers to try
% clu_vec = [2 3 4 5 10 20 30 40 50];
%     clu_vec = [100 200 300];  %MOST USED ONE
% clu_vec = 40;
% clu_vec = [18 20 22 24 26 28 30 32];
clu_vec = [5 10 20 30 50 70 100];
% clu_vec = [40 50 60];
% clu_vec = 30;
% clu_vec = [50 80 100 150 170];
replicates = 10;
[idx_clu,GMModel,clu_num,pcs,bic_vec] = sPCA_GMM_cluster_4(conc_trace2,bounds...
    ,K,t_bins,pca_vec,bic_name,clu_vec,replicates);
toc
%save the original cluster number for the surrogate analysis
ori_clu = clu_num;
%% Exclude clusters with snr below a threshold and with too few traces

%load the snr for the dataset
snr_mat = load(name_cell{1},'snr_mat');
snr_mat = snr_mat.snr_mat;

%calculate the average snr vector for each cluster

%allocate memory for the clusters
clu_snr = zeros(clu_num,size(snr_mat,2));

%for all the clusters
for clu = 1:clu_num
    
    %load the average matrix
    clu_snr(clu,:) = mean(snr_mat(idx_clu==clu,:),1);
end

% %plot a histogram of the averages
% figure
% for stim = 1:stim_num2
%     subplot(round(sqrt(stim_num2)),ceil(sqrt(stim_num2)),stim)
%     histogram(clu_snr(:,stim),100)
% end

%find a threshold by defining the 50th percentile on each stimulus
snr_thres = prctile(clu_snr,50,1);

%keep traces that pass the threshold in at least 1 stim
snr_thres_all = bsxfun(@gt,clu_snr,snr_thres);
snr_thres_vec = sum(snr_thres_all,2)>0;
%also modify the vector based on the number of traces per cluster
num_thres = 10;
%for all the clusters
for clu = 1:clu_num
    %kill the ones with less than num_thres traces
    if sum(idx_clu==clu)<num_thres
        snr_thres(idx_clu==clu) = 0;
        idx_clu(idx_clu==clu) = 0;
        snr_thres_vec(clu) = 0;
    end
end

%modify the idx_clu vector and clu_num to reflect this fact
%find all the clusters that stay
new_clu = find(snr_thres_vec);
%create a new idx_clu to overwrite the old one
new_idx = zeros(size(idx_clu));
%initialize a new cluster counter
c_count = 1;
%for all the clusters
for clu = 1:clu_num
    %check whether the cluster number stays
    if sum(new_clu==clu)>0
        %renumber all instances in the new idx vector
        new_idx(idx_clu==clu) = c_count;
        %update the counter
        c_count = c_count + 1;
        %clusters not found will be left as zeros
    end
end

%overwrite the old variables
idx_clu = new_idx;
clu_num = numel(unique(idx_clu))-1;
%% OFF Merge clusters that are too correlated
% close all
% % calculate the cluster averages
% %allocate memory for the averages
% clu_ave = zeros(clu_num,size(conc_trace,2));
% %and for the trace number
% clu_number = zeros(clu_num,1);
% %for all the clusters
% for clu = 1:clu_num
%     %calculate the cluster average
%     clu_ave(clu,:) = mean(conc_trace(idx_clu==clu,:),1);
%     %and store the number of traces going into each average
%     clu_number(clu) = sum(idx_clu==clu);
% end
% 
% 
% % calculate the correlation matrix across clusters
% [corr_matrix, pval] = corr(clu_ave');
% % corr_matrix = squareform(pdist(clu_ave));
% 
% % define the correlation threshold
% corr_threshold = 0.90;
% % define the pval threshold
% pval_threshold = 0.05;
% 
% figure
% subplot(1,2,1)
% imagesc(corr_matrix)
% subplot(1,2,2)
% imagesc(pval<pval_threshold&corr_matrix>corr_threshold)
% figure
% imagesc(normr_1(clu_ave,1))
% % produce a list of the correlated clusters
% [row, col] = find(tril(pval<pval_threshold&corr_matrix>corr_threshold,1));
% 
% % use the list to plot these clusters together
% figure
% % initialize a counter to lift the plots
% lift_count = 0;
% % get the increase interval
% lift_interval = 0.95*max(clu_ave(:));
% % for all the pairs
% for pairs = 1:length(row)
%     % plot the pairs together
%     plot(1:size(clu_ave,2), clu_ave(row(pairs),:) - lift_count,'k')
%     hold('on')
%     plot(1:size(clu_ave,2), clu_ave(col(pairs),:) - lift_count,'r')
%     % increase the counter
%     lift_count = lift_count + lift_interval;
% end
% 
% % merge the clusters by combining their indexes
% % for all the pairs
% for pairs = 1:length(row)
%     
% end
% 
% figure
% histogram(tril(corr_matrix,1))
%% OFF OLD FOR P6 FIGURE ONLY

% if exc_var == 1
%     exc_stim = [1:4,6,8,11,12,13,15];
%     
%     %get the new number of stimuli
%     new_stim = stim_num2-length(exc_stim);
%     %allocate memory for the new stim matrix
%     exc_trace = zeros(size(conc_trace,1),time_num.*new_stim);
%     %and for the new color code matrix
%     new_col = zeros(new_stim,time_num,8);
%     %initialize a time trace stim counter (to index the time trace)
%     stim_c = 1;
%     
%     %for all the stimuli
%     for stim = 1:stim_num2
%         %if the stimulus is not in the exclusion matrix
%         if ~any(stim==exc_stim)
%             %add it to the new matrix
%             exc_trace(:,stim_c:stim_c+time_num-1) = conc_trace(:,1+time_num*(stim-1):time_num*stim);
%             %and add the color code
%             new_col((stim_c-1+time_num)/time_num,:,:) = col_out(stim,:,:);
%             %update the counter
%             stim_c = stim_c + time_num;
%         end
%     end
%     
%     % figure
%     % imagesc(conc_trace)
%     % figure
%     % imagesc(exc_trace)
%     
%     %rewrite the conc_trace and stim_num2 variables
%     stim_num2 = new_stim;
%     conc_trace = exc_trace;
%     %also modify the color code variable
%     col_out = new_col;
% end
%% Save analysis output

%define the save path
save_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\Stage3_cluster\';

save_var = 1;
if save_var == 1

%     %get the root of the save name
%     [ori_name,~] = uiputfile(strcat(save_path,'*.*'));
    % get the original file name
    ori_name = name_whole(1:end-6);
    %save the clustering output
    save_clu = strcat(ori_name,'_clusters.mat');
    save(fullfile(save_path,save_clu),'pcs','GMModel','idx_clu','clu_num',...
        'conc_trace','bounds','K','t_bins','pca_vec','stim_num2','time_num',...
        'col_out','bic_vec','fish_ori')
end

load_var = 0;

if load_var == 1
    [f_name,f_path] = uigetfile(save_path);
    load(fullfile(f_path,f_name))
end
%% Cluster groups of fish together
%USE ONLY 1 GROUP AT A TIME
%% clean up
clearvars
close all force
addpath(genpath('E:\Behavioral data\Matlab'))
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
%% Load the traces (prep for clustering)

% for all the files
for files = 1:num_data
    
    %% Define/load constants
    % define the name to save the files
    [~,ori_name] = fileparts(name_cell{files,1});
    ori_name = ori_name(1:end-6);
    %path for the clustering BIC results
    bic_name = strcat('E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\bic_files\',...
    'GroupFile_',num2str(num_data),'_fish');
    %get the number of time points
    time_num = load(name_cell{files,1},'time_num');
    time_num = time_num.time_num;
    %get the number of stimuli
    stim_num2 = load(name_cell{files,1},'stim_num2');
    stim_num2 = stim_num2.stim_num2;
    
    %define the pre and post periods
    pre_time = false(time_num,1);
    pre_time(6:0.25*time_num) = 1;
    stim_time = false(time_num,1);
    stim_time(0.25*time_num+1:0.75*time_num) = 1;
    post_time = false(time_num,1);
    post_time(0.75*time_num+1:end) = 1;
    % load the color info
    col_out = load(name_cell{files,1},'col_out');
    col_out = col_out.col_out;

    %define a list of the stimuli to keep
    % keep_stim = [1 3 13 14 10 12];
    keep_stim = [];
    %load and concatenate all the traces

    %load the conc_trace variable into the trace_all cell
    conc_trace = load(name_cell{files},'conc_trace');
    %extract the matrix of values from the struct
    conc_trace = conc_trace.conc_trace;

    %load the fish origin information
    fish_ori = load(name_cell{files},'fish_ori');
    fish_ori = fish_ori.fish_ori;

    % load the anatomy info
    anatomy_info = load(name_cell{files},'cat_anatomy_all');
    anatomy_info = anatomy_info.cat_anatomy_all;
    
    % load the single reps
    single_reps = load(name_cell{files},'cat_reps');
    single_reps = single_reps.cat_reps;

    %get the number of fish
    fish_num = length(unique(fish_ori(:,1)));

    %get the number of traces
    trace_num = size(conc_trace,1);
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
    %% OFF Normalize by stimulus

%     %reshape the matrix to make the stimuli more accessible
%     conc_trace3 = reshape(conc_trace2,trace_num,time_num2,[]);
%     
%     %for all the stimuli
%     for stim = 1:stim_num2
%         %overwrite the values in conc_trace3
%         conc_trace3(:,:,stim) = normr_1(conc_trace3(:,:,stim),0);
%     end
%     
%     %overwrite conc_trace2
%     conc_trace2 = reshape(conc_trace3,trace_num,[]);
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
%     clu_vec = 30;
    % clu_vec = [50 80 100 150 170];
    replicates = 10;
    [idx_clu,GMModel,clu_num,pcs,bic_vec] = sPCA_GMM_cluster_Color(conc_trace2,bounds...
        ,K,t_bins,pca_vec,bic_name,clu_vec,replicates);
    toc
    %save the original cluster number for the surrogate analysis
    ori_clu = clu_num;
    %% Exclude clusters with snr below a threshold and with too few traces

    %load the snr for the dataset
    snr_mat = load(name_cell{files},'snr_mat');
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
    
    %% Calculate the cluster average
    
    %allocate memory for the averages
    clu_ave = zeros(clu_num,size(conc_trace,2));
    %and for the trace number
    clu_number = zeros(clu_num,1);
    %for all the clusters
    for clu = 1:clu_num
        %calculate the cluster average
        clu_ave(clu,:) = mean(conc_trace(idx_clu==clu,:),1);
        %and store the number of traces going into each average
        clu_number(clu) = sum(idx_clu==clu);
    end
    %% Load the info about seeds
    
    z_seed = load(name_cell{files,1},'cat_z_all');
    z_seed = z_seed.cat_z_all;
    
    xy_seed = load(name_cell{files,1},'cat_seed_all');
    xy_seed = xy_seed.cat_seed_all;
    
    ave_stack = load(name_cell{files,1},'cat_stack_all');
    ave_stack = ave_stack.cat_stack_all;
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


%     %get the root of the save name
%     [ori_name,~] = uiputfile(strcat(save_path,'*.*'));
    % get the original file name
%     ori_name = name_whole(1:end-6);
    %save the clustering output
    save_clu = strcat(ori_name,'_clusters.mat');
%     save(fullfile(save_path,save_clu),'pcs','GMModel','idx_clu','clu_num',...
%         'conc_trace','bounds','K','t_bins','pca_vec','stim_num2','time_num',...
%         'col_out','bic_vec','fish_ori','anatomy_info','single_reps')
    % assemble a structure with all of the data
    main_str = struct([]);
    main_str(1).name = ori_name;
    main_str(1).conc_trace = conc_trace;
    main_str(1).single_reps = single_reps;
    main_str(1).col_out = col_out;
    main_str(1).fish_ori = fish_ori;
    main_str(1).anatomy_info = anatomy_info;
    main_str(1).pcs = pcs;
    main_str(1).GMModel = GMModel;
    main_str(1).idx_clu = idx_clu;
    main_str(1).bounds = bounds;
    main_str(1).K = K;
    main_str(1).t_bins = t_bins;
    main_str(1).pca_vec = pca_vec;
    main_str(1).stim_num = stim_num2;
    main_str(1).time_num = time_num;
    main_str(1).bic_vec = bic_vec;
    main_str(1).clu_num = clu_num;
    main_str(1).clu_number = clu_number;
    main_str(1).clu_ave = clu_ave;
    main_str(1).z_seed = z_seed;
    main_str(1).xy_seed = xy_seed;
    main_str(1).ave_stack = ave_stack;
    main_str(1).snr_mat = snr_mat;
    
    % save the structure
    save(fullfile(save_path,save_clu),'main_str')
    
end
%% Cluster groups of fish together

%% clean up
clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))
%define whether to save the output or not
save_var = 1;
%define the fish to combine
% fish_comb = [1 2 3 4 5 6 7 8 9 10;1 1 2 2 3 4 4 5 5 6];
% fish_comb = [1 2;1 1];
%% Load the files and define paths

% set the rng
rng(1)

%get the folder where the image files are
tar_path_all = uipickfiles('FilterSpec',strcat(paths(1).stage1_path,'*.mat'));

%get the number of experiments selected
num_exp = length(tar_path_all);

% %if fish_comb wasn't specified, treat every file as an independent fish
% if ~exist('fish_comb','var')
%     fish_comb = [1:num_exp;1:num_exp];
% end
%define the list of labels to sort the files
label_list = {'_traces.mat','_plot.mat'};

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
time_num = load(name_cell{1,2},'time_num');
time_num = time_num.time_num;
%get the number of stimuli
stim_num2 = load(name_cell{1,2},'stim_num2');
stim_num2 = stim_num2.stim_num2;
%and the color info
col_out = load(name_cell{1,2},'col_out');
col_out = col_out.col_out;

%define the pre and post periods
pre_time = false(time_num,1);
pre_time(6:0.25*time_num) = 1;
stim_time = false(time_num,1);
stim_time(0.25*time_num+1:0.75*time_num) = 1;
post_time = false(time_num,1);
post_time(0.75*time_num+1:end) = 1;

% load the thresholding constants
snr_param = load(paths(1).param_path,'snr_param');
snr_param = snr_param.snr_param;

% select the corresponding thresholds based on the dataset
% decompose the name of the file
[~,fname] = fileparts(name_cell{1,1});
fname = strsplit(fname,'_');

% get the logical vectors to select the apropriate field
logic_vec = contains({snr_param(:).strain},fname{2}) & contains({snr_param(:).program},fname{3});
% get the values
perc_thres = snr_param(logic_vec).percentile;
stim_thres = snr_param(logic_vec).stimThreshold;
%% Get the snr threshold

%allocate memory to store the snr for the entire fish
snr_all = cell(num_data,1);
%for all the files
for files = 1:num_data
    %load and concatenate the snr for this file
    snr_temp = load(name_cell{files,1},'snr_cell');
    snr_all{files} = cat(1,snr_temp.snr_cell{:});
end
%concatenate the snr for all of the files
snr_mat = cat(1,snr_all{:});

% close all
% figure
% histogram(snr_mat(:),1000)
% 
% figure
% %for all of the stimuli
% for stim  = 1:stim_num2
%     subplot(round(sqrt(stim_num2)),ceil(sqrt(stim_num2)),stim)
%     histogram(snr_mat(:,stim),100)
% end

%allocate memory for storing the thresholds
thres_vec = zeros(stim_num2,1);

%for all of the stimuli
for stim = 1:stim_num2
    %calculate the 25th percentile and store as the threshold
    thres_vec(stim) = prctile(snr_mat(:,stim),perc_thres);
end

% figure
% plot(thres_vec)
%% Load the traces 

%define a list of the stimuli to keep
% keep_stim = [1 3 13 14 10 12];
keep_stim = [];
%load and concatenate all the traces

%allocate memory for the actual traces
trace_all = cell(num_data,1);
%for all the traces
for fish = 1:num_data
    %load the conc_trace variable into the trace_all cell
    trace_all{fish} = load(name_cell{fish,1},'conc_trace');
    %extract the matrix of values from the struct
    trace_all{fish} = trace_all{fish}.conc_trace;
    
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
end

%concatenate the responses
conc_trace = vertcat(trace_all{:});
%and create a vector with the fish of origin for each seed and the original
%seed number
fish_ori = zeros(size(conc_trace,1),2);
%initialize a counter for the index within the seeds
fish_count = 1;
%for all the experiments
for fish = 1:num_data
    %get the number of seeds in the fish
    seed_num = size(trace_all{fish},1);
    
    %insert it in the final vector
    fish_ori(fish_count:fish_count+seed_num-1,1) = fish;
    %also insert the number of each seed
    fish_ori(fish_count:fish_count+seed_num-1,2) = 1:seed_num;
    %update the counter
    fish_count = seed_num + fish_count;
end


% mix fish that are separated in 2 experiments
% allocate memory for the combination vector
fish_comb = zeros(num_data,2);
% parse the names
for fish = 1:num_data
    % get the filename
    [~,name,~] = fileparts(name_cell{fish,1});
    % split the filename
    name_parts = strsplit(name,'_');
    % if there's b in the number, put it in the second slot, otherwise the
    % first one
    if contains(name_parts{4},'b')
        fish_comb(fish,2) = fish;
    else
        fish_comb(fish,1) = fish;
    end
end
% get the number of fish
num_fish = sum(fish_comb(:,1)>0);
%modify fish ori to combine the fish on 2 different stacks based on the
%fish_comb vector

% allocate a copy of fish ori
fish_ori_copy = zeros(size(fish_ori));
fish_ori_copy(:,2) = fish_ori(:,2);
% initialize a counter to move through fish_comb
fish_counter = 1;
%for all the cols in fish_comb
% for fcol = 1:size(fish_comb,2)
for fish = 1:num_fish
    % replace the indexes of this fish
    fish_ori_copy(fish_ori(:,1)==fish_comb(fish_counter,1),1) = fish;
    % increment the counter
    fish_counter = fish_counter + 1;
    % if the next position exists and is populated, also use that indexing
    if size(fish_comb,1) >= fish_counter && fish_comb(fish_counter,2) > 0
        fish_ori_copy(fish_ori(:,1)==fish_comb(fish_counter,2),1) = fish;
        % increment the counter
        fish_counter = fish_counter + 1;
    end
%     %if the previous file also was from the same fish
%     if fcol>1 && fish_comb(2,fcol)==fish_comb(2,fcol-1)
%         curr_coord = fish_ori(:,1)==fish_comb(1,fcol);
%         fish_ori(curr_coord,2) = fish_ori(curr_coord,2) + file_c;
%     end
%     fish_ori(fish_ori(:,1)==fish_comb(1,fcol),1) = fish_comb(2,fcol);
%     file_c = max(fish_ori(fish_ori(:,1)==fish_comb(2,fcol),2));
end

% replace fish_ori
fish_ori = fish_ori_copy;

%detect if the stimulus is p17b. if it is, reorder the stimuli in the
%middle due to the projector switch
%get the program used
% [~,name_whole,~] = fileparts(tar_path_all{1});
% name_parts = strsplit(name_whole,'_');

% save the original stimulus number and file name
ori_fname = fname;
ori_stim_num2 = stim_num2;
% correct dataset specific mistakes
[conc_trace, stim_num2, col_out, fname] = dataset_fixes(fname,conc_trace,stim_num2,col_out,time_num);
%% Threshold the traces based on the snr calculation above

%the idea is to exclude any trace that has a value under the threshold,
%under any stimulus

%turn the NaNs (probably from traces with zeros across the board)into 0
snr_mat(isnan(snr_mat)) = 0;
%get the binary version of the traces
snr_bin = bsxfun(@gt,snr_mat,thres_vec');

thres_res = zeros(stim_num2,1);
for t = 1:stim_num2
    %to do this, AND the columns of the snr matrix
    snr_vec = logical(sum(snr_bin,2)<t);
    %get the number of traces excluded
    thres_res(t) = sum(snr_vec);
end

figure
plot(thres_res./length(snr_vec));
xlabel('Number of stimuli with significant signal required')
ylabel('Proportion of traces excluded')

%based on the results above, I'll set the threshold to be at least 10
%stimuli with significant signal, which excludes around 11% of the traces
snr_vec = logical(sum(snr_bin,2)>stim_thres);

% also remove the traces that have too high a signal (only a couple)
% define the dfof threshold
max_threshold = 100;
% identify these traces first
hightrace_idx = find(max(conc_trace,[],2)>max_threshold);
% modify the snr_vector to exclude these traces too
snr_vec(hightrace_idx) = 0;
% report the number of traces excluded
fprintf(strjoin({'Traces excluded because of max threshold:',...
    num2str(length(hightrace_idx)),'\r\n'},'_'))

%adapt the conc_trace and fish_ori matrices 
conc_trace = conc_trace(snr_vec,:);
fish_ori = fish_ori(snr_vec,:);

%allocate memory for concatenation of the matrices
cat_seed = cell(num_data,1);
cat_z = cell(num_data,1);
cat_stack = cell(num_data,1);
% allocate a cell for the anatomy info
cat_anatomy = cell(num_data,1);
% allocate a cell for the reps
cat_reps = cell(num_data,1);

%for all the fish
for fish = 1:num_data
    %load the seed xy information
    seed_concat = load(name_cell{fish,1},'seed_concat');
    cat_seed{fish} = seed_concat.seed_concat;
    
    %load the z position of each seed
    z_seed = load(name_cell{fish,1},'z_seed');
    
    % load the anatomy
    anatomy_info = load(name_cell{fish,1},'anatomy_info');
    anatomy_info = anatomy_info.anatomy_info;
    
    % load the single reps
    all_trace_reps = load(name_cell{fish,1},'all_trace_reps');
    all_trace_reps = all_trace_reps.all_trace_reps;
    % get the number of traces
    trace_num = size(all_trace_reps,1);
    % reshape to have stimuli in one dimension
    all_trace_reps = reshape(all_trace_reps,trace_num,time_num,stim_num2,[]);
    % fix the dataset specific issues
    [all_trace_reps, ~, ~, ~] = dataset_fixes(ori_fname,all_trace_reps,ori_stim_num2,[],time_num);

    
    %if it the first fish
    if fish == 1
        %just copy the z info
        cat_z{fish} = z_seed.z_seed;
        z_max = 0;
    else %otherwise, add the last z from the previous set to the current
        %so they match the ave_stack
        cat_z{fish} = z_max + z_seed.z_seed; 
    end
    %load the ave_stack of each seed
    ave_stack = load(name_cell{fish,1},'ave_stack');
    cat_stack{fish} = ave_stack.ave_stack;
    %rewrite the z_max variable for next iteration
    z_max = size(cat_stack{fish},3) + z_max;
    
    % if anatomy_info is empty, replace it with NaNs
    if isempty(anatomy_info)
        anatomy_info = nan(size(seed_concat.seed_concat,1),2);
    end
    cat_anatomy{fish} = anatomy_info;
    cat_reps{fish} = all_trace_reps;
end

%generate concatenated versions   
cat_seed_all = cat(1,cat_seed{:});
cat_z_all = cat(1,cat_z{:});
cat_stack_all = cat(3,cat_stack{:});
cat_anatomy_all = cat(1,cat_anatomy{:});
cat_reps = cat(1,cat_reps{:});

%filter them based on the threshold
cat_seed_all = cat_seed_all(snr_vec,:);
cat_z_all = cat_z_all(snr_vec,:);
% if the anatomy is empty, keep it empty
if ~isempty(cat_anatomy_all)
   cat_anatomy_all = cat_anatomy_all(snr_vec,:);
end
cat_reps = cat_reps(snr_vec,:);


%define the path for saving the concatenated files
thres_path = paths(1).stage2_path;
%define the saving path
%     [thres_name,thres_fpath] = uiputfile(thres_path);

thres_name = strcat(fname{3},'_',fname{2});
%save the corrected variables into a new file, already grouped
save_name = strcat(thres_name,'_thres.mat');
save(fullfile(thres_path,save_name),'conc_trace','fish_ori','cat_stack_all',...
    'cat_seed_all','cat_z_all','time_num','stim_num2','col_out','snr_mat','cat_anatomy_all','cat_reps')
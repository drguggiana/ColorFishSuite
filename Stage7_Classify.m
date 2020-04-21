% Run a classifier through the data

clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'Classify\');

data = load_clusters(cluster_path);
%% Train a classifier for each region

% define whether to run in parallel (need to also convert region for loop
% to parfor)
run_parallel = 0;
% define whether to classify on clusters or regions
cluster_flag = 0;
% define the region set to use
region_set = 1;
% define which regions to include in the analysis
switch region_set
    case 1
        tectum_regions = {'R-TcN','R-TcP'};
        af_regions = {'AF10'};
    case 2
        tectum_regions = {'R-TcN','R-TcP','R-Cb','R-Hb','R-Pt','L-TcN','L-TcP','L-Cb','L-Hb','L-Pt'};
        af_regions = {'AF4','AF5','AF8','AF9','AF10'};
end


% define the number of repeats to run each classification scheme
repeat_number = 10;
% define whether to subsample. If subsampling from each region within a
% dataset, use 1, if subsampling across datasets (usually with
% combination=1) then use 2r
subsample = 2;
% combine regions
region_combination = 1;
%define whether to shuffle labels (for neutral classification)
shuff_label = 0;
%define the number of classes per color (1,3,5,or 8) (or 10,11 and 12 for the
%p6p8 data)
classpcolor = 16;
%define the binning factor
bin_width = 10;
% define which portion of the trial to take (0 pre, 1 stim, 2 post,
% 3 pre and post)
portion = 1;

% define the subsampling constant for the data
sub_constant = 0.95;

tic

% set leave one out cross validation
loo = 1;
% set the fraction of traces
trace_frac_mult = 1;
% set the partitions
set_part = 1;
% set the repetition number (not used since repeating outside the
% classification function)
redec_num = 1;

% get the number of datasets
num_data = length(data);

% process the regions first
% allocate memory to store the region information
region_cell = cell(num_data,2);
% for all data sets
for datas = 1:num_data
    % load the anatomy info
    anatomy_info = data(datas).anatomy_info(:,1);
    
    % define the regions to be considered (depending on the stimulus protocol)
    if contains(data(datas).name,{'syn','Syn'})
        region_list = af_regions;
    else
        region_list = tectum_regions;
    end
    % separate the traces by region or cluster
    if cluster_flag
        [region_cell{datas,1},region_cell{datas,2}] = region_split(data(datas).single_reps,...
            data(datas).idx_clu,data(datas).name,region_combination,[],cluster_flag);
    else
        [region_cell{datas,1},region_cell{datas,2}] = region_split(data(datas).single_reps,...
            anatomy_info,data(datas).name,region_combination,region_list,cluster_flag);
    end
end

% if subsample is 2, subsample across datasets
if subsample == 2
    % allocate memory for all datasets
    subsample_vector = zeros(num_data,1);
    % for all datasets
    for datas = 1:num_data
        subsample_vector(datas) = min(cellfun(@size, region_cell{datas,1}(:,1), ...
            num2cell(ones(region_cell{datas,2},1))));
    end
    % take the overall min as the subsample number
    number_subsample = round(min(subsample_vector)*sub_constant);
end

% for all the data sets
for datas = 1:num_data
    % allocate memory for the output
    main_str = struct([]);
    % load the region info
    region_data = region_cell{datas,1};
    num_regions = region_cell{datas,2};
    % if subsample is on, determine how many traces to take per rep
    if subsample == 1
        number_subsample = round(min(cellfun(@size, region_data(:,1), ...
            num2cell(ones(num_regions,1))))*sub_constant);
    end
    % allocate memory
    class_cell = cell(num_regions,1);
    % get the numbers of time, stim and reps
    time_num = data(datas).time_num;
    stim_num = data(datas).stim_num;
    rep_num = size(region_data{1,1},2)/time_num/stim_num;
    % get the number of animals and the animal info
    fish_ori_all = data(datas).fish_ori;
    num_animals = length(unique(fish_ori_all(:,1)));
    
    % create the pool of workers
    if datas == 1 && run_parallel==1 && exist('worker_pool','var') == 0
        worker_pool = parpool;
    end
    % get just the coordinates to not broadcast the whole variable
    region_coord = region_data(:,3);
    % set the same fraction for all regions
    trace_frac = ones(size(region_data,1),1).*trace_frac_mult;

    
    % for all the regions
    for regions = 1:num_regions
        if isempty(region_data{regions,1})
            continue
        end
        % allocate a cell to store each repetition
        class_repeats = cell(3,repeat_number);
        % for all of the repeats
        for repeat = 1:repeat_number
            
            % get the traces
            region_traces = region_data{regions,1};
            % get the origin fish for these region traces
            region_ori = fish_ori_all(region_coord{regions}==1,1);
            % if subsampling, implement
            if subsample > 0
                % generate a random indexing vector with the desired number
                % of traces
                subsample_idx = randperm(size(region_traces,1),number_subsample);
                % take the determined number of traces from the total
                region_traces = region_traces(subsample_idx,:);
                region_ori = region_ori(subsample_idx); 
            end
            % get the number of fish present and their ID
            fish_list = unique(region_ori);
            num_fish = length(fish_list);
            % allocate memory to store the classifier results
            fish_classifiers = cell(3,num_fish);
            % for all the fish
            for fish = 1:num_fish

                % get the traces for this fish and reshape for normalization
                temp_stimuli = reshape(region_traces(region_ori==fish_list(fish),:),[],time_num*stim_num,rep_num);

                % for all the reps
                for reps = 1:rep_num
                    temp_stimuli(:,:,reps) = normr_1(temp_stimuli(:,:,reps),0)-0.5;
                    temp_stimuli(:,:,reps) = temp_stimuli(:,:,reps)./std(temp_stimuli(:,:,reps),0,2);
                    % zscore per neuron
%                     temp_stimuli(:,:,reps) = zscore(temp_stimuli(:,:,reps),0,2);
%                     temp_stimuli(:,:,reps) = (temp_stimuli(:,:,reps))./std(temp_stimuli(:,:,reps),0,2);
                end

                % reshape back
                temp_stimuli = reshape(temp_stimuli,size(temp_stimuli,1),[]);


                % run the classifier
                [fish_classifiers(:,fish)] = clss_things_color(temp_stimuli,loo,trace_frac(regions),...
                    set_part,redec_num,shuff_label,classpcolor,bin_width,stim_num,rep_num,portion);
            end
            % allocate memory for the average classifier
            average_classifier = cell(size(fish_classifiers,1),1);
    
        
            % for all the fields in the cell
            for field = 1:size(fish_classifiers,1)
                % average across the fish classifiers
                average_classifier{field} = squeeze(mean(cat(4,fish_classifiers{field,:}),4));
            end
            class_repeats(:,repeat) = average_classifier;
        end

        % allocate memory for the average classifier
        average_repeats = cell(size(class_repeats,1),1);

        % for all the fields in the cell
        for field = 1:size(class_repeats,1)
            % average across the fish classifiers
            average_repeats{field} = squeeze(cat(4,class_repeats{field,:}));
        end
        class_cell{regions} = average_repeats;
        
    end
    
    %eliminate the areas with empty results
    empty_vec = cellfun(@isempty,class_cell);
    class_cell = class_cell(~empty_vec);
    
    % store the region info and the classifier results
    main_str(1).region = region_data;
    main_str(1).class = class_cell;
    main_str(1).num_regions = num_regions;
    main_str(1).empty_vec = empty_vec;
    % also store the basic info
    main_str(1).name = data(datas).name;
    main_str(1).classpcolor = classpcolor;
    main_str(1).subsample = subsample;
    main_str(1).loo = loo;
    main_str(1).shuff_label = shuff_label;
    main_str(1).repeat_number = repeat_number;
    main_str(1).bin_width = bin_width;
    main_str(1).redec_num = redec_num;
    main_str(1).portion = portion;
    main_str(1).cluster_flag = cluster_flag;
    %% Save the classifier
    
    % assemble the file name
    file_name = strjoin({'class',data(datas).name,'classp',num2str(classpcolor),...
        'subsample',num2str(subsample),'loo',num2str(loo),'shuff',num2str(shuff_label),...
        'reps',num2str(repeat_number),'bin',num2str(bin_width),'portion',num2str(portion),...
        'cluster',num2str(cluster_flag),'.mat'},'_');
    % save the file
    save(fullfile(paths.classifier_path,file_name),'main_str')
    %% Plot the results

    Stage7b_plotClassifier
end

% delete the pool of workers if it exists
if exist('worker_pool', 'var')
    delete(worker_pool)
end

toc
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
run_parallel = 1;
% define which regions to include in the analysis
tectum_regions = {'L-TcN','R-TcN','L-TcP','R-TcP','L-Pt','R-Pt'};
af_regions = {'AF9','AF10'};
% tectum_regions = {'L-TcN','R-TcN','L-TcP','R-TcP'};
% af_regions = {'AF10'};
% define the number of repeats to run each classification scheme
repeat_number = 10;
% define whether to subsample
subsample = 1;
% combine regions
region_combination = 1;
%define whether to shuffle labels (for neutral classification)
shuff_label = 1;
%define the number of classes per color (1,3,5,or 8) (or 10,11 and 12 for the
%p6p8 data)
classpcolor = 10;
%define the binning factor
bin_width = 10;

tic
% allocate memory for both data sets
data_struct = struct([]);
% set leave one out cross validation
loo = 1;
% set the fraction of traces
trace_frac_mult = 1;
% set the partitions
set_part = 1;
% set the repetition number (not used since repeating outside the
% classification function)
redec_num = 1;
% define which portion of the trial to take (0 pre, 1 stim, 2 post,
% 3 pre and post)
portion = 1;

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
    if contains(data(datas).name,'syn')
        region_list = af_regions;
    else
        region_list = tectum_regions;
    end
    % separate the traces by region
    [region_cell{datas,1},region_cell{datas,2}] = region_split(data(datas).single_reps,...
        anatomy_info,data(datas).name,region_combination,region_list);
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
    number_subsample = min(subsample_vector);
end

% for all the data sets
for datas = 1:num_data
    
    % load the region info
    region_data = region_cell{datas,1};
    num_regions = region_cell{datas,2};
    % if subsample is on, determine how many traces to take per rep
    if subsample == 1
        number_subsample = min(cellfun(@size, region_data(:,1), ...
            num2cell(ones(num_regions,1))));
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
    parfor regions = 1:num_regions
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
                    temp_stimuli(:,:,reps) = normr_1(temp_stimuli(:,:,reps),0);
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
    data_struct(datas).region = region_data;
    data_struct(datas).class = class_cell;
    data_struct(datas).num_regions = num_regions;
    data_struct(datas).empty_vec = empty_vec;    
    
end

% delete the pool of workers if it exists
if exist('worker_pool', 'var')
    delete(worker_pool)
end

toc
%% Plot the classifier results

close all

% for all the data sets
for datas = 1:length(data)
    % load the variables of interest
    class_cell = data_struct(datas).class;
    num_regions = data_struct(datas).num_regions;
    reg_label = {data_struct(datas).region{:,2}};
    empty_vec = data_struct(datas).empty_vec;
    %get the number of non-empty classifiers
    reg_full = size(class_cell,1);
    
    %plot the confusion matrices
    figure
    %for all the regions
    %initialize a plot counter
    p_count = 1;
    for regs = 1:num_regions
        if empty_vec(regs) == 1
            continue
        end
        subplot(round(sqrt(reg_full)),ceil(sqrt(reg_full)),p_count)
        if redec_num == 1 && repeat_number == 1
            imagesc(class_cell{p_count}{1})
        else
            imagesc(mean(class_cell{p_count}{1},3))
        end
        % calculate the classifier accuracy
        acc = num2str(round(mean(class_cell{p_count}{2}),2));
        
        axis square
        title(strjoin({reg_label{regs}, 'A', acc},'_'),'Interpreter','None')
        % set the color scale to the max number of trials per category
        set(gca,'CLim',[0, sum(class_cell{p_count}{1}(:,1))])
        %update the counter
        p_count = p_count + 1;
    end
    % assemble the figure path 
    file_path = strjoin({'confMatrix',data(datas).name,'classp',num2str(classpcolor),...
        'subsample',num2str(subsample),'loo',num2str(loo),'shuff',num2str(shuff_label),...
        'reps',num2str(repeat_number),'bin',num2str(bin_width),'.png'},'_');
    saveas(gcf, fullfile(fig_path,file_path), 'png')
    
    %plot the diagonal values
    figure
    %initialize a plot counter
    p_count = 1;
    %for all the regions
    for regs = 1:num_regions
        if empty_vec(regs) == 1
            continue
        end
        plot(p_count,class_cell{p_count}{2},'ok')
        hold('on')
        % calculate the exact accuracy
        exact_acc = mean(class_cell{p_count}{2});
        % if there were repeats
        if redec_num > 1 || repeat_number > 1
            plot(p_count,exact_acc,'or')
        end
        p_count = p_count + 1;
    end
    set(gca,'XTick',1:reg_full,'XTickLabels',reg_label(~empty_vec),'FontSize',20,...
        'XTickLabelRotation',45,'XLim',[0 reg_full+1])
    ylabel('Classifier Accuracy (a.u.)','FontSize',20)
    set(gca,'YLim',[0 1])
    % assemble the figure path 
    file_path = strjoin({'classAccuracy',data(datas).name,'classp',num2str(classpcolor),...
        'subsample',num2str(subsample),'loo',num2str(loo),'shuff',num2str(shuff_label),...
        'reps',num2str(repeat_number),'bin',num2str(bin_width),'.png'},'_');
    saveas(gcf, fullfile(fig_path,file_path), 'png')
    
%     %plot the diagonal values, normalized by number of traces
%     figure
%     %initialize a plot counter
%     p_count = 1;
%     %for all the regions
%     for regs = 1:num_regions
%         if empty_vec(regs) == 1
%             continue
%         end
%         plot(p_count,class_cell{p_count}{2}/size(data_struct(datas).region{regs},1),'ok')
%         hold('on')
%         p_count = p_count + 1;
%     end
%     set(gca,'XTick',1:reg_full,'XTickLabels',reg_label(~empty_vec),'FontSize',20,...
%         'XTickLabelRotation',45,'XLim',[0 reg_full+1])
%     ylabel('Cluster Norm. Classifier Accuracy (a.u.)','FontSize',20)
    
%     %plot the diagonal values, normalized by the area
%     figure
%     %initialize a plot counter
%     p_count = 1;
%     %for all the regions
%     for regs = 1:num_regions
%         if empty_vec(regs) == 1
%             continue
%         end
%         plot(p_count,class_cell{p_count}{2}/area_all(regs),'ok')
%         hold('on')
%         p_count = p_count + 1;
%     end
%     set(gca,'XTick',1:reg_full,'XTickLabels',reg_label(~empty_vec),'FontSize',20,...
%         'XTickLabelRotation',45,'XLim',[0 reg_full+1])
%     ylabel('Area Norm. Classifier Accuracy (a.u.)','FontSize',20)
end
autoArrangeFigures
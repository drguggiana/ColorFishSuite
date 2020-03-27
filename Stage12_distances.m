% Calculate and plot metrics about location of the seeds in the brain
%% clean up
clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'Distances\');

data = load_clusters(cluster_path);
%% Define which sorting index to use
% get the number of datasets
num_data = size(data,2);
% define the index
sorting_constant = 2;
%define the number of bins
n_bins = 20;
% define the number of random reps
randomization = 100;
% allocate memory for the data
sorting_cell = cell(num_data,3);

switch sorting_constant
    case 1
        % for all of the datasets
        for datas = 1:num_data
            sorting_cell{datas,1} = data(datas).idx_clu;
            sorting_cell{datas,2} = data(datas).clu_num;
            sorting_cell{datas,3} = data(datas).clu_number;
        end
    case 2 % raw color classes
        % for all of the datasets
        for datas = 1:num_data
            % get the gains
            gains = data(datas).delta_norm;
            % get the color with the max for each seed
            [~,sorting_cell{datas,1}] = max(abs(gains),[],2);
            % set the number of groups
            sorting_cell{datas,2} = 4;
            % get the number of seeds per group
            numbers = zeros(size(gains,1),1);
            % for all the groups
            for group = 1:4
                numbers(group) = sum(sorting_cell{datas,1}==group);
            end
            sorting_cell{datas,3} = numbers;
        end
    case 3 % gain clusters

end

% index out the non-tectal and non-AF10 traces
% for all the datasets
for datas = 1:num_data
    region_combination = 1;
    
    % define which regions to keep depending on the dataset
    if contains(data(datas).name, {'Syn','syn'})
        region_list = {'AF10'};
    else
        region_list = {'R-TcN','R-TcP'};
    end
        % load the anatomy info
    anatomy_info = data(datas).anatomy_info(:,1);
    
    % separate the traces by region
    [region_cell,~] = region_split(data(datas).single_reps,...
        anatomy_info,data(datas).name,region_combination,region_list);
    % rewrite the index vector
    sorting_cell{datas,1}(region_cell{3}==0) = 0;
    
end
%% Get the distance distributions for each cluster


%allocate memory for the results
files_dist = cell(num_data,1);

%for all the files
for datas = 1:num_data

    %and the cluster info
    sorting_index = sorting_cell{datas,1};
    %get the number of clusters
    group_num = sorting_cell{datas,2};
    % calculate the distances
    distances = distance_calculation(data(datas),sorting_index);
    %assemble a matrix with the cum distributions per cluster

    %allocate memory for the matrix
    cdf_mat = zeros(group_num,n_bins);
    %for all the clusters
    for group = 1:group_num
        % concatenate the distances for all fish
        fish_dist = vertcat(distances{group,:});
        %get the histogram counts
        cdf_mat(group,:) = histcounts(fish_dist,n_bins,'Normalization','cdf');
    end
    %store the info for this file
    files_dist{datas} = cdf_mat;
    
end
%% Generate surrogate distributions for each cluster

% get the number of datasets
num_data = size(data,2);
%allocate memory for the results
random_dist = cell(num_data,1);

%for all the files
for datas = 1:num_data

    %and the cluster info
    sorting_index = sorting_cell{datas,1};
    group_num = sorting_cell{datas,2};
    % allocate a variable to accumulate the cdf
    cdf_cell = cell(randomization,1);
    % for all the randomizations
    for randomized = 1:randomization 
        %store the info for this file
        temp_dist = distance_calculation(data(datas),sorting_index,1);

        %allocate memory for the matrix
        cdf_mat = zeros(group_num,n_bins);
        %for all the clusters
        for group = 1:group_num
            % concatenate the distances for all fish
            fish_dist = vertcat(temp_dist{group,:});
            %get the histogram counts
            cdf_mat(group,:) = histcounts(fish_dist,n_bins,'Normalization','cdf');
        end
        % save the cdf
        cdf_cell{randomized} = cdf_mat;
    end
    % store the average and standard error in the output
    random_dist{datas} = {mean(cat(3,cdf_cell{:}),3),std(cat(3,cdf_cell{:}),0,3)};
end
%% Plot the distributions compare to the surrogate
close all
%for all the files
for datas = 1:num_data
    figure
    % get the cluster number
    group_num = sorting_cell{datas,2};
    members = sorting_cell{datas,3};
    
    % for all the clusters
    for group = 1:group_num
        subplot(round(sqrt(group_num)),ceil(sqrt(group_num)),group)
        plot(1:n_bins,files_dist{datas}(group,:))
        hold on
        shadedErrorBar(1:n_bins,random_dist{datas}{1}(group,:),random_dist{datas}{2}(group,:))
        set(gca,'XTick',[],'YTick',[])
        title(strjoin({num2str(group),num2str(members(group))},'_'),'Interpreter','None')
    end
    t = sgtitle(data(datas).name);
    t.Interpreter = 'None';
    % assemble the figure path
    file_path = strjoin({'distanceCDF',data(datas).name,'sorting',num2str(sorting_constant),...
        'nbins',num2str(n_bins),'rand',num2str(randomization),'.png'},'_');
    saveas(gcf, fullfile(fig_path,file_path), 'png')
end
autoArrangeFigures
%% Plot cluster anatomy
close all

% define the target clusters
target_groups = [1,4];

% for both datasets
for datas = 1:num_data
    % get the idx of the groups
    sorting_index = sorting_cell{datas,1};
    % get the number of groups
    group_num = sorting_cell{datas,2};
    % get the fish information
    fish_ori = data.fish_ori;
    % get the coordinates
    xy_mat = data.xy_seed;
    z_mat = data.z_seed;
    % get the number of fish
    fish_num = size(unique(fish_ori(:,1)),1);

    % for all the fish
    for fish = 1%:fish_num
        figure
        % get the z sections for this animal
        fish_z = z_mat(fish_ori(:,1)==fish);
        % get the number of z sections
        z_num = length(unique(fish_z));
        % get the slices of the stack for this fish
        fish_stack = data(datas).ave_stack(:,:,unique(fish_z));
        % for all the target groups
        for group = 1:length(target_groups)
            % get the coordinates of the seeds in this cluster and fish
            fish_xy = xy_mat(fish_ori(:,1)==fish&sorting_index==target_groups(group));
            % for all the z
            for z = 1:z_num
                % get the current frame
                curr_frame = uint8(fish_stack(:,:,z));
                % run through the seeds
                for seeds = 1:size(fish_xy,1)
                    % plot the seed in the frame
                     curr_frame(fish_xy(seeds).pxlist) = 255;
                end
                subplot(round(sqrt(z_num)),ceil(sqrt(z_num)),z)
                imshow(curr_frame)
            end
        end
    end
    
    
    
end


% Model the interaction between regions

clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'Model\');

data = load_clusters(cluster_path);
%% Run the modelling 
% TODO: use the clusters per region
close all
% define the regions to use
tectum_regions = {'L-TcN','R-TcN','L-TcP','R-TcP','L-Cb','R-Cb','L-Hb','R-Hb','L-Pt','R-Pt'};
tectum_numbers = 1:10;
af_regions = {'AF4','AF5','AF8','AF9','AF10'};
af_numbers = [1 2 5 6 7];

num_datasets = size(data,2);

% allocate memory to store the model data
model_cell = cell(size(data,2),1);
% allocate memory to store the region information
region_cell = cell(num_datasets,2);
% for all the data sets
for datas = 1:num_datasets
    % get the region info
    region_info = data(datas).anatomy_info(:,1);
    % define the regions to be considered (depending on the stimulus protocol)
    if contains(data(datas).name,{'syn','Syn'})
        region_list = af_regions;
        region_numbers = af_numbers;
    else
        region_list = tectum_regions;
        region_numbers = tectum_numbers;
    end
    % get the region numbers 
    num_data = length(region_numbers);

    %get all the pairwise combinations of the regions
    region_comb = [nchoosek(region_numbers,2);fliplr(nchoosek(region_numbers,2))];

    %get the number of combs
    num_comb = size(region_comb,1);
    
    % allocate memory for the model data within this dataset
    current_models = cell(num_comb,1);

    % define the period of interest (0 pre, 1 stim, 2 post, 3 pre-post)
    period = 1;
    % get the target period labeled with ones
    rest_all = period_of_interest(period,data(datas).stim_num,1);
    %for all the combinations
    for combs = 1:num_comb
        %concatenate the clusters from both animals involved
        tar1 = data(datas).region_clusters(region_comb(combs,1)).clu_ave;
        tar2 = data(datas).region_clusters(region_comb(combs,2)).clu_ave;
        % if either of them is empty, put a nan in the cell and skip
        if isempty(tar1) || isempty(tar2)
            current_models{combs} = NaN;
            continue
        end
        
        %for all the averages in 1, calculate models from the raw traces in 2
        %allocate memory to store the model results
        model_para = cell(size(tar1,1),1);
        %for all the averages
        for clu = 1:size(tar1,1)
            fprintf(strcat('Current comb:',num2str(combs),'Current clu: ',num2str(clu),'\r\n'))

            model_para{clu} = fitrlinear(tar2',tar1(clu,:)','CrossVal','on');

        end
        % fill up the dataset cell
        current_models{combs} = model_para;

    end
    % fill up the overall cell
    model_cell{datas} = current_models;
end
%% Plot the results

close all

% for all the data sets
for datas = 1:size(data,2)
    
    % define the regions to be considered (depending on the stimulus protocol)
    if contains(data(datas).name,{'syn','Syn'})
        region_list = af_regions;
        region_numbers = af_numbers;
    else
        region_list = tectum_regions;
        region_numbers = tectum_numbers;
    end
    % get the region numbers 
    num_data = length(region_numbers);
    %get all the pairwise combinations of the regions
    region_comb = [nchoosek(1:num_data,2);fliplr(nchoosek(1:num_data,2))];

    %get the number of combs
    num_comb = size(region_comb,1);
    % allocate memory for the combination matrix
    combination_matrix = zeros(num_data);
    % for all the combinations
    for combs = 1:num_comb
        % get the coordinates of the combination
        target = region_comb(combs,1);
        source = region_comb(combs,2);
        % get the number of clusters in the x_coord
        clu_num = data(datas).region_clusters(region_numbers(target)).clu_num;
        % if the model of interest is a nan, put a nan and skip the
        % iteration
        if ~iscell(model_cell{datas}{combs}) && isnan(model_cell{datas}{combs})
            combination_matrix(target,source) = nan;
            continue
        end
        % calculate the average of the cluster losses for this combination
        % for all the clusters
        for clu = 1:clu_num
            combination_matrix(source,target) = ...
                combination_matrix(source,target) + ...
                kfoldLoss(model_cell{datas}{combs}{clu})/clu_num;
        end
    end
    % plot the matrix
    imagesc(1-combination_matrix)
    if contains(data(datas).name,{'syn','Syn'})
        set(gca,'XTick',1:num_data,'XTickLabel',region_list,'XTickLabelRotation',45)
    else
        set(gca,'XTick',1:2:num_data,'XTickLabel',{'TcN','TcP','Pt','Hb','Cb'},'XTickLabelRotation',45)
    end
    set(gca,'YTick',1:num_data,'YTickLabel',region_list)
    set(gca,'CLim',[0.65 1])
    
    axis square

    % create the settings
    fig_set = struct([]);
    
    fig_set(1).fig_path = fig_path;
    fig_set(1).fig_name = strjoin({'modelMatrix',data(datas).name,'.eps'},'_');
    fig_set(1).fig_size = 5;
    fig_set(1).colorbar = 1;
    fig_set(1).colorbar_label = 'Goodness of Fit';
    fig_set(1).box = 'on';
    
    h = style_figure(gcf,fig_set);
    
end
autoArrangeFigures
%% Simulate seed extraction with seed mixing

%% Load the seeds from the downsampled data and the extracted data
% load the data
clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'mixROI\');

data = load_clusters(cluster_path);
% define the color scheme depending on the stimulus type
if contains(data(1).name,'p17')
    color_scheme = [1 0 0;0 1 0;0 0 1;1 0 1];
    cone_color_scheme = [0.5 0 0;0 0.5 0;0 0 0.5;0.5 0 0.5];
    stim_labels = {'Red','Green','Blue','UV'};
else
%     color_scheme = distinguishable_colors(6);
    color_scheme = [1 0 0;1 0 1;1 0 0;1 0 1;1 0 0;1 0 1];
    cone_color_scheme = [];
    stim_labels = [];
end
% get the number of data sets
num_data = size(data,2);
%% Calculate correlation matrix and sort to look for correspondence

% close all
% 
% figure
% [rho,pval] = corr(data(1).conc_trace',data(2).conc_trace');
% 
% imagesc(rho)
%% Allocate to clusters with GM model and compare

downsampled_traces = data(2).conc_trace;
roi_number = size(downsampled_traces,1);
%define the sPCA parameters to use
bounds_top = 1:data(1).time_num:size(downsampled_traces,2);
bounds_bottom = [bounds_top(2:end)-1,size(downsampled_traces,2)];
bounds = [bounds_top;bounds_bottom];
K = ones(1,data(1).stim_num).*4;
t_bins = ones(data(1).stim_num,1).*10;
pca_vec = ones(data(1).stim_num,1).*1;

%define the vector of cluster numbers to try
clu_vec = [];

replicates = 20;

[~,~,~,~,~,f_data] = sPCA_GMM_cluster_Color(downsampled_traces,bounds...
    ,K,t_bins,pca_vec,[],clu_vec,replicates,[],2);

% use the GMM to allocate the downsampled traces
convolved_clusters = cluster(data(1).GMModel,f_data);
%% Infer the valid clusters by comparing original and cleaned up idx

[~,~,~,~,~,f_data] = sPCA_GMM_cluster_Color(data(1).conc_trace,bounds...
    ,K,t_bins,pca_vec,[],clu_vec,replicates,[],2);

% get the original idx
original_idx = cluster(data(1).GMModel,f_data);
% get the final idx
final_idx = data(1).idx_clu;
% get a vector with the cluster number
original_clunum = data(1).GMModel.NumComponents;

% allocate memory for the LUT
idx_LUT = zeros(original_clunum,2);
% for all the original clusters
for clu = 1:original_clunum
    idx_LUT(clu,1) = clu;
    idx_LUT(clu,2) = mode(final_idx(original_idx==clu));
end
%% Correct the idx for the kernels

for roi = 1:roi_number
    convolved_clusters(roi) = idx_LUT(convolved_clusters(roi)==idx_LUT(:,1),2);
end
%% Compare the clusters allocation between downsampled and the actual data

close all

% assemble a comparison matrix

% allocate memory for the matrix
cluster_matrix = zeros(data(1).clu_num,data(2).clu_num);

% for all the traces
for roi = 1:roi_number
%     % if the Zhou cluster number is not there, skip
%     if any(kernels_valid,kernels_idx(roi)) == 0
%         continue
%     end
    % get the row coordinate (local clusters)
    x = convolved_clusters(roi);
    % get the Zhou cluster number
    y = data(2).idx_clu(roi);
    % add the position of the cluster as a counter
    cluster_matrix(x,y) = cluster_matrix(x,y) + 1;
end

% plot the matrix
imagesc(cluster_matrix)
cmap = magma(256);
cmap(1,:) = [1 1 1];
colormap(cmap)
xlabel('Downsampled clusters')
ylabel('Main clusters')
%% Calculate types and compare

if contains(data(1).name,'p17b')
    
    close all
    % allocate memory to store the matrices
    type_cell = cell(num_data+1,3);

    % for all of the datasets
    for datas = 1:num_data
        
        % use the gains
        delta_norm = data(datas).delta_norm;
        % get the 10th percentile
        zero_threshold = prctile(abs(delta_norm),10,1);
        % zero the values below a threshold
        delta_norm(abs(delta_norm)<zero_threshold&abs(delta_norm)>0) = 0;
        % turn negatives into -1 and positives into 1
        delta_norm(delta_norm>0) = 1;
        delta_norm(delta_norm<0) = -1;
        
        % quantify the occurrence of each pattern
        [pattern,ia,ic] = unique(delta_norm,'rows');
        
        % get the number of patterns
        pattern_num = length(ia);
        
        % allocate vector for the number
        pattern_counts = zeros(pattern_num,1);
        % count the occurrences
        % for all the patterns
        for pat = 1:pattern_num
            pattern_counts(pat) = sum(ic==pat);
        end
        
        % sort by abundance
        [pattern_counts,sort_idx] = sort(pattern_counts,'descend');
        
        pattern = pattern(sort_idx,:);

        % allocate memory for the colors
        pattern_full = zeros(size(pattern,1),4,3);
        % transform the indexes into colors
        for channel = 1:3
            pattern_full(pattern(:,channel)==1,channel,channel) = 1;
            pattern_full(pattern(:,channel)==0,channel,:) = 1;
            if channel == 1
                pattern_full(pattern(:,4)==1,4,[1 3]) = 1;
                pattern_full(pattern(:,4)==0,4,:) = 1;
            end
 
        end
        
        % store the matrix
        type_cell{datas,1} = pattern;
        type_cell{datas,2} = pattern_counts./sum(pattern_counts);
        type_cell{datas,3} = pattern_full;
        
        % eliminate the patterns with only 1 instance
        elim_vector = pattern_counts<2;
        pattern_counts = pattern_counts(~elim_vector);
        pattern_full = pattern_full(~elim_vector,:,:);
        
        figure
        set(gcf,'Color','w')
        subplot(2,1,2)
        image(permute(pattern_full,[2 1 3]))
%         hold on
%         line([1 1],[0 0],'Color','w')
%         set(gca,'XLim',[-0.2 size(pattern_full,1)])
        set(gca,'YScale','linear','XTick',[],'Visible','off')

        subplot(2,1,1)
        bar((pattern_counts))
%         BarPlotBreak(pattern_counts,pattern_counts(2)*1.8,pattern_counts(1)*0.9,'Line',0.6,2)
%         breakplot(1:length(pattern_counts),pattern_counts,400,1900,'Line')
        set(gca,'YScale','linear','XTick',[],'Visible','off')
%         break_axis = breakyaxis([400 1900],0.05, 0.1);
        

        axis tight
   
        % create the settings
        fig_set = struct([]);
        
        fig_set(1).fig_path = fig_path;
        fig_set(1).fig_name = strjoin({'responseTypes',data(datas).name,'.eps'},'_');
        fig_set(1).fig_size = 3.6;
        fig_set(2).fig_size = 3.6;
        fig_set(1).painters = 1;
        fig_set(3).fig_size = 3.6;
        fig_set(4).fig_size = 3.6;
        fig_set(5).fig_size = 3.6;
        fig_set(6).fig_size = 3.6;
        
%         h = style_figure(gcf,fig_set);
        
    end
end

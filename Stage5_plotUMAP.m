clearvars
close all

gains_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\Stage4_gain\';
% cluster_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\Stage3_cluster\';

% data = load_clusters(gains_path,{'_gains'});
data = load_clusters(gains_path);
%% Extract the PCA features from the data

% get the PC info from the structure
pcs = data.pcs(:,1);
% allocate memory for the vector
pc_matrix = zeros(size(data.conc_trace,1),size(pcs{1},2),data.stim_num);
% define the stimulus time
stim_time = 21:60;
% get the data and reshape to separate stim and time
stim_data = reshape(data.conc_trace,[],data.time_num,data.stim_num);
% for all the stimuli
for stim = 1:data.stim_num
    pc_matrix(:,:,stim) = stim_data(:,stim_time,stim)*pcs{stim};
end
% concatenate and column normalize the matrix for use
pc_matrix = normr_1(reshape(pc_matrix,[],size(pcs{1},2)*data.stim_num),2);
%% Use UMAP to embed the gains


[reduced_data, umap] = run_umap(cat(2,pc_matrix,normr_1(data.delta_norm,2)), 'n_neighbors', 10, 'min_dist', 0.1);
%% Plot the results
close all
num_points = size(reduced_data,1);
% close all

for i = 1:4
    color_raw = round(normr_1(data.delta_norm(:,i),1)*200);
    color_edges = [min(color_raw),max(color_raw)];
    cmap = jet(diff(color_edges)+1);

    figure
%     for points = 1:num_points
%         plot(reduced_data(points,1),reduced_data(points,2),'o',...
%             'MarkerFaceColor',cmap(color_raw(points)-color_edges(1)+1,:),...
%             'MarkerEdgeColor',cmap(color_raw(points)-color_edges(1)+1,:))
%         hold on
%     end

    scatter(reduced_data(:,1),reduced_data(:,2),[],cmap(color_raw-color_edges(1)+1,:),'filled')
end

figure
% color_raw = round(normr_1(data.fish_ori(:,1),1)*200);
% color_raw = round(normr_1(data.idx_clu,1)*200);
color_raw = data.idx_clu;
color_edges = [min(color_raw),max(color_raw)];
cmap = jet(diff(color_edges)+1);
scatter(reduced_data(:,1),reduced_data(:,2),[],cmap(color_raw-color_edges(1)+1,:))

autoArrangeFigures
%%
% close all
% 
% plot(delta_norm(:,2),delta_norm(:,4),'o')

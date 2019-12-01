clearvars
close all

gains_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\Stage3_cluster\';
% cluster_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\Stage3_cluster\';

% data = load_clusters(gains_path,{'_gains'});
data = load_clusters(gains_path);
%define the stim labels based on the paradigm
if contains(data(1).name,'syn')
    %define the region labels
    reg_label = {'AF4','AF5','AF6','AF7','AF8','AF9','AF10'};
    reg_map = [0 4 5 6 7 8 9 10];
else
    %define the region labels
    reg_label = {'L-TcN','R-TcN','L-TcP','R-TcP','L-Cb','R-Cb','L-Hb','R-Hb','L-Pt','R-Pt'};
    reg_map = [0:10];
end
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


% [reduced_data, umap] = run_umap(cat(2,pc_matrix,normr_1(data.delta_norm,2)), 'n_neighbors', 10, 'min_dist', 0.1);

color_raw = data.anatomy_info(:,1);
color_raw(isnan(color_raw)) = 0;
[reduced_data, umap] = run_umap(cat(2,pc_matrix,normr_1(data.delta_norm,2),normr_1(color_raw,2)), 'n_neighbors', 10, 'min_dist', 0.1);
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
% color_raw = data.idx_clu;
color_raw = data.anatomy_info(:,1);
color_raw(isnan(color_raw)) = 0;
% color_raw = round(normr_1(pc_matrix(:,4),1)*200);
color_edges = [nanmin(color_raw),nanmax(color_raw)];
cmap = hsv(diff(color_edges)+1);
% s = scatter(reduced_data(:,1),reduced_data(:,2),[],cmap(color_raw-color_edges(1)+1,:));
s = gscatter(reduced_data(:,1),reduced_data(:,2),color_raw,cmap,'.',10);
legend(reg_label)

autoArrangeFigures
%%
% close all
% 
% plot(delta_norm(:,2),delta_norm(:,4),'o')

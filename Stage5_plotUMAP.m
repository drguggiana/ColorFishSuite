clearvars
close all

load('paths.mat')
addpath(genpath(paths(1).main_path))

main_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'UMAP\');

% reset the rng
rng(1);

% data = load_clusters(gains_path,{'_gains'});
data = load_clusters(main_path);

% get the number of dataset
num_data = size(data,2);
% allocate memory for the index
index_cell = cell(num_data,1);
% allocate memory for the region list
region_list = cell(num_data,1);
% define the region set to use (1 tc vs rgc, 2 all vs all)
region_set = 1;

% for all the datasets
for datas = 1:num_data
    region_combination = 1;
    
    % define which regions to keep depending on the dataset
    if contains(data(datas).name, {'Syn','syn'})
        if region_set == 1
            region_list{datas} = {'AF10'};
        elseif region_set == 2
            region_list{datas} = {'AF4','AF5','AF6','AF7','AF8','AF9','AF10'};
        end
    else
        if region_set == 1
            region_list{datas} = {'R-TcN','R-TcP'};
        elseif region_set == 2
            region_list{datas} = {'L-TcN','R-TcN','L-TcP','R-TcP','L-Cb','R-Cb','L-Hb','R-Hb','L-Pt','R-Pt'};
        end
    end
    % load the anatomy info
    anatomy_info = data(datas).anatomy_info(:,1);
    
    % separate the traces by region
    [region_cell,~] = region_split(data(datas).single_reps,...
        anatomy_info,data(datas).name,region_combination,region_list{datas});
    % rewrite the index vector
    index_cell{datas} = region_cell{3}==1;
    
end

% Allocate memory to store the embeddings
UMAP_cell = cell(num_data,1);
% for all the datasets selected
for datas = 1:length(data)
    close all
    %scan for the p17b
    if contains(data(datas).name,'p17b')
        %if it's p17b
        stim_labels = {'Red','Green','Blue','UV'};
        %define the plot colors
        plot_col = [1 0 0;0 1 0;0 0 1;1 0 1];
        % set gain to on
        gains = 1;
        % get the gains
        delta_norm = data(datas).delta_norm(index_cell{datas}==1,:);
    else %if it's p6p8 instead
        %define the stim labels (for p6p8 data)
        stim_labels = {'Red CK','UV CK','Red GR','UV GR','Red FL','UV FL'};
        %define the plot colors
%         plot_col = [1 0 0;0 0 1;0.8 0 0;0 0 0.8;0.4 0 0;0 0 0.4];
        plot_col = [1 0 0;1 0 1;1 0 0;1 0 1;1 0 0;1 0 1];
        gains = 0;
    end
    %% Extract the PCA features from the data
    
    % get the PC info from the structure
    pcs = data(datas).pcs(:,1);
    % get the data
    conc_trace = data(datas).conc_trace(index_cell{datas}==1,:);

    % allocate memory for the vector
    pc_matrix = zeros(size(conc_trace,1),size(pcs{1},2),data(datas).stim_num);
    % define the stimulus time
    stim_time = 21:60;
    % get the data(datas) and reshape to separate stim and time
    stim_data = reshape(conc_trace,[],data(datas).time_num,data(datas).stim_num);
    % for all the stimuli
    for stim = 1:data(datas).stim_num
        pc_matrix(:,:,stim) = stim_data(:,stim_time,stim)*pcs{stim};
    end
    % concatenate and column normalize the matrix for use
    pc_matrix = normr_1(reshape(pc_matrix,[],size(pcs{1},2)*data(datas).stim_num),2);
    %% Use UMAP to embed the gains
    
    [reduced_data, umap] = run_umap(pc_matrix, 'n_neighbors', 10, 'min_dist', 0.1);
%     [reduced_data, umap] = run_umap(cat(2,pc_matrix,normr_1(data(datas).delta_norm,2)), 'n_neighbors', 10, 'min_dist', 0.1);
    % store the embedding
    UMAP_cell{datas} = reduced_data;
end
%% Plot the results based on stimulus
close all

for datas = 1:num_data
    
    % load the embedding
    reduced_data = UMAP_cell{datas};
    num_points = size(reduced_data,1);
    % close all
    % get the number of traces, time and stimuli
    conc_trace = data(datas).conc_trace(index_cell{datas}==1,:);
    
    trace_num = size(conc_trace,1);
    time_num = data(datas).time_num;
    stim_num = data(datas).stim_num;
    
    
%     % get the reshaped activity
%     average_levels = reshape(conc_trace,trace_num,time_num,stim_num);
%     % take only the stimulation time
%     average_levels = average_levels(:,stim_time,:);
    average_levels = abs(data(datas).delta_norm(index_cell{datas}==1,:));
    % take the absolute average
%     average_levels = squeeze(mean(abs(average_levels),2));
    if gains
        for i = 1:4
            % scale and center the gain values
            color_raw = round(log(normr_1(delta_norm(:,i),1))*200);
            % remove the inf values
            color_raw(isinf(color_raw)) = 0;
            % get the range of values
            color_edges = [min(color_raw),max(color_raw)];
            % get the intensity values for the colormap
            intensities = linspace(0, 1, diff(color_edges)+1)';
            % get the vector of zeros to fill up the other channels
            filler = zeros(diff(color_edges)+1,1);
            % depending on the channel, define the colormap
            switch i
                case 1
                    cmap = [intensities,filler,filler];
                case 2
                    cmap = [filler,intensities,filler];
                case 3
                    cmap = [filler,filler,intensities];
                case 4
                    cmap = [intensities,filler,intensities];
            end

%             figure
%             scatter(reduced_data(:,1),reduced_data(:,2),[],cmap(color_raw-color_edges(1)+1,:),'filled')
%             title(strjoin({'UMAP',data(datas).name,'Gain',num2str(i)},'_'),'Interpreter','None')
%             set(gca,'TickLength',[0 0])
%             file_path = strjoin({'UMAP',data(datas).name,'Gain',num2str(i),'set',num2str(region_set),'.png'},'_');
%             saveas(gcf, fullfile(fig_path,file_path), 'png')
        end
        figure
%         % use the max value of the parameter (either average resp or gain)
%         [~,max_values] = max(average_levels,[],2);
        % use the top percentile of the gain
        perc_values = prctile(average_levels,75,1);
        % allocate memory for the max_values
        max_values = zeros(size(average_levels,1),1);
        temp_values = zeros(size(average_levels));
        % for all colors
        for stim = 1:stim_num
            % get the ROIs passing threshold in each color
            temp_values(average_levels(:,stim)>perc_values(stim),stim) = 1;
            
            max_values(average_levels(:,stim)>perc_values(stim)) = stim;
        end
        
        % get the ROIs with more than 1 selectivity
        opponents = sum(temp_values,2)+3;
        opponents(opponents<5) = 0;
        % join it with the max_values
        max_values(opponents~=0) = opponents(opponents~=0);
        max_values = max_values + 1;
        
%         cmap = [255,0,0; 0,255,0; 0,0,255; 255,0,255]/255;
        cmap = [0.8 0.8 0.8;1 0 0;0 1 0;0 0 1;1 0 1;0 1 1;0.8 0.8 0;0 0 0];
%         scatter(reduced_data(:,1),reduced_data(:,2),[],cmap(max_values,:),'filled','o')
        gscatter(reduced_data(:,1),reduced_data(:,2),max_values,cmap,'.',10)
%         legend({'Red', 'Green', 'Blue', 'UV'})
%         legend({'NR','Red', 'Green', 'Blue', 'UV','2-C','3-C','4-C'},...
%             'Location','best','Orientation','horizontal')
        legend('off')
        set(gca,'TickLength',[0 0],'visible','off')
        title(strjoin({'UMAP',data(datas).name,'Combined'},'_'),'Interpreter','None')
        file_path = strjoin({'UMAP',data(datas).name,'Combined','set',num2str(region_set),'.png'},'_');
%         saveas(gcf, fullfile(fig_path,file_path), 'png')
        print(fullfile(fig_path,file_path),'-dpng','-r600')
    else

%         for i = 1:data(datas).stim_num
%             average_perstim = average_levels(:,i);
%             color_raw = round(log(normr_1(average_perstim,1))*200);
%             color_raw(isinf(color_raw)) = 0;
%             color_edges = [min(color_raw),max(color_raw)];
%             % get the intensity values for the colormap
%             intensities = linspace(0, 1, diff(color_edges)+1)';
%             % get the vector of zeros to fill up the other channels
%             filler = zeros(diff(color_edges)+1,1);
%             % depending on the channel, define the colormap
%             switch i
%                 case 1
%                     cmap = [intensities,filler,filler];
%                     marker = 'o';
%                 case 2
%                     cmap = [intensities,filler,intensities];
%                     marker = 'o';
%                 case 3
%                     cmap = [intensities,filler,filler];
%                     marker = 's';
%                 case 4
%                     cmap = [intensities,filler,intensities];
%                     marker = 's';
%                 case 5
%                     cmap = [intensities,filler,filler];
%                     marker = 'd';
%                 case 6
%                     cmap = [intensities,filler,intensities];
%                     marker = 'd';
%             end
% %             cmap = jet(diff(color_edges)+1);
% 
% %             figure
% % 
% %             scatter(reduced_data(:,1),reduced_data(:,2),[],cmap(color_raw-color_edges(1)+1,:),'filled',marker)
% %             title(strjoin({data(datas).name,'Stim',stim_labels{i}},'_'),'Interpreter','None')
% %             set(gca,'TickLength',[0 0],'visible','off')
% %             file_path = strjoin({'UMAP',data(datas).name,'Stim',stim_labels{i},'set',num2str(region_set),'.png'},'_');
% %             saveas(gcf, fullfile(fig_path,file_path), 'png')
%         end
        
        figure
        [~,max_values] = max(average_levels,[],2);
        % color order: maroon (R CK), navy (U CK), orange (R GR), blue (U
        % GR), pink (R FL), lavender (U FL)
%         cmap = [128,0,0; 0,0,128; 245,130,48; 0,130,200; 250,190,190; 255,0,255]/255;
        cmap = [1 0 0;1 0 1;1 0 0;1 0 1;0.5 0 0;0.5 0 0.5];
        emap = [1 0 0;1 0 1;0 0 0;0 0 0;0.5 0 0;0.5 0 0.5];
        markers = {'o','o','s','s','d','d'};
        for stim = 1:stim_num
            
            scatter(reduced_data(max_values==stim,1),reduced_data(max_values==stim,2),...
                [],markers{stim},'MarkerEdgeColor',emap(stim,:),'MarkerFaceColor',cmap(stim,:))
            hold on
        end
        legend({'Red CK', 'UV CK', 'Red GR', 'UV GR', 'Red FL', 'UV FL'},'Location','northwest')
        colormap(cmap)
        colorbar('Ticklabels',stim_labels,'TickLabelInterpreter','None','Ticks',linspace(0.08,0.92,6))
        title(strjoin({data(datas).name,'Combined'},'_'),'Interpreter','None')
        set(gca,'TickLength',[0 0],'visible','off')
        file_path = strjoin({'UMAP',data(datas).name,'Combined','set',num2str(region_set),'.png'},'_');
%         saveas(gcf, fullfile(fig_path,file_path), 'png')
        print(fullfile(fig_path,file_path),'-dpng','-r600')
    end
end
%% Plot the results based on region
for datas = 1:num_data

    % load the embedding
    reduced_data = UMAP_cell{datas};
    if contains(data(datas).name,'p17') && region_set == 2
        figure
        % color_raw = round(normr_1(data.fish_ori(:,1),1)*200);
        % color_raw = round(normr_1(data.idx_clu,1)*200);
        % color_raw = data.idx_clu;
        color_raw = data(datas).anatomy_info(index_cell{datas}==1,1);
        color_raw(isnan(color_raw)) = 0;
        % color_raw = round(normr_1(pc_matrix(:,4),1)*200);
        color_edges = [nanmin(color_raw),nanmax(color_raw)];
%         cmap = lines(diff(color_edges)+1); 
        cmap = distinguishable_colors(diff(color_edges)+1);
        if contains(data(datas).name,'syn')
            cmap(end,:) = [0.8 0.8 0.8];
        end
        % s = scatter(reduced_data(:,1),reduced_data(:,2),[],cmap(color_raw-color_edges(1)+1,:));
        s = gscatter(reduced_data(:,1),reduced_data(:,2),color_raw,cmap,'.',10);
        legend(region_list{datas})
        set(gca,'TickLength',[0 0],'visible','off')
        % assemble the figure path
        file_path = strjoin({'UMAP',data(datas).name,'Region','set',num2str(region_set),'.png'},'_');
%         saveas(gcf, fullfile(fig_path,file_path), 'png')
        print(fullfile(fig_path,file_path),'-dpng','-r600')
    end
    
    autoArrangeFigures
end
%% Plot the clusters in the embedding
for datas = 1:num_data
    
    % load the embedding
    reduced_data = UMAP_cell{datas};
%     close all
    figure
    % get the cluster indexes
    idx_clu = data(datas).idx_clu(index_cell{datas}==1);
    
    % get the number of clusters
    clu_num = data(datas).clu_num;
    
    % set a colormap for them
    cluster_cmap = distinguishable_colors(clu_num,[0 0 0;1 1 1]);
    % plot the scatter
    s = gscatter(reduced_data(:,1),reduced_data(:,2),idx_clu,cluster_cmap,'.',10);
    legend(string(1:clu_num))
    set(gca,'TickLength',[0 0],'visible','off')
    % assemble the figure path
    file_path = strjoin({'UMAP',data(datas).name,'Cluster','set',num2str(region_set),'.png'},'_');
%     saveas(gcf, fullfile(fig_path,file_path), 'png')
    print(fullfile(fig_path,file_path),'-dpng','-r600')
    
end
%% Plot the p8 UMAPs separately by stimulus modality
close all
if contains(data(datas).name,'p8')
    for datas = 1:num_data

        % load the embedding
        reduced_data = UMAP_cell{datas};
        num_points = size(reduced_data,1);
        % close all
        % get the number of traces, time and stimuli
        trace_num = size(conc_trace,1);
        time_num = data(datas).time_num;
        stim_num = data(datas).stim_num;

        % get the reshaped activity
        average_levels = reshape(conc_trace,trace_num,time_num,stim_num);
        % take only the stimulation time
        average_levels = average_levels(:,stim_time,:);
        % take the absolute average
        average_levels = squeeze(mean(abs(average_levels),2));

        [~,max_values] = max(average_levels,[],2);

        % define the colormap
        cmap = [0.8 0.8 0.8;1 0 0;1 0 1];
        % for all stimulus types
        for stim = 1:2:5
    %         fig('units','centimeters','width',10,'height',8,'fontsize',5)
            figure

            % create a temporary color scheme to plot the other stim gray
            color_raw = zeros(size(reduced_data,1),1);
            color_raw(max_values==stim) = 2;
            color_raw(max_values==stim+1) = 3;
            color_raw(color_raw==0) = 1;
            % plot the scatter
            s = gscatter(reduced_data(:,1),reduced_data(:,2),color_raw,cmap,'.',...
                20);
    %         scatter(reduced_data(max_values==stim,1),reduced_data(max_values==stim,2),...
    %             [],'o','MarkerFaceColor',cmap(1,:))
    %         hold on
    %         scatter(reduced_data(max_values==stim+1,1),reduced_data(max_values==stim+1,2),...
    %             [],'o','MarkerFaceColor',cmap(2,:))

            legend(s(2:3),stim_labels([stim,stim+1]),'Location','northwest','Fontsize',20)
            colormap(cmap)
    %         colorbar('Ticklabels',stim_labels,'TickLabelInterpreter','None','Ticks',linspace(0.08,0.92,6))
            title(strjoin({data(datas).name,'Combined'},'_'),'Interpreter','None')
            set(gca,'TickLength',[0 0],'visible','off')
            file_path = strjoin({'UMAP',data(datas).name,'Stim',num2str(stim),'set',num2str(region_set),'.png'},'_');
    %         saveas(gcf, fullfile(fig_path,file_path), 'png')
            print(fullfile(fig_path,file_path),'-dpng','-r600')
        end

    end
end
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
hist_colors = zeros(num_data,3);
% for all the datasets
for datas = 1:num_data
    region_combination = 1;
    
    % define which regions to keep depending on the dataset
    if contains(data(datas).name, {'Syn','syn'})
        if region_set == 1
            region_list{datas} = {'AF10'};
            hist_colors(datas,:) = paths.afOT_colors(datas,:);
        elseif region_set == 2
            region_list{datas} = {'AF4','AF5','AF6','AF7','AF8','AF9','AF10'};
        end
    else
        if region_set == 1
            region_list{datas} = {'R-TcN','R-TcP'};
            hist_colors(datas,:) = paths.afOT_colors(datas,:);
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

% define the labels for the gain histogram
gain_labels = {'NR','Red','Green','Blue','UV','2-chrom','3-chrom','4-chrom'};

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
histo = figure;

% define whether to use gains or average levels
metric = 'average';
% allocate memory to store the max values for stats
max_cell = cell(num_data,1);
for datas = num_data:-1:1
    
    % load the embedding
    reduced_data = UMAP_cell{datas};
    num_points = size(reduced_data,1);
    % close all
    % get the number of traces, time and stimuli
    conc_trace = data(datas).conc_trace(index_cell{datas}==1,:);
    
    trace_num = size(conc_trace,1);
    time_num = data(datas).time_num;
    stim_num = data(datas).stim_num;
    
    switch metric
        case 'average'
            % get the reshaped activity
            average_levels = reshape(conc_trace,trace_num,time_num,stim_num);
            % take only the stimulation time
            average_levels = average_levels(:,stim_time,:);
            % take the absolute average
            average_levels = squeeze(mean(abs(average_levels),2));
        case 'gains'
            % use the actual gains
            average_levels = abs(data(datas).delta_norm(index_cell{datas}==1,:));
    end

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
%         % use the max value of the parameter (either average resp or gain)
%         [~,max_values] = max(average_levels,[],2);
        % use the top percentile of the gain
        perc_values = prctile(average_levels,75,1);
        % allocate memory for the max_values
        max_values = zeros(size(average_levels,1),1);
        temp_values = zeros(size(average_levels));
        % for all colors
        for stim = 1:size(average_levels,2)
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
        % store the value for stats
        max_cell{datas} = max_values;
%         histogram(max_values,'Normalization','probability','LineWidth',1,'FaceColor',hist_colors(3-datas,:))
%         hold on

        if datas == 1
            figure(histo)
            
            % allocate memory for the bin counts
            bin_counts = cell(2,1);
            % get the histcounts
            for i = 1:2
                bin_counts{i} = histcounts(max_cell{i},'Normalization','probability');
            end
%             histogram(max_values,'Normalization','probability','LineWidth',1,'FaceColor',hist_colors(3-datas,:))
            bars = bar(vertcat(bin_counts{:})');
            tint_colors = tint_colormap(hist_colors,0);
            set(bars(1),'FaceColor',tint_colors(2,:))
            set(bars(2),'FaceColor',tint_colors(1,:))
            
            box off
            set(gca,'XTick',1:8,'XTickLabels',gain_labels,'XTickLabelRotation',45)
            set(gca,'FontSize',7,'LineWidth',1,'TickLength',[0 0])
%             ylabel('Proportion ROIs')
%             title('Response Gains')  
%             legend({'AF10','Tectum'})
%             set(gcf,'Color','w')
%             file_path = fullfile(fig_path,strjoin({'histoUMAP',metric,data(1).name,data(2).name,...
%                 'Combined','set',num2str(region_set),'.png'},'_'));
%             export_fig(file_path,'-r600')
            
            % create the settings
            fig_set = struct([]);
            
            fig_set(1).fig_path = fig_path;
            fig_set(1).fig_name = strjoin({'histoUMAP',metric,data(1).name,data(2).name,...
                'Combined','set',num2str(region_set),'.eps'},'_');
            fig_set(1).fig_size = [4.2 3.5];
            fig_set(1).painters = 1;
%             fig_set(1).linewidth = 1;
            
            h = style_figure(gcf,fig_set);
        end
        
        figure
        
%         cmap = [255,0,0; 0,255,0; 0,0,255; 255,0,255]/255;
        cmap = [0.8 0.8 0.8;1 0 0;0 1 0;0 0 1;1 0 1;0 0 0;0 0 0;0 0 0];
%         scatter(reduced_data(:,1),reduced_data(:,2),[],cmap(max_values,:),'filled','o')
        gscatter(reduced_data(:,1),reduced_data(:,2),max_values,cmap,'.',10)
%         legend({'Red', 'Green', 'Blue', 'UV'})
%         legend({'NR','Red', 'Green', 'Blue', 'UV','2-C','3-C','4-C'},...
%             'Location','best','Orientation','horizontal')
        legend('off')
        set(gca,'TickLength',[0 0],'visible','off','LineWidth',2)
        set(gcf,'Color','w')
        title(strjoin({'UMAP',data(datas).name,'Combined'},'_'),'Interpreter','None')
        file_path = fullfile(fig_path,strjoin({'UMAP',metric,data(datas).name,...
            'Combined','set',num2str(region_set),'.png'},'_'));
        export_fig(file_path,'-r600')
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
        file_path = fullfile(fig_path,strjoin({'UMAP',data(datas).name,...
            'Combined','set',num2str(region_set),'.png'},'_'));
        export_fig(file_path,'-r600')
    end
%     autoArrangeFigures
end
%% Run stats on the histogram of types

% allocate memory for the counts for both data sets
full_counts = cell(num_data,1);
% for both data sets
for datas = 1:num_data
    % get the fish info
    fish_info = data(datas).fish_ori(index_cell{datas}==1,1);
    % get the number of fish
    num_fish = length(unique(fish_info));
    
    % get the types of bins and their number
    bin_types = unique(max_cell{datas});
    type_number = length(bin_types);
    
    % allocate memory for the count
    count_per_fish = zeros(num_fish,type_number);
    
    % for the types of bins
    for bintype = 1:type_number
        % for all the fish
        for fish = 1:num_fish
            % get the current max values
            current_values = max_cell{datas}(fish_info==fish);
            count_per_fish(fish,bintype) = sum(current_values==bin_types(bintype));
        end
    end
    % store the results
    full_counts{datas} = count_per_fish;
end

% allocate memory for the test results
test_result = zeros(type_number,1);
% run wilcoxon's sign rank for each type
for bintype = 1:type_number
    test_result(bintype) = signrank(full_counts{1}(:,bintype),full_counts{2}(:,bintype));
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
        [leg,obj_list] = legendflex(region_list{datas},'ncol',1,...
            'box','off','anchor',[4 8],'padding',[0 0 0],'buffer',[0 0]);
        
        % run through the markers to make them bigger
        temp_obj = findobj(obj_list,'type','Line');
        temp_text = findobj(obj_list,'type','Text');
        % for all the markers
        for marks = 2:2:length(temp_obj)
            set(temp_obj(marks),'MarkerSize',20)
            set(temp_text(marks/2),'FontSize',12)
        end
        set(gca,'TickLength',[0 0],'visible','off','LineWidth',2)
        % adapt the figure size to encompass the legend
        set(gca,'Position',get(gca,'Position')-[0 0 0.15 0.15])
        
        set(gca,'TickLength',[0 0],'visible','off','LineWidth',2)


        % create the settings
        fig_set = struct([]);
        
        fig_set(1).fig_path = fig_path;
        fig_set(1).fig_name = strjoin({'UMAP',data(datas).name,...
            'Region','set',num2str(region_set),'.png'},'_');
        fig_set(1).fig_size = 3.6;
        fig_set(2).skip = 1;
        
        h = style_figure(gcf,fig_set);
    end
    
%     autoArrangeFigures
end
%% Plot the clusters in the embedding
close all
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
    [leg,obj_list] = legendflex(cellstr(string(1:clu_num)),'ncol',2,...
        'box','off','anchor',[4 8],'padding',[0 0 0],'buffer',[0 0]);

    % run through the markers to make them bigger
    temp_obj = findobj(obj_list,'type','Line');
    temp_text = findobj(obj_list,'type','Text');
    % for all the markers
    for marks = 2:2:length(temp_obj)
        set(temp_obj(marks),'MarkerSize',20)
        set(temp_text(marks/2),'FontSize',12)
    end
    set(gca,'TickLength',[0 0],'visible','off','LineWidth',2)
    % adapt the figure size to encompass the legend
    set(gca,'Position',get(gca,'Position')-[0 0 0.15 0.15])

    
    % create the settings
    fig_set = struct([]);
    
    fig_set(1).fig_path = fig_path;
    fig_set(1).fig_name = strjoin({'UMAP',data(datas).name,...
        'Cluster','set',num2str(region_set),'.png'},'_');
    fig_set(1).fig_size = 3.6;
%     % invert the figure list to put the legend last
%     fig_list = get(gcf,'Children');
%     set(gcf,'Children',flipud(fig_list));
    % skip editing the flex legend
    fig_set(2).skip = 1;
    
    h = style_figure(gcf,fig_set);
    
end
%% Plot the p8 UMAPs separately by stimulus modality
close all
if contains(data(datas).name,'p8')
    % define the labels of the patterns
    pattern_label = {'Checker','Gratings','Flash'};
    % define the dataset names
    fig_name = {'Tectum','AF10'};
    % allocate memory to store the color raw vectors
    color_cell = cell(num_data,stim_num/2);
    % for all datasets
    for datas = 1:num_data

        % load the embedding
        reduced_data = UMAP_cell{datas};
        num_points = size(reduced_data,1);
        % close all
        % get the number of traces, time and stimuli
        conc_trace = data(datas).conc_trace;
        trace_num = size(conc_trace,1);
        time_num = data(datas).time_num;
        stim_num = data(datas).stim_num;

        % get the reshaped activity
        average_levels2 = reshape(conc_trace,trace_num,time_num,stim_num);
        % take only the stimulation time
        average_levels2 = average_levels2(:,stim_time,:);
        % take the absolute average
        average_levels2 = squeeze(max(abs(average_levels2),[],2));

        % define the colormap
        cmap = [0.8 0.8 0.8;1 0 0;1 0 1;0 0 0];
        % for all stimulus types
        for stim = 1:2:5
            figure
            
            % copy the corresponding average_levels
            average_levels = average_levels2(:,stim:stim+1);
            
            perc_values = prctile(average_levels,75,1);
            % allocate memory for the max_values
            max_values = zeros(size(average_levels,1),1);
            temp_values = zeros(size(average_levels));
            % for all colors
            for stim2 = 1:2
                % get the ROIs passing threshold in each color
                temp_values(average_levels(:,stim2)>perc_values(stim2),stim2) = 1;
                
                max_values(average_levels(:,stim2)>perc_values(stim2)) = stim2;
            end
            
            % get the ROIs with more than 1 selectivity
            opponents = sum(temp_values,2)+1;
            opponents(opponents<3) = 0;
            % join it with the max_values
            max_values(opponents~=0) = opponents(opponents~=0);
            color_raw = max_values + 1;

            % save the color_raw vector for quantification
            color_cell{datas,(stim+1)/2} = color_raw;
            % plot the scatter
            s = gscatter(reduced_data(:,1),reduced_data(:,2),color_raw,cmap,'.',...
                20,'off');

            colormap(cmap)
            set(gcf,'Color','w')
            title(strjoin({data(datas).name,'Combined'},'_'),'Interpreter','None')
            set(gca,'TickLength',[0 0],'visible','off')
%             file_path = fullfile(fig_path,strjoin({'UMAP',data(datas).name,...
%                 'Stim',num2str(stim),'set',num2str(region_set),'.png'},'_'));
%             export_fig(file_path,'-r600')
            
            % create the settings
            fig_set = struct([]);
            
            fig_set(1).fig_path = fig_path;
            fig_set(1).fig_name = strjoin({'UMAP',data(datas).name,...
                'Stim',num2str(stim),'set',num2str(region_set),'.png'},'_');
            fig_set(1).fig_size = 2.46;
            
            h = style_figure(gcf,fig_set);
        end
        %% Produce a ternary plot
        close all
        figure
        % collapse the vectors for this dataset
        color_all = horzcat(color_cell{datas,:}); 
        
        % calculate a red-UV index
        % allocate memory for it
        red_uv_index = zeros(size(average_levels2,1),3);
        for pair = 1:2:6
            stim1 = average_levels2(:,pair);
            stim2 = average_levels2(:,pair+1);
%             red_uv_index(:,(pair+1)/2) = (stim1-stim2)./(stim1+stim2);
            red_uv_index(:,(pair+1)/2) = stim1+stim2;
        end
        
        % get the sizes
        total_size = sum(average_levels2,2);
        total_size = total_size.*50;
        
        % normalize the index across rows
%         red_uv_index = log(normr_1(red_uv_index,1));
        red_uv_index = normr_1(red_uv_index,1);
%         % get the unique patterns
%         [unique_seq,ia,ic] = unique(color_all,'rows');
%         % allocate memory for the amounts of ROIs
%         roi_counts = zeros(size(unique_seq,1),1);
%         % get the counts per pattern
%         for pattern = 1:size(unique_seq)
%             roi_counts(pattern) = sum(ic==pattern);
%         end
%         
%         
%         roi_counts = log(roi_counts);
%         roi_counts(roi_counts==0) = 1;
%         roi_counts = roi_counts*20;
        
%         imagesc(sortrows(color_all))
%         ternplot(unique_seq(:,1),unique_seq(:,2),unique_seq(:,3),roi_counts,'o','filled','majors',4)
        h = ternplot(red_uv_index(:,1),red_uv_index(:,2),red_uv_index(:,3),200,...
            cmap(color_raw,:),'.','majors',4);
        children = get(gca,'Children');
        % for all but the first (i.e. all the text objects)
        for i = 2:length(children)
            set(children(i),'FontSize',15)
            if i > 5 && i < 10
                % get the current position
                current_position = get(children(i),'Position');
                % displace the text
                current_position(1) = current_position(1) - 0.09;
                set(children(i),'Position',current_position)
            elseif i > 9
                % get the current position
                current_position = get(children(i),'Position');
                % displace the text
                current_position(2) = current_position(2) - 0.03;
                set(children(i),'Position',current_position)
            end
        end
        
%         h = ternlabel(pattern_label{:});
%         set(gca,'TickLength',[0 0],'XTick',1:3,'XTickLabels',{'Checker','Grating','Flash'})
%         set(gca,'YDir','normal','FontSize',15)
%         ylabel('ROIs')
%         colormap(cmap)
        title(fig_name{datas})
        set(gca,'FontSize',20)
        set(gcf,'Color','w')
        file_path = fullfile(fig_path,strjoin({'barsUMAP',data(datas).name,'.png'},'_'));
        
        export_fig(file_path,'-r600')
        %% Produce a graph diagram
        close all
        figure
        
        fontsize = 15;
        % turn the color matrix into binary, excluding mixed selectivity
        % cells
        binary_matrix = color_all(~any(color_all>3,2),:);
        binary_matrix = binary_matrix>1;
        % get the number of multicolor cells
        multi_number = size(color_all,1)-size(binary_matrix,1);
        
        % get the unique patterns
        [unique_seq,ia,ic] = unique(binary_matrix,'rows');
        % allocate memory for the amounts of ROIs
        roi_counts = zeros(size(unique_seq,1),1);
        % get the counts per pattern
        for pattern = 1:size(unique_seq)
            roi_counts(pattern) = sum(ic==pattern);
        end
        % reorder the pattern
        roi_counts = roi_counts([1 2 3 5 4 6 7 8]);
        unique_seq = unique_seq([1 2 3 5 4 6 7 8],:);
        % draw the venn diagram
%         [H,S] = venn(roi_counts(2:end));
%         schemaball


        plot([-1 1],[-1 -1],'-k','LineWidth',1)
        hold on
        plot([1 0],[-1 1],'-k','LineWidth',1)
        plot([-1 0],[-1 1],'-k','LineWidth',1)
        
        plot([-1 0],[-1 0],'-k','LineWidth',1)
        plot([1 0],[-1 0],'-k','LineWidth',1)
        plot([0 0],[1 0],'-k','LineWidth',1)

        % define the vectors of coordinates and the colors
        x_points = [-1 1 0 0 -0.5 0.5 0];
        y_points = [-1 -1 1 -1 0 0 0];
        colors = [2 3 4 5 6 7 8];
        scatter(x_points,y_points,roi_counts(colors).*10,...
            (roi_counts(colors)),'o','filled','MarkerEdgeColor',[0 0 0])
        scatter(-0.8,0.8,roi_counts(1).*10,[0.8 0.8 0.8],'filled')
        scatter(0.8,0.8,multi_number.*10,[0 0 0],'filled')
        text(-0.8,0.8,{num2str(roi_counts(1))},'Color',[0 0 0],...
            'FontWeight','bold','HorizontalAlignment','center','FontSize',fontsize)
        text(0.8,0.8,{num2str(multi_number)},'Color',[1 1 1],...
            'FontWeight','bold','HorizontalAlignment','center','FontSize',fontsize)
        cba = colorbar;
        set(cba,'TickLength',0,'LineWidth',2,'FontSize',fontsize)
        ylabel(cba,'Nr ROIs','FontSize',fontsize)
        
        set(gca,'TickLength',[0 0],'XTick',[],'YTick',[],'visible','off')
        axis equal
        colormap(magma)
%         set(gcf,'Color','w')
%         file_path = fullfile(fig_path,strjoin({'vennUMAP',data(datas).name,'.png'},'_'));
%         print(file_path,'-dpng','-r600','-painters')
        
        
        % create the settings
        fig_set = struct([]);
        
        fig_set(1).fig_path = fig_path;
        fig_set(1).fig_name = strjoin({'graphUMAP',data(datas).name,'.png'},'_');
        fig_set(1).fig_size = 2.91;
        fig_set(1).painters = 1;
        
        h = style_figure(gcf,fig_set);
        
%         export_fig(file_path,'-r600','-nocrop')
%         roi_counts = log(roi_counts);
%         roi_counts(roi_counts==0) = 1;
%     end
    %% Produce a venn diagram
            close all
%     for datas = 1:num_data
        figure
        
        fontsize = 15;
        % turn the color matrix into binary, excluding mixed selectivity
        % cells
        binary_matrix = color_all(~any(color_all>3,2),:);
        binary_matrix = binary_matrix>1;
        % get the number of multicolor cells
        multi_number = size(color_all,1)-size(binary_matrix,1);
        
        % get the unique patterns
        [unique_seq,ia,ic] = unique(binary_matrix,'rows');
        % allocate memory for the amounts of ROIs
        roi_counts = zeros(size(unique_seq,1),1);
        % get the counts per pattern
        for pattern = 1:size(unique_seq)
            roi_counts(pattern) = sum(ic==pattern);
        end
        % reorder the pattern
        roi_counts = roi_counts([1 2 3 5 4 6 7 8]);
        unique_seq = unique_seq([1 2 3 5 4 6 7 8],:);
        % draw the venn diagram
        
        colors = magma(3);
        [H,S] = venn(roi_counts(2:end),'FaceColor',{colors(1,:), colors(2,:), colors(3,:)});
        hold on
        
        
        % get the normalized numbers as percentages
        roi_percentages = round(roi_counts*100./sum(roi_counts));
        
        % place text on the circles
        % for all the areas
        for areas = 1:7
            text(S.ZoneCentroid(areas,1),S.ZoneCentroid(areas,2),...
                strcat(num2str(roi_percentages(areas+1)),'%'),...
                'FontSize',7,'HorizontalAlignment','center')
        end
        
        set(gca,'Visible','off')
        axis equal
        fig_set = struct([]);
        
        fig_set(1).fig_path = fig_path;
        fig_set(1).fig_name = strjoin({'vennUMAP',data(datas).name,'.png'},'_');
        fig_set(1).fig_size = 3.8;
        fig_set(1).painters = 1;
%         fig_set(1).visible = 'off';
        
        h = style_figure(gcf,fig_set);
    end

    
end
%% Calculate distance distributions in PCA space

if contains(data(1).name,'p17b')
    close all
    figure
    
    % get the number and list of types
    type_list = unique(max_cell{1});
    type_number = length(type_list);
    % allocate memory for the medians and mads
    median_mad = zeros(type_number,num_data,2);
    % for both data sets
    for datas = 1:num_data
       % get the individual fish information
       fish_ori = data(datas).fish_ori(:,1);
       % get the number of fish
       num_fish = length(unique(fish_ori));
       
       % get the traces and get only the stimulation period
       conc_trace = reshape(data(datas).conc_trace,[],time_num,stim_num);
       conc_trace  = reshape(conc_trace(:,21:60,:),[],40*stim_num);

       % allocate memory for the distances
       distance_cell = cell(num_fish,type_number);
       % for all the fish
       for fish = 1:num_fish
           % get the traces from that animal
           fish_traces = conc_trace(fish_ori==fish&index_cell{datas},:);
           % get the percentile assignment for that animal
           fish_partial = fish_ori(index_cell{datas});
           perc_values = max_cell{datas}(fish_partial==fish);

           % for all the types
           for types = 1:type_number
               % perform PCA on the traces
               [~,score] = pca(fish_traces(perc_values==type_list(types),:));
               % if there's not enough dimensionality, skip
               if size(score,2) < 5
                   continue
               end
               % calculate the distance between points in 3d PCA space
               distance_cell{fish,types} = pdist(score(:,1:5));
           end

       end
       % concatenate the results across fish and plot

       % for all the types
       for types = 1:type_number
            distribution = horzcat(distance_cell{:,types});
            subplot(round(sqrt(type_number)),ceil(sqrt(type_number)),types)
            histogram(distribution,50,'Normalization','pdf');
%             N = histcounts(distribution,50,'Normalization','cdf');
%             plot(N)
            hold on
            
            median_mad(types,datas,1) = median(distribution);
            median_mad(types,datas,2) = mad(distribution);
            
       end
    end
    
    % plot the median and mad
    figure
    bar(median_mad(:,:,1))
    for datas = 1:num_data
%         bar(median_mad(:,datas,1),'FaceAlpha',0.5)
        hold on
        errorbar(1:type_number,median_mad(:,datas,1),median_mad(:,datas,2),'.k')
    end
end
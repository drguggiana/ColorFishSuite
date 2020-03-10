clearvars
close all

load('paths.mat')
addpath(genpath(paths(1).main_path))

main_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'UMAP\');

% data = load_clusters(gains_path,{'_gains'});
data = load_clusters(main_path);

% for all the datasets selected
for datas = 1:length(data)
    close all
    
    %define the stim labels based on the paradigm
    if contains(data(datas).name,'syn')
        %define the region labels
        reg_label = {'AF4','AF5','AF6','AF7','AF8','AF9','AF10'};
        reg_map = [0 4 5 6 7 8 9 10];
    else
        %define the region labels
        reg_label = {'L-TcN','R-TcN','L-TcP','R-TcP','L-Cb','R-Cb','L-Hb','R-Hb','L-Pt','R-Pt'};
        reg_map = [0:10];
    end
    %scan for the p17b
    if contains(data(datas).name,'p17b')
        %if it's p17b
        stim_labels = {'Red','Green','Blue','UV'};
        %define the plot colors
        plot_col = [1 0 0;0 1 0;0 0 1;1 0 1];
        % set gain to on
        gains = 1;
    else %if it's p6p8 instead
        %define the stim labels (for p6p8 data)
        stim_labels = {'Red CK','UV CK','Red GR','UV GR','Red FL','UV FL'};
        %define the plot colors
        plot_col = [1 0 0;0 0 1;0.8 0 0;0 0 0.8;0.4 0 0;0 0 0.4];
        gains = 0;
    end
    %% Extract the PCA features from the data
    
    % get the PC info from the structure
    pcs = data(datas).pcs(:,1);
    % allocate memory for the vector
    pc_matrix = zeros(size(data(datas).conc_trace,1),size(pcs{1},2),data(datas).stim_num);
    % define the stimulus time
    stim_time = 21:60;
    % get the data(datas) and reshape to separate stim and time
    stim_data = reshape(data(datas).conc_trace,[],data(datas).time_num,data(datas).stim_num);
    % for all the stimuli
    for stim = 1:data(datas).stim_num
        pc_matrix(:,:,stim) = stim_data(:,stim_time,stim)*pcs{stim};
    end
    % concatenate and column normalize the matrix for use
    pc_matrix = normr_1(reshape(pc_matrix,[],size(pcs{1},2)*data(datas).stim_num),2);
    %% Use UMAP to embed the gains
    
    [reduced_data, umap] = run_umap(pc_matrix, 'n_neighbors', 10, 'min_dist', 0.1);
%     [reduced_data, umap] = run_umap(cat(2,pc_matrix,normr_1(data(datas).delta_norm,2)), 'n_neighbors', 10, 'min_dist', 0.1);
    
%     color_raw = data.anatomy_info(:,1);
%     color_raw(isnan(color_raw)) = 0;
%     [reduced_data, umap] = run_umap(cat(2,pc_matrix,log(normr_1(data.delta_norm,2)),normr_1(color_raw,2)), 'n_neighbors', 10, 'min_dist', 0.1);
    %% Plot the results
    close all
    num_points = size(reduced_data,1);
    % close all
    % get the number of traces, time and stimuli
    trace_num = size(data(datas).conc_trace,1);
    time_num = data(datas).time_num;
    stim_num = data(datas).stim_num;
    
    % get the reshaped activity
    average_levels = reshape(data(datas).conc_trace,trace_num,time_num,stim_num);
    % take only the stimulation time
    average_levels = average_levels(:,stim_time,:);
    % take the absolute average
    average_levels = squeeze(mean(abs(average_levels),2));
    if gains
        for i = 1:4
            % scale and center the gain values
            color_raw = round(log(normr_1(data(datas).delta_norm(:,i),1))*200);
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

            figure
            scatter(reduced_data(:,1),reduced_data(:,2),[],cmap(color_raw-color_edges(1)+1,:),'filled')
            title(strjoin({'UMAP',data(datas).name,'Gain',num2str(i)},'_'),'Interpreter','None');
            file_path = strjoin({'UMAP',data(datas).name,'Gain',num2str(i),'.png'},'_');
            saveas(gcf, fullfile(fig_path,file_path), 'png')
        end
        figure
        [~,max_values] = max(average_levels,[],2);
        cmap = [255,0,0; 0,255,0; 0,0,255; 255,0,255]/255;
        scatter(reduced_data(:,1),reduced_data(:,2),[],cmap(max_values,:),'filled','o')
        title(strjoin({'UMAP',data(datas).name,'Combined'},'_'),'Interpreter','None');
        file_path = strjoin({'UMAP',data(datas).name,'Combined','.png'},'_');
        saveas(gcf, fullfile(fig_path,file_path), 'png')
    else

        for i = 1:data(datas).stim_num
            average_perstim = average_levels(:,i);
            color_raw = round(log(normr_1(average_perstim,1))*200);
            color_raw(isinf(color_raw)) = 0;
            color_edges = [min(color_raw),max(color_raw)];
            % get the intensity values for the colormap
            intensities = linspace(0, 1, diff(color_edges)+1)';
            % get the vector of zeros to fill up the other channels
            filler = zeros(diff(color_edges)+1,1);
            % depending on the channel, define the colormap
            switch i
                case 1
                    cmap = [intensities,filler,filler];
                    marker = 'o';
                case 2
                    cmap = [intensities,filler,intensities];
                    marker = 'o';
                case 3
                    cmap = [intensities,filler,filler];
                    marker = 's';
                case 4
                    cmap = [intensities,filler,intensities];
                    marker = 's';
                case 5
                    cmap = [intensities,filler,filler];
                    marker = 'd';
                case 6
                    cmap = [intensities,filler,intensities];
                    marker = 'd';
            end
%             cmap = jet(diff(color_edges)+1);

            figure

            scatter(reduced_data(:,1),reduced_data(:,2),[],cmap(color_raw-color_edges(1)+1,:),'filled',marker)
            title(strjoin({data(datas).name,'Stim',stim_labels{i}},'_'),'Interpreter','None')
            file_path = strjoin({'UMAP',data(datas).name,'Stim',stim_labels{i},'.png'},'_');
            saveas(gcf, fullfile(fig_path,file_path), 'png')
        end
        
        figure
        [~,max_values] = max(average_levels,[],2);
        % color order: maroon (R CK), navy (U CK), orange (R GR), blue (U
        % GR), pink (R FL), lavender (U FL)
        cmap = [128,0,0; 0,0,128; 245,130,48; 0,130,200; 250,190,190; 255,0,255]/255;
        scatter(reduced_data(:,1),reduced_data(:,2),[],cmap(max_values,:),'filled','o')
%         legend({'Red CK', 'UV CK', 'Red GR', 'UB GR', 'Red FL', 'UV FL'})
        colormap(cmap)
        colorbar('Ticklabels',stim_labels,'TickLabelInterpreter','None','Ticks',linspace(0.08,0.92,6))
        title(strjoin({data(datas).name,'Combined'},'_'),'Interpreter','None')
        file_path = strjoin({'UMAP',data(datas).name,'Combined','.png'},'_');
        saveas(gcf, fullfile(fig_path,file_path), 'png')
    end
    
    
%     figure
%     % color_raw = round(normr_1(data.fish_ori(:,1),1)*200);
%     % color_raw = round(normr_1(data.idx_clu,1)*200);
%     % color_raw = data.idx_clu;
%     color_raw = data(datas).anatomy_info(:,1);
%     color_raw(isnan(color_raw)) = 0;
%     % color_raw = round(normr_1(pc_matrix(:,4),1)*200);
%     color_edges = [nanmin(color_raw),nanmax(color_raw)];
%     cmap = hsv(diff(color_edges)+1);
%     % s = scatter(reduced_data(:,1),reduced_data(:,2),[],cmap(color_raw-color_edges(1)+1,:));
%     s = gscatter(reduced_data(:,1),reduced_data(:,2),color_raw,cmap,'.',10);
%     legend(reg_label)
%     
%     % assemble the figure path
%     file_path = strjoin({'UMAP',data(datas).name,'Region','.png'},'_');
%     saveas(gcf, fullfile(fig_path,file_path), 'png')
    
    autoArrangeFigures
end
%%
% close all
% 
% plot(delta_norm(:,2),delta_norm(:,4),'o')

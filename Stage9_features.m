% Calculate calcium response standard features per color

% load the data
clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'Features\');

data = load_clusters(cluster_path);
% define the color scheme depending on the stimulus type
if contains(data(1).name,'p17')
    color_scheme = [1 0 0;0 1 0;0 0 1;1 0 1];
    cone_color_scheme = [0.5 0 0;0 0.5 0;0 0 0.5;0.5 0 0.5];
else
    color_scheme = distinguishable_colors(6);
    cone_color_scheme = [];
end
% get the number of data sets
num_data = size(data,2);
%% OFF Define a subset of areas to work with


% % get the number of dataset
% num_data = size(data,2);
% % allocate memory for the index
% index_cell = cell(num_data,1);
% % for all the datasets
% for datas = 1:num_data
%     region_combination = 1;
%     
%     % define which regions to keep depending on the dataset
%     if contains(data(datas).name, {'Syn','syn'})
%         region_list = {'AF10'};
%     else
%         region_list = {'R-TcN','R-TcP'};
%     end
%     % load the anatomy info
%     anatomy_info = data(datas).anatomy_info(:,1);
%     
%     % separate the traces by region
%     [region_cell,~] = region_split(data(datas).single_reps,...
%         anatomy_info,data(datas).name,region_combination,region_list);
%     % rewrite the index vector
%     index_cell{datas} = region_cell{3}==1;
%     
% end
%% Calculate overall average features

close all

% allocate memory for the averages
ave_cell = cell(num_data,2);
figure
% for all the data sets
for datas = 1:num_data
    % calculate the average and std of the response
    ave_resp = nanmean(data(datas).conc_trace,1);
    std_resp = nanstd(data(datas).conc_trace,0,1);
    % store in the cell
    ave_cell{datas} = cat(1,ave_resp,std_resp);
    % plot
    shadedErrorBar(1:size(ave_resp,2),ave_resp + datas*5,std_resp)
    hold on
    
end
title(strjoin({'Average Trace',data(:).figure_name},' '),'Interpreter','None')
file_path = strjoin({'AverageTrace',data(:).name,'.png'},'_');
saveas(gcf, fullfile(fig_path,file_path), 'png')
%% Calculate averages per area

close all
% define the trace offset
trace_offset = 5;
% define the framerate
framerate = 1/0.952;
if contains(data(1).name,'p17b')

    % allocate memory for the output
    average_perarea = cell(num_data,1);

    % for all of the datasets
    for datas = 1:num_data
        % allocate memory for the 
        % get the number of areas
        area_list = unique(data(datas).anatomy_info(:,1));
        area_list = area_list(~isnan(area_list));
        area_number = size(area_list,1);

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

        figure
        % initialize an area counter
        a_count = 1;
        % initialize a vector to store the areas plotted
        plotted_area = zeros(area_number,1);
        % for all the areas
        for area = 1:area_number
            % if it's AF7, skip
            if contains(reg_label{area},'AF7')
                continue
            end
            % get the average trace
            ave_trace = nanmean(data(datas).conc_trace(data(datas).anatomy_info(:,1)==area_list(area),:),1);
            std_trace = nanstd(data(datas).conc_trace(data(datas).anatomy_info(:,1)==area_list(area),:),0,1);
            ave_perstim = reshape(ave_trace,[],data(datas).stim_num);
            std_perstim = reshape(std_trace,[],data(datas).stim_num);
            time_vector = (0:size(ave_trace,2)-1)./framerate;
            time_perstim = reshape(time_vector,[],data(datas).stim_num);
            % split by stimulus
            for stim = 1:data(datas).stim_num
                % plot it
                shadedErrorBar(time_perstim(:,stim),...
                    ave_perstim(:,stim)+(a_count-1)*trace_offset,std_perstim(:,stim),...
                    {'color',color_scheme(stim,:),'LineWidth',1})
                hold on
            end

            % update the counter
            plotted_area(area) = 1;
            a_count = a_count + 1;

        end
        
        % plot the intermediate lines
        for stim = 1:data(datas).stim_num-1
            plot([time_perstim(end,stim),time_perstim(end,stim)],...
                get(gca,'YLim'),'k','LineWidth',2)
        end
        pbaspect([2,1,1])
        axis tight
        set(gca,'YTick',0:trace_offset:(a_count-2)*trace_offset,'YTickLabels',reg_label(plotted_area==1))
        set(gca,'TickLength',[0 0])
        xlabel('Time (s)')
        sgtitle(strcat(data(datas).figure_name,'+'),'Interpreter','None')
%         set(gca,'XLim',[0,time_perstim(end,end)])
        box off

        
        % save the figure
        file_path = strjoin({'Average_Trace_perArea',data(datas).name,'.png'},'_');
        saveas(gcf, fullfile(fig_path,file_path), 'png')
    end
    autoArrangeFigures
end
%% Calculate max response, delay to peak

close all
% define the number of parameters
param_num = 3;
param_label = {'LogMax response','Delay to peak','Abs Logmean response'};
% define the interval to look at
target_interval = 5:35;

% allocate memory to store the numbers per trace
calcium_cell = cell(num_data,1);
% for all of the datasets
for datas = 1:num_data
    % get the traces
    conc_trace = data(datas).conc_trace;
    % get the number of stimuli
    stim_num = data(datas).stim_num;
    
    % get the number of traces
    num_traces = size(conc_trace,1);
    % reshape the traces per stimulus
    conc_trace = reshape(conc_trace,num_traces,[],stim_num);
    % allocate memory for the calculations
    calcium_matrix = zeros(num_traces,stim_num,param_num);
    % for all the traces
    for trace = 1:num_traces
        % for all the stimuli
        for stim = 1:stim_num
            % get the max response and the delay to peak
            [calcium_matrix(trace,stim,1),calcium_matrix(trace,stim,2)] = nanmax(conc_trace(trace,target_interval,stim));
            % log the max response
            calcium_matrix(trace,stim,1) = log(calcium_matrix(trace,stim,1));
            calcium_matrix(trace,stim,3) = nanmean(abs(log(conc_trace(trace,target_interval,stim))),2);
        end
    end
    % store the results for later use
    calcium_cell{datas} = calcium_matrix;
    % plot the results
    for param = 1:param_num
        figure
        [h,L,MX,MED] = violin(squeeze(calcium_matrix(:,:,param)),'mc',[],'medc','k');
        set(L,'visible','off')
        % for all the stimuli
        for stim = 1:data(datas).stim_num
            set(h(stim),'facecolor',color_scheme(stim,:))
        end
        hold on
        set(gca,'XTick',[],'TickLength',[0 0])
        axis tight
        box off
    %         plotSpread(squeeze(calcium_matrix(:,:,param)))
        title(data(datas).figure_name)
        ylabel(param_label{param},'Interpreter','None')
        file_path = strjoin({param_label{param},data(datas).name,'.png'},'_');
        saveas(gcf, fullfile(fig_path,file_path), 'png')
    end
    
end

autoArrangeFigures
%% Plot the above results per area
close all
if contains(data(1).name,'p17b')

    % for all of the datasets
    for datas = 1:num_data

        % get the corresponding calcium data
        calcium_matrix = calcium_cell{datas};
        % get the number of areas
        area_list = unique(data(datas).anatomy_info(:,1));
        area_list = area_list(~isnan(area_list));
        area_number = size(area_list,1);

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

        % for all the parameters
        for param = 1:param_num
            figure
            % for all the areas
            for area = 1:area_number
                % get the average trace
                plot_matrix = calcium_matrix(data(datas).anatomy_info(:,1)==area_list(area),:,param);
                % plot it
                subplot(round(sqrt(area_number)),ceil(sqrt(area_number)),area)
    %             violin(plot_matrix);
                plotSpread(plot_matrix,'distributionColors',color_scheme);
                ylabel(reg_label{area})
                set(gca,'XTick',[])
            end
            sgtitle(strjoin({param_label{param},'perRegion',data(datas).figure_name},' '),'Interpreter','None')
            file_path = strjoin({param_label{param},'perRegion',data(datas).name,'.png'},'_');
            saveas(gcf, fullfile(fig_path,file_path), 'png')
        end
    end
    autoArrangeFigures
end
%% Plot the parameters per cluster

close all
if contains(data(1).name,'p17b')

    % for all of the datasets
    for datas = 1:num_data

        % get the number of areas
        area_list = unique(data(datas).anatomy_info(:,1));
        area_list = area_list(~isnan(area_list));
        area_number = size(area_list,1);

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

        % allocate a cell to store the data
        clu_cell = cell(area_number,1);
        % for all the areas
        for area = 1:area_number
            % get the cluster averages
            clu_ave = data(datas).region_clusters(area).clu_ave;
            % if there are not clusters, skip
            if isempty(clu_ave)
                continue
            end
            % reshape the matrix accordingly
            clu_ave = reshape(clu_ave,size(clu_ave,1),[],stim_num);
            % allocate memory for the results
            clu_param = zeros(size(clu_ave,1),stim_num,param_num);
            % for all of the stimuli
            for stim = 1:stim_num
                % calculate the parameters of interest and store
                [clu_param(:,stim,1),clu_param(:,stim,2)] = nanmax(clu_ave(:,target_interval,stim),[],2);
                clu_param(:,stim,3) = nanmean(abs(clu_ave(:,target_interval,stim)),2);
            end
            % store in the cell
            clu_cell{area} = clu_param;
        end
        % plot the results
        % for all the parameters
        for param = 1:param_num
            figure
            % for all the areas
            for area = 1:area_number

                if isempty(clu_cell{area})
                    continue
                end
                % get the data
                plot_matrix = clu_cell{area}(:,:,param);

                % plot it
                subplot(round(sqrt(area_number)),ceil(sqrt(area_number)),area)
    %             violin(plot_matrix);
                plotSpread(plot_matrix);
                ylabel(reg_label{area})
            end
            sgtitle(strjoin({param_label{param},'perRegionCluster',data(datas).figure_name},' '),'Interpreter','None')
            file_path = strjoin({param_label{param},'perRegionCluster',data(datas).name,'.png'},'_');
            saveas(gcf, fullfile(fig_path,file_path), 'png')
        end
    end
    autoArrangeFigures
end
%% Average Gain plots
close all
if contains(data(1).name,'p17b')
    % for all of the datasets
    for datas = 1:num_data

        % get the corresponding calcium data
        plot_matrix = data(datas).delta_norm;

        figure
  
        % get the average trace
%         plot_matrix = calcium_matrix(data(datas).anatomy_info(:,1)==area_list(area),:);
        % plot it
%         subplot(round(sqrt(area_number)),ceil(sqrt(area_number)),area)
%         violin(plot_matrix);
        plotSpread(plot_matrix,'distributionColors',cone_color_scheme);
%         ylabel(reg_label{area})
        set(gca,'XTick',[],'TickLength',[0 0])
        axis tight
        sgtitle(data(datas).figure_name,'Interpreter','None')
        ylabel('Gain (a.u.)')
        file_path = strjoin({'Gain',data(datas).name,'.png'},'_');
        saveas(gcf, fullfile(fig_path,file_path), 'png')
    end
    autoArrangeFigures
end
%% Gain plots per region

close all
if contains(data(1).name,'p17b')
    % for all of the datasets
    for datas = 1:num_data

        % get the corresponding calcium data
        calcium_matrix = data(datas).delta_norm;
        % get the number of areas
        area_list = unique(data(datas).anatomy_info(:,1));
        area_list = area_list(~isnan(area_list));
        area_number = size(area_list,1);

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


        figure
        % for all the areas
        for area = 1:area_number
            % get the average trace
            plot_matrix = calcium_matrix(data(datas).anatomy_info(:,1)==area_list(area),:);
            % plot it
            subplot(round(sqrt(area_number)),ceil(sqrt(area_number)),area)
            %             violin(plot_matrix);
            plotSpread(plot_matrix);
            ylabel(reg_label{area})
        end
        sgtitle(strjoin({'GainPerRegion',data(datas).figure_name},' '),'Interpreter','None')
        file_path = strjoin({'GainPerRegion',data(datas).name,'.png'},'_');
        saveas(gcf, fullfile(fig_path,file_path), 'png')
    end
    autoArrangeFigures
end
%% Plot the gains per cluster
if contains(data(1).name,'p17b')
    close all

    % for all of the datasets
    for datas = 1:num_data

        % get the number of areas
        area_list = unique(data(datas).anatomy_info(:,1));
        area_list = area_list(~isnan(area_list));
        area_number = size(area_list,1);

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

        % allocate a cell to store the data
        clu_cell = cell(area_number,1);
        % for all the areas
        for area = 1:area_number
            % get the cluster averages
            clu_idx = data(datas).region_clusters(area).idx_clu;
            clu_num = data(datas).region_clusters(area).clu_num;
            area_idx = data(datas).anatomy_info(:,1)==area_list(area);
            stim_num = data(datas).stim_num;
            gains = data(datas).delta_norm;
            % calculate the average gains
            % allocate memory for the gains
            gain_ave = zeros(clu_num,stim_num);
            % for all the clusters
            for clu = 1:clu_num
                temp = gains(area_idx,:);
                gain_ave(clu,:) = nanmean(temp(clu_idx==clu,:));
            end

            % if there are not clusters, skip
            if isempty(clu_ave)
                continue
            end
    %         % reshape the matrix accordingly
    %         clu_ave = reshape(clu_ave,size(clu_ave,1),[],stim_num);
    %         % allocate memory for the results
    %         clu_param = zeros(size(clu_ave,1),stim_num,param_num);
    %         % for all of the stimuli
    %         for stim = 1:stim_num
    %             % calculate the parameters of interest and store
    %             [clu_param(:,stim,1),clu_param(:,stim,2)] = nanmax(clu_ave(:,target_interval,stim),[],2);
    %             clu_param(:,stim,3) = nanmean(abs(clu_ave(:,target_interval,stim)),2);
    %         end
            % store in the cell
            clu_cell{area} = gain_ave;
        end
        % plot the results

        figure
        % for all the areas
        for area = 1:area_number

            if isempty(clu_cell{area})
                continue
            end
            % get the data
            plot_matrix = clu_cell{area};

            % plot it
            subplot(round(sqrt(area_number)),ceil(sqrt(area_number)),area)
            %             violin(plot_matrix);
            plotSpread(plot_matrix);
            ylabel(reg_label{area})
        end
        sgtitle(strjoin({'GainPerRegionCluster',data(datas).figure_name},' '),'Interpreter','None')
        file_path = strjoin({'GainPerRegionCluster',data(datas).name,'.png'},'_');
        saveas(gcf, fullfile(fig_path,file_path), 'png')
    end
    autoArrangeFigures
end
%% Plot the region clusters in table format
if contains(data(1).name,'p17b')

    close all
    
    % for all of the datasets
    for datas = 1:num_data
        % get the cluster indexes
        idx_clu = data(datas).idx_clu;
        % get the gains
        delta_norm = data(datas).delta_norm;
        % get the cluster number
        clu_num = data(datas).clu_num;
        % get the traces per cluster
        trace_number = data(datas).clu_number;
        % allocate memory for the cluster gains
        cluster_gains = zeros(clu_num,4);
        % for all the clusters
        for clu = 1:clu_num
            %calculate the average gain
            cluster_gains(clu,:) = mean(delta_norm(idx_clu==clu,:),1);
        end
        % sort the gains by trace number
        [sort_trace,sort_idx] = sort(trace_number);
        cluster_gains = cluster_gains(sort_idx,:);
        % plot the gains
        fig('height',14,'width',7)
        [X,Y] = meshgrid(1:4,1:clu_num);
        % for all colors
        for color = 1:4
            scatter(X(:,color),Y(:,color),normr_1(abs(cluster_gains(:,color)),1).*150+1,...
                cone_color_scheme(color,:),'o','filled')
            hold on
            % select the negative gains
            negative_idx = cluster_gains(:,color)<0;
            % draw a black outline outside the negative circles
            scatter(X(negative_idx,color),Y(negative_idx,color),...
                normr_1(abs(cluster_gains(negative_idx,color)),1).*150+1,'w','*')
        end
        set(gca,'YTick',1:clu_num,'YTickLabel',sort_trace,'TickLength',[0 0])
        set(gca,'XLim',[0.8,4.2],'YLim',[0,clu_num+1],'XTick',[])
        pbaspect([1,2,1])
        title(data(datas).figure_name,'Interpreter','None')
        
        file_path = strjoin({'AverageGainsCluster',data(datas).name,'.png'},'_');
        saveas(gcf, fullfile(fig_path,file_path), 'png')
    end
    autoArrangeFigures
end
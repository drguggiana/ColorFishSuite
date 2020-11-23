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
    stim_labels = {'Red','Green','Blue','UV'};
else
%     color_scheme = distinguishable_colors(6);
    color_scheme = [1 0 0;1 0 1;1 0 0;1 0 1;1 0 0;1 0 1];
    cone_color_scheme = [];
    stim_labels = [];
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

% define the framerate
framerate = 1/0.952;
if contains(data(1).name,'p17b')

    % allocate memory for the output
    average_perarea = cell(num_data,1);

    % for all of the datasets
    for datas = 1:num_data
        % get the anatomy info
        anatomy_info = data(datas).anatomy_info(:,1);


        %define the stim labels based on the paradigm
        if contains(data(datas).name,'syn')
            %define the region labels
            reg_label = {'AF4','AF5','AF6','AF7','AF8','AF9','AF10'};
            reg_map = [0 4 5 6 7 8 9 10];
            % define the trace offset
            trace_offset = 5;
        else
            %define the region labels
%             reg_label = {'L-TcN','R-TcN','L-TcP','R-TcP','L-Cb','R-Cb','L-Hb','R-Hb','L-Pt','R-Pt'};
            reg_label = {'TcN','TcP','Cb','Hb','Pt'};
            anatomy_info(anatomy_info==2) = 1;
            anatomy_info(anatomy_info==4) = 3;
            anatomy_info(anatomy_info==6) = 5;
            anatomy_info(anatomy_info==8) = 7;
            anatomy_info(anatomy_info==10) = 9;
            % define the trace offset
            trace_offset = 3;
%             reg_map = [0:10];
        end
        
        % get the number of areas
        area_list = unique(anatomy_info);
        area_list = area_list(~isnan(area_list));
        area_number = size(area_list,1);

        figure
        % initialize an area counter
        a_count = 1;
        % initialize a vector to store the areas plotted
        plotted_area = zeros(area_number,1);
        % for all the areas
        for area = 1:area_number
            % if it's AF7, skip
            if contains(reg_label{area},{'AF7','AF6'})
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
                [-1 (a_count-1)*trace_offset],'k','LineWidth',2)
        end
%         pbaspect([2,1,1])
        axis tight
        set(gca,'YTick',0:trace_offset:(a_count-2)*trace_offset,'YTickLabels',reg_label(plotted_area==1))
        set(gca,'TickLength',[0 0],'LineWidth',2)
        if datas == 2
            set(gca,'XTick',[])
        else
            xlabel('Time (s)')
        end
%         title(data(datas).figure_name,'Interpreter','None')
        ylabel(data(datas).figure_name)

%         set(gca,'YLim',[-0.5,])
%         box off
        axis tight
        
%         set(gca,'FontSize',20)
%         % save the figure
%         file_path = strjoin({'Average_Trace_perArea',data(datas).name,'.png'},'_');
% %         saveas(gcf, fullfile(fig_path,file_path), 'png')
%         print(fullfile(fig_path,file_path),'-dpng','-r600')
        % create the settings
        fig_set = struct([]);
        
        fig_set(1).fig_path = fig_path;
        fig_set(1).fig_name = strjoin({'Average_Trace_perArea',data(datas).name,'.png'},'_');
        fig_set(1).fig_size = 3.6;
        
        h = style_figure(gcf,fig_set);
    end
    autoArrangeFigures
end
%% Calculate max response, delay to peak

close all

% define the number of parameters
param_num = 3;
param_label = {'LogMax response','Delay to peak','Abs Logmean response'};
param_plot_label = {'LogMax response (delta F/F)','Delay to max (s)','Abs Logmean response (delta F/F)'};
% define the plot axis limits
y_lim = [-4 4;-1.7 13.5;0 5];
% define the interval to look at
% target_interval = 5:35;
target_interval = 20:30;

% define the panel label
if contains(data(1).name,'p17b')
    fig_name = {data.figure_name};
else
    fig_name = {'Tectum','AF10'};
end

% allocate memory to store the numbers per trace
calcium_cell = cell(num_data,1);
% also allocate memory for the kw test
kw_cell = cell(num_data,param_num);
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
            [calcium_matrix(trace,stim,1),calcium_matrix(trace,stim,2)] = ...
                nanmax(abs(conc_trace(trace,target_interval,stim)));
            % get the time units
            calcium_matrix(trace,stim,2) = calcium_matrix(trace,stim,2)./data(datas).framerate;
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
        [h,L,MX,MED] = violin(squeeze(calcium_matrix(:,:,param)),'mc',[],'medc','k','facealpha',1);
        set(L,'visible','off')
        % for all the stimuli
        for stim = 1:data(datas).stim_num
            set(h(stim),'facecolor',color_scheme(stim,:))
        end
        hold on
        set(gca,'XTick',[],'TickLength',[0 0],'LineWidth',3)
        
        axis tight
        box off
    %         plotSpread(squeeze(calcium_matrix(:,:,param)))
%         title(fig_name{datas})
%         ylabel(param_plot_label{param},'Interpreter','None')
%         set(gca,'FontSize',20)
        set(gca,'YLim',y_lim(param,:))
        
        % create the settings
        fig_set = struct([]);
        
        fig_set(1).fig_path = fig_path;
        fig_set(1).fig_name = strjoin({param_label{param},data(datas).name,'.png'},'_');
        if contains(data(datas).name,'p17b')
            fig_set(1).fig_size = 1.5;
        else
            fig_set(1).fig_size = 1.8;
        end
        h = style_figure(gcf,fig_set);
%         set(gcf,'Color','w')
%         file_path = fullfile(fig_path,strjoin({param_label{param},data(datas).name,'.png'},'_'));
%         export_fig(file_path,'-r600')
        
%         [kw_cell{datas,param},tbl,stats] = kruskalwallis(squeeze(calcium_matrix(:,:,param)),[],'off');
%         [kw_cell{datas,param},tbl,stats] = friedman(squeeze(calcium_matrix(:,:,param)),1,'off');
%         s = multcompare(stats,'CType','bonferroni','Display','on');


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
% define the fontsize
fontsize = 20;
if contains(data(1).name,'p17b')
    % for all of the datasets
    for datas = 1:num_data

        % get the corresponding calcium data
        plot_matrix = data(datas).delta_norm;
        
%         plot_matrix(:,4) = plot_matrix(:,4)./10;

%         fig('units','centimeters','height',6,'width',9,'fontsize',20)
        figure
        
%         violin(plot_matrix);
        plotSpread(plot_matrix,'distributionColors',cone_color_scheme.*2);
%         boxplot(plot_matrix)
        set(gca,'XTick',[],'TickLength',[0 0])
        axis tight
%         title(data(datas).figure_name,'FontName','Arial')
        ylabel('Gain (a.u.)','FontName','Arial')
%         set(gca,'FontSize',fontsize,'LineWidth',2)
%         file_path = fullfile(fig_path,strjoin({'Gain',data(datas).name,'.png'},'_'));
% %         print(fullfile(fig_path,file_path),'-dpng','-r600')
%         set(gcf,'Color','w')
%         export_fig(file_path,'-r600')
        
        
        % create the settings
        fig_set = struct([]);
        
        fig_set(1).fig_path = fig_path;
        fig_set(1).fig_name = strjoin({'Gain',data(datas).name,'.eps'},'_');
        fig_set(1).fig_size = 4;
        
        h = style_figure(gcf,fig_set);
        
        
%         [~,tbl,stats] = friedman(squeeze(plot_matrix),1,'off');
        [~,tbl,stats] = anova2(squeeze(plot_matrix),1,'off');
        s = multcompare(stats,'CType','bonferroni','Display','off');
    end
%     autoArrangeFigures
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
%         fig('height',14,'width',7)
        figure
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
                normr_1(abs(cluster_gains(negative_idx,color)),1).*150+10,'w','.')
        end
        set(gca,'YTick',1:clu_num,'YTickLabel',sort_trace,'TickLength',[0 0],'LineWidth',2)
        set(gca,'XLim',[0.8,4.2],'YLim',[0,clu_num+1],'XTick',[])
        set(gca,'XTick',1:4,'XTickLabels',{'L','M','B','UV'})
%         pbaspect([1,2,1])
        title(data(datas).figure_name,'Interpreter','None')
%         axis tight
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 5 10])
        file_path = strjoin({'AverageGainsCluster',data(datas).name,'.png'},'_');
        print(fullfile(fig_path,file_path),'-dpng','-r600')
    end
    autoArrangeFigures
end
%% Find the gain patterns

if contains(data(1).name,'p17b')
    
    close all
    % allocate memory to store the matrices
    type_cell = cell(num_data+1,3);
    
%     %define the zero threshold
%     zero_threshold = 1e-3;
    % 
    % for all of the datasets
    for datas = 1:num_data
        
        % use the gains
        delta_norm = data(datas).delta_norm;
        % get the 10th percentile
        zero_threshold = prctile(abs(delta_norm),10,1);
%         zero_threshold = prctile(abs(delta_norm),5,1);
        % zero the values below a threshold
        delta_norm(abs(delta_norm)<zero_threshold&abs(delta_norm)>0) = 0;
        % turn negatives into -1 and positives into 1
        delta_norm(delta_norm>0) = 1;
        delta_norm(delta_norm<0) = -1;

%         % Use the raw data
%         delta_norm = reshape(data(datas).conc_trace,[],data(datas).time_num,data(datas).stim_num);
%         % get the p2p deflection
%         [max_val,max_idx] = max(delta_norm,[],2);
%         [min_val,min_idx] = min(delta_norm,[],2);
%         
%         p2p = max_val-min_val;
%         % determine the null kernels
%         sd_matrix = squeeze(std(delta_norm(:,1:20,:),0,2));
%         null_kernels = squeeze(p2p)<(10.*sd_matrix);
%         % get the on-off classification
%         % allocate memory for the allocation
%         on_off_matrix = zeros(size(max_idx,1),size(max_idx,3));
%         on_off_matrix(squeeze(max_idx)>squeeze(min_idx)) = 1;
%         on_off_matrix(squeeze(max_idx)<squeeze(min_idx)) = -1;
%         on_off_matrix(null_kernels) = 0;
%         delta_norm = on_off_matrix;
        
        
        
%         % sort rows and plot
%         figure
%         imagesc(sortrows(unique(delta_norm,'rows')))
        
        % quantify the occurrence of each pattern
%         delta_norm(:,3) = 0;
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
        set(gca,'XLim',[-0.2 size(pattern_full,1)])
        set(gca,'YScale','linear','XTick',[],'Visible','off')

        subplot(2,1,1)
%         bar((pattern_counts))
        BarPlotBreak(pattern_counts,pattern_counts(2)*1.8,pattern_counts(1)*0.9,'Line',0.6,2)
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
        
        h = style_figure(gcf,fig_set);
        
    end
end
%% Load the Zhou et al. data

% load the Zhou data
ref_data = load(paths.reference_path);
%% Plot the ref results
% close all
% isolate the kernels for the color channels
kernel_matrix = cat(3,ref_data.AK_R_Mat,ref_data.AK_G_Mat,ref_data.AK_B_Mat,ref_data.AK_UV_Mat);

% get only the dorsal retina ones (ventral FOV)
% dorsal_bool = ref_data.Pos_Mat(:,1)>=1.5&ref_data.Pos_Mat(:,1)<=2.5;
dorsal_bool = ref_data.Pos_Mat(:,1)>=1&ref_data.Pos_Mat(:,1)<=3;
% dorsal_bool = ones(size(ref_data.Pos_Mat(:,1),1),1)==1;
kernel_matrix = kernel_matrix(dorsal_bool,:,:);

% calculate the 0 kernels based on the 10SD criterion used in the ref
% calculate the SD
sd_matrix = squeeze(std(kernel_matrix(:,1:150,:),0,2));

% classify in on and off

% get the maxima and minima
[max_val,max_idx] = max(kernel_matrix,[],2);
[min_val,min_idx] = min(kernel_matrix,[],2);

% calculate the null kernels
null_kernels = squeeze(max_val-min_val)<(5.*sd_matrix);
% null_kernels = sd_matrix<10;

% allocate memory for the allocation
on_off_matrix = zeros(size(max_idx,1),size(max_idx,3));
on_off_matrix(squeeze(max_idx)>squeeze(min_idx)) = 1;
on_off_matrix(squeeze(max_idx)<squeeze(min_idx)) = -1;
on_off_matrix(null_kernels) = 0;

[pattern_ref,ia,ic] = unique(on_off_matrix,'rows');
        
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
pattern_ref = pattern_ref(sort_idx,:);

% allocate memory for the colors
pattern_full = zeros(size(pattern_ref,1),4,3);
% transform the indexes into colors
for channel = 1:3
    pattern_full(pattern_ref(:,channel)==1,channel,channel) = 1;
    pattern_full(pattern_ref(:,channel)==0,channel,:) = 1;
    if channel == 1
        pattern_full(pattern_ref(:,4)==1,4,[1 3]) = 1;
        pattern_full(pattern_ref(:,4)==0,4,:) = 1;
    end
    
end

% store the matrix
type_cell{3,1} = pattern_ref;
type_cell{3,2} = pattern_counts./sum(pattern_counts);
type_cell{3,3} = pattern_full; 
figure
subplot(2,1,2)
image(permute(pattern_full,[2 1 3]))
set(gca,'TickLength',[0 0],'YTick',1:4,'YTickLabels',{'R','G','B','UV'})
%         imagesc(pattern')
subplot(2,1,1)
bar(pattern_counts,'FaceColor',[0.8 0.8 0.8])
% BarPlotBreak(pattern_counts,pattern_counts(2)*1.8,pattern_counts(1)*0.9,'Line',0.6,2)
set(gca,'YScale','linear')
set(gca,'TickLength',[0 0])
%         set(gca,'XDir','reverse')
axis tight

   
% create the settings
fig_set = struct([]);

fig_set(1).fig_path = fig_path;
fig_set(1).fig_name = 'responseTypesZhou.eps';
fig_set(1).fig_size = 4.4;
fig_set(2).fig_size = 4.4

h = style_figure(gcf,fig_set);

autoArrangeFigures
%% Compare the datasets
if contains(data(1).name,'p17b')
    close all
    
    % get the combinations
    comb_vector = nchoosek(1:3,2);
    
    % get the number of combinations
    num_comb = size(comb_vector,1);
    
    % for all the combinations
    for combs = 1:num_comb
        figure
        % get the correlation components
        [corr_1,idx1] = sortrows(type_cell{comb_vector(combs,1),1});
        [corr_2,idx2] = sortrows(type_cell{comb_vector(combs,2),1});
        % get the correlation
        corr_matrix = corr(corr_1',corr_2');
        % filter the matrix
        sorted_count1 = log(type_cell{comb_vector(combs,1),2}(idx1));
        sorted_count2 = log(type_cell{comb_vector(combs,2),2}(idx2));
        
        % sort the patterns
        sorted_pattern1 = type_cell{comb_vector(combs,1),3}(idx1,:,:);
        sorted_pattern2 = type_cell{comb_vector(combs,2),3}(idx2,:,:);
        
        % multiply the rows and columns based on their quantities
        corr_matrix = corr_matrix.*sorted_count1;
        corr_matrix = (corr_matrix'.*sorted_count2)';
        % then filter
        corr_matrix(corr_matrix<20) = 0;
        
        
        % calculate the correlation matrix and plot
        subplot(20,1,1:18)
        imagesc(corr_matrix)
        set(gca,'TickLength',[0 0])
        %     colorbar
        
        axis equal
        axis tight
        
        subplot(20,1,19:20)
        image(permute(sorted_pattern2,[2 1 3]))
        axis equal
        axis tight
    end
    autoArrangeFigures
end
%% Plot ROIs per fish

close all

% allocate memory for the data
fish_rois = cell(num_data,2);
% for both datasets
for datas = 1:num_data
    % get the fish ori info
    fish_ori = data(datas).fish_ori(:,1);
    % get the number of fish
    fish_num = length(unique(fish_ori));
    % get the counts
    fish_rois{datas,2} = histcounts(fish_ori,fish_num);
    % store the name
    fish_rois{datas,1} = data(datas).figure_name;
    
end

tb = table(fish_rois(:,2));
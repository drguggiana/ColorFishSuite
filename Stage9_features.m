% Calculate calcium response standard features per color

% load the data
clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'Features\');

data = load_clusters(cluster_path);
%% Calculate overall average features

close all
% get the number of data sets
num_data = size(data,2);
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
title(strjoin({'Average_Trace',data(:).name},'_'),'Interpreter','None')
file_path = strjoin({'AverageTrace',data(:).name,'.png'},'_');
saveas(gcf, fullfile(fig_path,file_path), 'png')
%% Calculate averages per area

close all
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
        % for all the areas
        for area = 1:area_number
            % get the average trace
            ave_trace = nanmean(data(datas).conc_trace(data(datas).anatomy_info(:,1)==area_list(area),:),1);
            std_trace = nanstd(data(datas).conc_trace(data(datas).anatomy_info(:,1)==area_list(area),:),0,1);
            % plot it
            subplot(area_number,1,area)
            shadedErrorBar(1:size(ave_trace,2),ave_trace,std_trace)
            ylabel(reg_label{area})
        end
        sgtitle(strjoin({'Average_Trace_perArea',data(datas).name},'_'),'Interpreter','None')
        file_path = strjoin({'Average_Trace_perArea',data(datas).name,'.png'},'_');
        saveas(gcf, fullfile(fig_path,file_path), 'png')
    end
    autoArrangeFigures
end
%% Calculate max response, delay to peak

close all
% define the number of parameters
param_num = 3;
param_label = {'LogMax_response','Delay_to_peak','Abs_Logmean_response'};
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
        violin(squeeze(calcium_matrix(:,:,param)));
        hold on
%         plotSpread(squeeze(calcium_matrix(:,:,param)))
        sgtitle(strjoin({param_label{param},data(datas).name},'_'),'Interpreter','None')
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
                plotSpread(plot_matrix);
                ylabel(reg_label{area})
            end
            sgtitle(strjoin({param_label{param},'perRegion',data(datas).name},'_'),'Interpreter','None')
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
            sgtitle(strjoin({param_label{param},'perRegionCluster',data(datas).name},'_'),'Interpreter','None')
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
        plotSpread(plot_matrix);
%         ylabel(reg_label{area})
        
        sgtitle(strjoin({'Gain',data(datas).name},'_'),'Interpreter','None')
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
        sgtitle(strjoin({'GainPerRegion',data(datas).name},'_'),'Interpreter','None')
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
        sgtitle(strjoin({'GainPerRegionCluster',data(datas).name},'_'),'Interpreter','None')
        file_path = strjoin({'GainPerRegionCluster',data(datas).name,'.png'},'_');
        saveas(gcf, fullfile(fig_path,file_path), 'png')
    end
    autoArrangeFigures
end
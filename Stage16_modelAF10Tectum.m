%% Model the tectal data with the AF10 data

clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'Model\');

data = load_clusters(cluster_path);

% define the color scheme
color_scheme = [1 0 0;0 1 0;0 0 1;1 0 1];
%% Run the modelling 
% TODO: use the clusters per region
close all
% % define the regions to use
% tectum_regions = {'L-TcN','R-TcN','L-TcP','R-TcP','L-Cb','R-Cb','L-Hb','R-Hb','L-Pt','R-Pt'};
% tectum_numbers = 1:10;
% af_regions = {'AF4','AF5','AF8','AF9','AF10'};
% af_numbers = [1 2 5 6 7];
% 
% num_datasets = size(data,2);

% % allocate memory to store the model data
% model_cell = cell(size(data,2),1);
% % allocate memory to store the region information
% region_cell = cell(num_datasets,2);
% % for all the data sets
% for datas = 1:num_datasets
%     % get the region info
%     region_info = data(datas).anatomy_info(:,1);
%     % define the regions to be considered (depending on the stimulus protocol)
%     if contains(data(datas).name,{'syn','Syn'})
%         region_list = af_regions;
%         region_numbers = af_numbers;
%     else
%         region_list = tectum_regions;
%         region_numbers = tectum_numbers;
%     end
%     % get the region numbers 
%     num_data = length(region_numbers);
% 
%     %get all the pairwise combinations of the regions
%     region_comb = [nchoosek(region_numbers,2);fliplr(nchoosek(region_numbers,2))];
% 
%     %get the number of combs
%     num_comb = size(region_comb,1);
% 
%     
%     % allocate memory for the model data within this dataset
%     current_models = cell(num_comb,1);
% 
%     % define the period of interest (0 pre, 1 stim, 2 post, 3 pre-post)
%     period = 1;
%     % get the target period labeled with ones
%     rest_all = period_of_interest(period,data(datas).stim_num,1);
%     %for all the combinations
%     for combs = 1:num_comb
%         %concatenate the clusters from both animals involved
% 
%         tar1 = data(datas).region_clusters(region_comb(combs,1)).clu_ave;
%         tar2 = data(datas).region_clusters(region_comb(combs,2)).clu_ave;
%         % if either of them is empty, put a nan in the cell and skip
%         if isempty(tar1) || isempty(tar2)
%             current_models{combs} = NaN;
%             continue
%         end
        % define the fitter (AF10 in this case)
        % get the vector with region info
        region_vector2 = data(2).anatomy_info(:,1);
        region_vector2 = region_vector2==10;
        idx_clu2 = data(2).idx_clu;
        idx_clu2(~region_vector2) = 0;
        % get the traces
        traces2 = data(2).conc_trace;
        traces2 = reshape(traces2,[],data(2).time_num,data(2).stim_num);
        traces2 = reshape(traces2(:,21:60,:),size(traces2,1),[]);
        % allocate memory for the cluster averages
        tar2 = zeros(data(2).clu_num,size(traces2,2));
        % for all the clusters
        for clu = 1:data(2).clu_num
            tar2(clu,:) = mean(traces2(idx_clu2==clu,:),1);
        end
        
%         tar2 = data(2).region_clusters(7).clu_ave;
        % define the fitted (Tectum)
%         tar1 = cat(1,data(1).region_clusters([1,3]).clu_ave);
        % get the vector with region info
        region_vector1 = data(1).anatomy_info(:,1);
        region_vector1 = region_vector1==1|region_vector1==3;
        idx_clu1 = data(1).idx_clu;
        idx_clu1(~region_vector1) = 0;
        % get the traces
        traces1 = data(1).conc_trace;
        traces1 = reshape(traces1,[],data(1).time_num,data(1).stim_num);
        traces1 = reshape(traces1(:,21:60,:),size(traces1,1),[]);
        % allocate memory for the cluster averages
        tar1 = zeros(data(1).clu_num,size(traces1,2));
        % for all the clusters
        for clu = 1:data(1).clu_num
            tar1(clu,:) = mean(traces1(idx_clu1==clu,:),1);
        end
        
        %for all the averages in 1, calculate models from the raw traces in 2
        %allocate memory to store the model results
        model_para = cell(size(tar1,1),1);
        %for all the averages
        for clu = 1:size(tar1,1)
            disp(strcat('Current clu: ',num2str(clu)))


            model_para{clu} = fitrlinear(tar2',tar1(clu,:)','CrossVal','on');

        end
%         % fill up the dataset cell
%         current_models{combs} = model_para;
%     end
%     % fill up the overall cell
%     model_cell{datas} = current_models;
% end
%% Plot the fit

close all
% get the fit qualities

% get the number of tectal clusters
tectum_clunum = size(tar1,1);
% allocate memory for the losses
loss_vector = zeros(tectum_clunum,1);
% for all the clusters
for clu = 1:tectum_clunum
    loss_vector(clu) = 1-kfoldLoss(model_para{clu});
end

% sort by quality
[sorted_loss,sort_idx] = sort(loss_vector,'descend');
% sort the models
sorted_model = model_para(sort_idx);

% plot the fit quality
figure
plot(sorted_loss,'-ko','LineWidth',1)
% create the settings
fig_set = struct([]);

fig_set(1).fig_path = fig_path;
fig_set(1).fig_name = strjoin({'fitqualAF10Tectum.eps'},'_');
fig_set(1).fig_size = 3.6;

fig_set(1).box = 'off';
fig_set(1).painters = 1;
h = style_figure(gcf,fig_set);

% % sort the traces too
% sorted_traces = tar1(sort_idx,:);

% plot the traces and the fit quality
figure
%     fig('units','centimeters','width',10,'height',20)
a_count = 1;
% get the clusters indexes
idx_clu = cat(1,data(1).region_clusters(1).idx_clu,...
    data(1).region_clusters(3).idx_clu+data(1).region_clusters(1).clu_num);

% allocate memory for the new idx
new_idx = zeros(size(idx_clu));
% reorder the indexes according to the sorting
for clu = 1:tectum_clunum
    new_idx(idx_clu==clu) = find(sort_idx==clu);
end
% rewrite the old idx
idx_clu = new_idx;

% get the indexes of the regions from the main conc_trace
region_idx = (data(1).anatomy_info(:,1)==1)|(data(1).anatomy_info(:,1)==3);
% get the raw traces
% conc_trace = data(1).conc_trace;
conc_trace = data(1).conc_trace(region_idx,:);

% get the number of stimuli
stim_num = data(1).stim_num;
% % get the number of clusters
% clu_num = data(1).clu_num;
% define the trace offset
trace_offset = 3;
% calculate the top of the plot
plot_top = trace_offset*(tectum_clunum-1);
% framerate added manually, need to fix this
framerate = data(1).framerate;
% for all the clusters
for clu = 1:tectum_clunum
    % get the average trace
    ave_trace = normr_1(nanmean(conc_trace(idx_clu==clu,:),1),1);
    std_trace = nanstd(conc_trace(idx_clu==clu,:),0,1);
    ave_perstim = reshape(ave_trace,[],stim_num);
    std_perstim = reshape(std_trace,[],stim_num);
    time_vector = (0:size(ave_trace,2)-1)./framerate;
    time_perstim = reshape(time_vector,[],stim_num);
    
    % get the fitted traces
    fit_trace = normr_1(kfoldPredict(sorted_model{clu}),1);
    fit_perstim = reshape(fit_trace,[],stim_num);
    % add the blank periods as they are not fitted
    fit_perstim = cat(1,zeros(20,4),fit_perstim,zeros(20,4));
    
    % split by stimulus
    for stim = 1:stim_num
        % plot it
        shadedErrorBar(time_perstim(:,stim),...
            ave_perstim(:,stim)+(plot_top-(a_count-1)*trace_offset),std_perstim(:,stim),...
            {'color',color_scheme(stim,:),'LineWidth',1})
%         plot(time_perstim(:,stim),...
%             ave_perstim(:,stim)+(plot_top-(a_count-1)*trace_offset),'Color',color_scheme(stim,:),'LineWidth',1)
        hold on
        plot(time_perstim(:,stim),fit_perstim(:,stim)+(plot_top-(a_count-1)*trace_offset),...
            '-','Color',[0 0 0],'LineWidth',1)
%             '--','Color',tint_colormap(color_scheme(stim,:),0.5),'LineWidth',1)
        
    end
    
    % update the counter
    a_count = a_count + 1;
    
end
axis tight
%     % plot the intermediate lines
%     for stim = 1:stim_num-1
%         plot([time_perstim(end,stim),time_perstim(end,stim)],...
%             get(gca,'YLim'),'k','LineWidth',2)
%     end
%     pbaspect([1,2,1])
axis tight
set(gca,'YTick',0:trace_offset:(a_count-2)*trace_offset,'YTickLabels',string(tectum_clunum:-1:1))
set(gca,'TickLength',[0 0])
xlabel('Time (s)')
%     sgtitle(strjoin({'Average_Trace_perArea',data(datas).name},'_'),'Interpreter','None')
%         set(gca,'XLim',[0,time_perstim(end,end)])
box off



% title(data(1).figure_name,'Interpreter','None')
%     set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 5 10])
% assemble the figure path
set(gca,'FontSize',10,'LineWidth',2)
%     file_path = strjoin({'clusterTraces',data(datas).name,'.png'},'_');
%     print(fullfile(fig_path,file_path),'-dpng','-r600')
% create the settings
fig_set = struct([]);

fig_set(1).fig_path = fig_path;
fig_set(1).fig_name = strjoin({'modelAF10Tectum.eps'},'_');
if contains(data(1).name,'syn')
    fig_set(1).fig_size = [5.3 9.8];
else
    fig_set(1).fig_size = [5.3 13.8];
end

fig_set(1).box = 'off';
fig_set(1).painters = 1;
h = style_figure(gcf,fig_set);


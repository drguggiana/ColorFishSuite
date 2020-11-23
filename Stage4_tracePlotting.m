clearvars
close all

load('paths.mat')
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'Traces\');


data = load_clusters(cluster_path);
% define the color scheme depending on the stimulus type
if contains(data(1).name,'p17')
    color_scheme = [1 0 0;0 1 0;0 0 1;1 0 1];
    dataset_labels = {'RAs','RGCs'};
else
    color_scheme = distinguishable_colors(6);
    dataset_labels = {'Tectum','AF10'};
end
%% Plot all the trial averaged traces

close all

for datas = 1:length(data)
   
%     sort_traces = sort_by_average(data(datas), 'conc_trace');
    [sort_idx,sorted_traces] = sort(data(datas).idx_clu);
    figure
    
    subplot(1,20,1:19)
    set(gcf,'color','w')
    imagesc(normr_1(data(datas).conc_trace(sorted_traces,:),0))
%     title(dataset_labels{datas},'Interpreter','None')
    set(gca,'YTick',[1,length(data(datas).idx_clu)],...
        'XTick',0:(100/0.952):size(data(datas).conc_trace,2),...
        'XTickLabel',0:100:size(data(datas).conc_trace,2)/0.952)
%     xlabel('Time (s)')
%     ylabel('ROI')
    set(gca,'TickLength',[0 0],'LineWidth',0.5,'FontSize',12)
    axis tight
    
    % create the settings
    fig_set = struct([]);
    
    fig_set(1).fig_path = fig_path;
    fig_set(1).fig_name = strjoin({'traces',data(datas).name,'.eps'},'_');
    fig_set(1).box = 'on';
    fig_set(1).LineWidth = 0.05;
%     fig_set(1).xlabel = 'Wavelength (nm)';
%     fig_set(1).ylabel = 'Abs/Em  ';

    switch data(datas).name
        case {'p17b_syngc6s','p17b_gc6s'}
            fig_set(1).fig_size = [11 3.5];
            fig_set(1).colorbar_label = 'Normalized Delta F/F';
            fig_set(1).colorbar = 1;
        case {'p8_gc6s','p8_SynG6s'}
            fig_set(1).fig_size = [8.5 3.2];
            fig_set(1).colorbar_label = 'Normalized Delta F/F';
            fig_set(1).colorbar = 1;
        otherwise
            fig_set(1).fig_size = [5 3];
%             fig_set(1).colorbar_label = 'Norm. Delta F/F';
    end

%     if contains(data(datas).name,'p17')
%         fig_set(1).fig_size = [11 3.5];
%     else
%         fig_set(1).fig_size = [8.5 3.2];
%     end
%     fig_set(1).crop = 0;
    
    h = style_figure(gcf,fig_set);

end


% autoArrangeFigures
%% Plot the cluster averages

close all

for datas = 1:length(data)
   
    sort_clu = sort_by_average(data(datas), 'clu_ave');

    figure
    imagesc(normr_1(data(datas).clu_ave(sort_clu,:),0))
    title(strjoin({'Clusters',data(datas).figure_name},' '),'Interpreter','None')
    set(gca,'FontSize',15,'TickLength',[0 0])
    % assemble the figure path 
    file_path = fullfile(fig_path,strjoin({'clusters',data(datas).name,'.png'},'_'));
%     saveas(gcf, fullfile(fig_path,file_path), 'png')
%     set(gcf,'PaperUnits','centimeters','Position',)
    
    set(gcf,'Color','w')
    export_fig(file_path,'-r600')
end


% autoArrangeFigures
%% Plot the cluster averages as traces

close all

for datas = 1:length(data)
    figure
%     fig('units','centimeters','width',10,'height',20)
    a_count = 1;
    % get the clusters indexes
    idx_clu = data(datas).idx_clu;
    % get the raw traces 
    conc_trace = data(datas).conc_trace;
    % get the number of stimuli
    stim_num = data(datas).stim_num;
    % get the number of clusters
    clu_num = data(datas).clu_num;
    % define the trace offset
    trace_offset = 5;
    % calculate the top of the plot
    plot_top = trace_offset*(clu_num-1);
    % framerate added manually, need to fix this
    framerate = data(datas).framerate;
    % for all the clusters
    for clu = 1:clu_num
        % get the average trace
        ave_trace = nanmean(conc_trace(idx_clu==clu,:),1);
        std_trace = nanstd(conc_trace(idx_clu==clu,:),0,1);
        ave_perstim = reshape(ave_trace,[],stim_num);
        std_perstim = reshape(std_trace,[],stim_num);
        time_vector = (0:size(ave_trace,2)-1)./framerate;
        time_perstim = reshape(time_vector,[],stim_num);
        % split by stimulus
        for stim = 1:stim_num
            % plot it
            shadedErrorBar(time_perstim(:,stim),...
                ave_perstim(:,stim)+(plot_top-(a_count-1)*trace_offset),std_perstim(:,stim),...
                {'color',color_scheme(stim,:),'LineWidth',1})
            hold on
        end
        text(280,(plot_top-(a_count-1)*trace_offset)+3,num2str(sum(idx_clu==clu)),...
            'FontSize',7,'FontName','Arial','Color',[0.5 0.5 0.5])
        % update the counter
        a_count = a_count + 1;
        
    end
    
    
    % assemble the labels based on the dataset
    ylabels = cell(clu_num,1);
    
    % for all the clusters
    for clu = 1:clu_num
        % select the prefix depending on dataset
        switch data(datas).name
            case 'p17b_syngc6s'
                prefix = 'RGC';
            case 'p17b_gc6s'
                prefix = 'RA';
            otherwise
                prefix = '';
        end
        ylabels{clu} = strcat(prefix,num2str(clu_num-clu+1));
    end
%     % plot the intermediate lines
%     for stim = 1:stim_num-1
%         plot([time_perstim(end,stim),time_perstim(end,stim)],...
%             get(gca,'YLim'),'k','LineWidth',2)
%     end
%     pbaspect([1,2,1])
    axis tight
    set(gca,'YTick',0:trace_offset:(a_count-2)*trace_offset,'YTickLabels',ylabels)
    set(gca,'TickLength',[0 0])
    xlabel('Time (s)')
%     sgtitle(strjoin({'Average_Trace_perArea',data(datas).name},'_'),'Interpreter','None')
    %         set(gca,'XLim',[0,time_perstim(end,end)])
    box off
    
    

%     title(data(datas).figure_name,'Interpreter','None')
%     set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 5 10])
    % assemble the figure path 
    set(gca,'FontSize',10,'LineWidth',2)
%     file_path = strjoin({'clusterTraces',data(datas).name,'.png'},'_');
%     print(fullfile(fig_path,file_path),'-dpng','-r600')
        % create the settings
    fig_set = struct([]);
    
    fig_set(1).fig_path = fig_path;
    fig_set(1).fig_name = strjoin({'clusterTraces',data(datas).name,'.eps'},'_');
    if contains(data(datas).name,'syn')
        fig_set(1).fig_size = [5.3 9.8];
    else
        fig_set(1).fig_size = [5.3 13.8];
    end

    fig_set(1).box = 'off';
    fig_set(1).painters = 1;
    h = style_figure(gcf,fig_set);

end


% autoArrangeFigures
%% Plot a selected trace
close all
% define the target traces
n = 660;
target_vector = [n n+1 n+2];

% initialize a variable to lift the traces
trace_lift = 5;
% define the frame rate
framerate = data(1).framerate;

% define the color scheme depending on the stimulus type
if contains(data(1).name,'p17')
    color_scheme = [1 0 0;0 1 0;0 0 1;1 0 1];
else
    color_scheme = distinguishable_colors(6);
end

figure

% initialize a trace counter
trace_counter = 0;

% for all the datasets
for datas = 1:length(data)
    % for all the target traces
    for target_trace = target_vector
        % get the trace
        trace = data(datas).conc_trace(target_trace,:);
        trace_perstim = reshape(trace,[],data(datas).stim_num);
        % get the time vector
        time_perstim = reshape((0:length(trace)-1)./framerate,[],data(datas).stim_num);
        % for all the stimuli
        for stim = 1:data(datas).stim_num
            plot(time_perstim(:,stim),trace_perstim(:,stim)+(trace_counter)*trace_lift,...
                'color',color_scheme(stim,:),'LineWidth',1)
            hold on
        end
        % update the trace counter
        trace_counter = trace_counter + 1;
    end
end
% % plot the intermediate lines
% for stim = 1:data(1).stim_num-1
%     plot([time_perstim(end,stim),time_perstim(end,stim)],...
%         get(gca,'YLim'),'k','LineWidth',1)
% end
set(gca,'XColor','w','YColor','w')
set(gca,'XLim',[0,time_perstim(end,end)])
set(gca,'YTick',linspace(0,(length(data)-1)*trace_lift,length(data)),...
    'YTickLabels',{data.figure_name},'TickLabelInterpreter','None')
% box off
% set(gca,'TickLength',[0 0],'LineWidth',2)
xlabel('Time (s)')
% pbaspect([2,1,1])
axis tight
% set(gca,'FontSize',15)

% create the settings
fig_set = struct([]);

fig_set(1).fig_path = fig_path;
fig_set(1).fig_name = strjoin({'SingleTrace',data.name,'.eps'},'_');
fig_set(1).fig_size = [5.7 4.3];

h = style_figure(gcf,fig_set);
%% Clusters proportion per area

% plot a matrix indicating how many instances of a cluster are in each
% region

% only do this for the p17 dataset
if contains(data(1).name,'p17b')
    close all
%     % define the region labels
%     region_labels = {'Red','Green','Blue','UV'};
    % define the fontsize
    fontsize = 10;
    
    % for all the datasets
    for datas = 1:length(data)
        figure
        
        % get the anatomy info
        anatomy_info = data(datas).anatomy_info;

        % define the labels and exclude AF 6 and 7 cause too few terminals
        switch data(datas).figure_name
            case {'RAs','Tectum'}
                region_labels = {'L-TcN','R-TcN','L-TcP','R-TcP','L-Cb','R-Cb','L-Hb','R-Hb','L-Pt','R-Pt'};
            case {'RGCs','AF10'}
                region_labels = {'AF4','AF5','AF8','AF9','AF10'};
                anatomy_info(anatomy_info(:,1)==6) = NaN;
                anatomy_info(anatomy_info(:,1)==7) = NaN;
        end
        % get the list and number of regions
        region_list = unique(anatomy_info(:,1));
        region_list = region_list(~isnan(region_list));
        region_num = length(region_list);
        % get the clusters
        idx_clu = data(datas).idx_clu;
        % get the number of clusters
        clu_num = data(datas).clu_num;
        
        % allocate memory for the matrix
        cluster_perregion = zeros(clu_num,region_num);
        % go through all of the regions and clusters filling up the matrix
        for clu = 1:clu_num
            for region = 1:region_num
                cluster_perregion(clu,region) = sum(anatomy_info(:,1)==region_list(region) & ...
                    idx_clu==clu);
            end
        end
        imagesc(log(normr_1(cluster_perregion,2)))
        set(gca,'TickLength',[0 0])
        switch data(datas).figure_name
            case {'RAs','Tectum'}
%                 set(gca,'XTick',1:2:region_num,'XTickLabels',{'TcN','TcP','Pt','Hb','Cb'},'FontSize',fontsize,...
%                     'XTickLabelRotation',45)
%                 set(gca,'XTick',1:region_num,'XTickLabels',region_labels,'FontSize',fontsize,...
%                     'XTickLabelRotation',45)
                set(gca,'XTick',[])
            case {'RGCs','AF10'}
                set(gca,'XTick',[])
%                 set(gca,'XTick',1:region_num,'XTickLabels',region_labels,'FontSize',fontsize,...
%                     'XTickLabelRotation',45)
        end
        set(gca,'YTick',1:2:clu_num,'FontSize',fontsize)
        axis square
%         title(data(datas).figure_name)
%         set(gca,'FontSize',15,'LineWidth',2)
%         set(gcf,'Color','w')
%         colormap(magma)
%         cba = colorbar;
%         set(cba,'TickLength',0,'LineWidth',2)
%         ylabel(cba,'Log Fraction of Traces')
        % assemble the figure path
%         file_path = strjoin({'clusterPerArea',data(datas).name,'.png'},'_');
%         file_path = fullfile(fig_path,file_path);
%         saveas(gcf, fullfile(fig_path,file_path), 'png')
%         print(file_path,'-dpng','-r600') 
%         export_fig(file_path,'-r600')
        
        % create the settings
        fig_set = struct([]);
        
        fig_set(1).fig_path = fig_path;
        fig_set(1).fig_name = strjoin({'clusterPerArea',data(datas).name,'.eps'},'_');
        fig_set(1).fig_size = 5;
        fig_set(1).colorbar = 1;
        fig_set(1).colorbar_label = 'Log(Fraction Traces)';
        fig_set(1).box = 'on';
        fig_set(1).LineWidth = 0.05;
        
        h = style_figure(gcf,fig_set);

    end
%     outerbar(file_path,gca,'Fraction traces/cluster');

    autoArrangeFigures
end
%% Plot the raw traces per area

close all

% for all datasets
for datas = 1:length(data)
    % load the traces
    conc_trace = data(datas).conc_trace;
    
    % load the anatomy info
    anatomy_info = data(datas).anatomy_info;
    
    % define the labels and exclude AF 6 and 7 cause too few terminals
    switch data(datas).figure_name
        case {'RAs','Tectum'}
            region_labels = {'L-TcN','R-TcN','L-TcP','R-TcP','L-Cb','R-Cb','L-Hb','R-Hb','L-Pt','R-Pt'};
        case {'RGCs','AF10'}
            region_labels = {'AF4','AF5','AF8','AF9','AF10'};
            anatomy_info(anatomy_info(:,1)==6) = NaN;
            anatomy_info(anatomy_info(:,1)==7) = NaN;
    end
    % load the cluster info
    idx_clu = data(datas).idx_clu;
    % get a list of regions
    region_list = unique(anatomy_info(:,1));
    % exclude NaNs
    region_list = region_list(~isnan(region_list));
    % get the number of regions
    region_number = length(region_list);
    
    % for all the regions
    for region = 1:region_number
        
        figure
        % get the indexes and traces for this region
        idx_region = idx_clu(region_list(region)==anatomy_info(:,1));
        conc_region = conc_trace(region_list(region)==anatomy_info(:,1),:);
        
        % sort the traces by the index
%         [~,sort_idx] = sort(idx_region);
%         sorted_traces = conc_region(sort_idx,:);
        sorted_traces = sort_traces(conc_region);
        
        % plot
        imagesc(normr_1(sorted_traces,0))
        colormap(magma)
        set(gca,'LineWidth',2,'TickLength',[0 0],'FontSize',15)
        set(gca,'YTick',[1,length(idx_region)],'XTick',[])
%             'XTick',0:(100/data(datas).framerate):size(sorted_traces,2),...
%             'XTickLabel',0:100:size(sorted_traces,2)/data(datas).framerate)
        set(gca,'YTickLabelRotation',90)
%         xlabel('Time (s)')
%         ylabel('ROIs')
%         title(region_labels{region},'FontSize',7)
        axis tight
%         file_path = strjoin({'regionTraces',data(datas).name,region_labels{region},'.png'},'_');
%         print(fullfile(fig_path,file_path),'-dpng','-r600')
        
        % create the settings
        fig_set = struct([]);
        
        fig_set(1).fig_path = fig_path;
        fig_set(1).fig_name = strjoin({'regionTraces',data(datas).name,region_labels{region},'.eps'},'_');
        fig_set(1).box = 'on';
        fig_set(1).colorbar = 1;
        fig_set(1).LineWidth = 0.05;
%         fig_set(1).colorbar_label = 'Delta F/F';
        fig_set(1).fig_size = [4 3];
        
        h = style_figure(gcf,fig_set);
    end
    
    

    
end
autoArrangeFigures
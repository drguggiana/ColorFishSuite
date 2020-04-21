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
else
    color_scheme = distinguishable_colors(6);
end
%% Plot all the trial averaged traces

close all

for datas = 1:length(data)
   
%     sort_traces = sort_by_average(data(datas), 'conc_trace');
    [sort_idx,sort_traces] = sort(data(datas).idx_clu);
    figure
    
    subplot(1,20,1:19)
    set(gcf,'color','w')
    imagesc(normr_1(data(datas).conc_trace(sort_traces,:),0))
    title(data(datas).figure_name,'Interpreter','None')
    set(gca,'YTick',1:1000:length(data(datas).idx_clu),...
        'XTick',0:(30/0.952):size(data(datas).conc_trace,2),...
        'XTickLabel',0:30:size(data(datas).conc_trace,2)/0.952)
    xlabel('Time (s)')
    ylabel('ROI')
    set(gca,'TickLength',[0 0])
    
    subplot(1,20,20)
    imagesc(sort_idx)
    colormap(gca,'colorcube')
    set(gca,'XTick',[],'YTick',[])
    
    file_path = strjoin({'traces',data(datas).name,'.png'},'_');
    saveas(gcf, fullfile(fig_path,file_path), 'png')
end


autoArrangeFigures
%% Plot the cluster averages

close all

for datas = 1:length(data)
   
    sort_clu = sort_by_average(data(datas), 'clu_ave');

    figure
    imagesc(normr_1(data(datas).clu_ave(sort_clu,:),0))
    title(strjoin({'Clusters',data(datas).figure_name},' '),'Interpreter','None')
    % assemble the figure path 
    file_path = strjoin({'clusters',data(datas).name,'.png'},'_');
    saveas(gcf, fullfile(fig_path,file_path), 'png')
end


autoArrangeFigures
%% Plot the cluster averages as traces

close all

for datas = 1:length(data)
   
    fig('units','centimeters','width',10,'height',20)
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
                ave_perstim(:,stim)+(a_count-1)*trace_offset,std_perstim(:,stim),...
                {'color',color_scheme(stim,:),'LineWidth',1})
            hold on
        end
        
        % update the counter
        a_count = a_count + 1;
        
    end
    axis tight
    % plot the intermediate lines
    for stim = 1:stim_num-1
        plot([time_perstim(end,stim),time_perstim(end,stim)],...
            get(gca,'YLim'),'k','LineWidth',2)
    end
    pbaspect([1,2,1])
    axis tight
    set(gca,'YTick',0:trace_offset:(a_count-2)*trace_offset,'YTickLabels',string(1:clu_num))
    set(gca,'TickLength',[0 0])
    xlabel('Time (s)')
%     sgtitle(strjoin({'Average_Trace_perArea',data(datas).name},'_'),'Interpreter','None')
    %         set(gca,'XLim',[0,time_perstim(end,end)])
    box off
    
    

    title(strcat(data(datas).figure_name,'+'),'Interpreter','None')
    % assemble the figure path 
    file_path = strjoin({'clusterTraces',data(datas).name,'.png'},'_');
    saveas(gcf, fullfile(fig_path,file_path), 'png')
end


autoArrangeFigures
%% Plot a selected trace
close all
% define the target trace
target_trace = 1050;

% initialize a variable to lift the traces
trace_lift = 5;
% define the frame rate
% TODO: add to structure
framerate = 1/0.952;

% define the color scheme depending on the stimulus type
if contains(data(1).name,'p17')
    color_scheme = [1 0 0;0 1 0;0 0 1;1 0 1];
else
    color_scheme = distinguishable_colors(6);
end

figure
% for all the datasets
for datas = 1:length(data)
    % get the trace
    trace = data(datas).conc_trace(target_trace,:);
    trace_perstim = reshape(trace,[],data(datas).stim_num);
    % get the time vector
    time_perstim = reshape((0:length(trace)-1)./framerate,[],data(datas).stim_num);
    % for all the stimuli
    for stim = 1:data(datas).stim_num
        plot(time_perstim(:,stim),trace_perstim(:,stim)+(datas-1)*trace_lift,...
            'color',color_scheme(stim,:),'LineWidth',2)
        hold on
    end
end

% plot the intermediate lines
for stim = 1:data(datas).stim_num-1
    plot([time_perstim(end,stim),time_perstim(end,stim)],...
        get(gca,'YLim'),'k','LineWidth',2)
end
set(gca,'XLim',[0,time_perstim(end,end)])
set(gca,'YTick',linspace(0,(length(data)-1)*trace_lift,length(data)),...
    'YTickLabels',{data.figure_name},'TickLabelInterpreter','None')
box off
set(gca,'TickLength',[0 0])
xlabel('Time (s)')
pbaspect([2,1,1])
axis tight

% assemble the figure path
file_path = strjoin({'SingleTrace',data.name,'.png'},'_');
saveas(gcf, fullfile(fig_path,file_path), 'png')


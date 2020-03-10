clearvars
close all

load('paths.mat')
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'Traces\');


data = load_clusters(cluster_path);
%% Plot all the stim average traces and the clusters

close all

for datas = 1:length(data)
   
    sort_traces = sort_by_average(data(datas), 'conc_trace');
    sort_clu = sort_by_average(data(datas), 'clu_ave');
    figure
%     subplot(1,2,1)
    imagesc(normr_1(data(datas).conc_trace(sort_traces,:),0))
    title(strcat('Traces_',data(datas).name),'Interpreter','None')
    file_path = strjoin({'traces',data(datas).name,'.png'},'_');
    saveas(gcf, fullfile(fig_path,file_path), 'png')
%     subplot(1,2,2)
    figure
    imagesc(normr_1(data(datas).clu_ave(sort_clu,:),0))
    title(strcat('Clusters_',data(datas).name),'Interpreter','None')
    % assemble the figure path 
    file_path = strjoin({'clusters',data(datas).name,'.png'},'_');
    saveas(gcf, fullfile(fig_path,file_path), 'png')
end


autoArrangeFigures
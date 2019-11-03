clearvars


close all

cluster_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\Stage3_cluster\';
fig_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Figures\Traces\';


data = load_clusters(cluster_path);
%%

close all

for datas = 1:length(data)
    figure
    subplot(1,2,1)
    imagesc(normr_1(data(datas).conc_trace,0))
    subplot(1,2,2)
    imagesc(normr_1(data(datas).clu_ave,0))
    % assemble the figure path 
    file_path = strjoin({'tracesClusters',data(datas).name,'.png'},'_');
    saveas(gcf, fullfile(fig_path,file_path), 'png')
end


autoArrangeFigures
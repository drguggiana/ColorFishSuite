%% Clean up and load data
clearvars

close all

cluster_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\Stage3_cluster\';

data = load_clusters(cluster_path);
%% Calculate PCA trajectories, Niessing-style

close all

%define the stim labels based on the paradigm
%extract the actual file name
% [~,stim_name,~] = fileparts(name_cell{1});
%scan for the p17b
if contains(data(1).file_name,'p17b')
    %if it's p17b
    stim_labels = {'Red','Green','Blue','UV'};
else %if it's p6p8 instead
    %define the stim labels (for p6p8 data)
    stim_labels = {'Red CK','UV CK','Red GR','UV GR','Red FL','UV FL'};
end

%define the plot colors
plot_col = [1 0 0;0 1 0;0 0 1;1 0 1];

%define the sets of time regions to correlate
% time_corr = (1:40)';
time_corr = [1:10;11:20;21:30;31:40];
% time_corr = [1:5;6:10;11:15;16:20;21:25;26:30;31:35;36:40];

%get the time axis
%get the number of time points
timep_num = size(time_corr,1);
%allocate memory for the axis
timep_axis = zeros(timep_num,1);
%for all the time points
for timep = 1:timep_num
    timep_axis(timep) = mean(time_corr(timep,:));
end
%get the number of time regions
num_times = size(time_corr,1);

% %allocate memory to store the correlations
% tcorr_mat = zeros(num_data+1,num_times,stim_num2,stim_num2);

num_data = size(data,2);

%for both data sets
for datas = 1:num_data
    
    % get the time and stim nums
    time_num = data(datas).time_num;
    stim_num2 = data(datas).stim_num;
    % get the number of animals and the animal info
    fish_ori = data(datas).fish_ori;
    num_animals = length(unique(fish_ori(:,1)));
    %if a normal data set
    if datas <= num_data
        %load the clusters
        conc_trace = data(datas).clu_ave;
        %get the number of traces
        trace_num = size(conc_trace,1);
        %reshape the matrix to calculate correlations between stimuli
        resh_trace = reshape(conc_trace,trace_num,time_num,stim_num2);
        %extract only the stim period
        resh_trace = resh_trace(:,21:60,:);
    else
        %use the stimuli
        resh_trace = cone_stim;
        
        %get the number of traces
        trace_num = size(resh_trace,1);
        %reshape the matrix 
        resh_trace = reshape(resh_trace,trace_num,size(resh_trace,2)/stim_num2,stim_num2);
    end
    
    pca_mat = zeros(size(resh_trace,2),3,stim_num2);
    figure
    %for all the times
    for stim = 1:stim_num2
        % split between animals
        % for all the animals
        for animals = 1:num_animals
            
    %         %get the number of time points per stimulus
    %         t_perstim = size(time_corr,2);
    %         %now reshape again for calculating the correlation across traces for
    %         %each stimulus
    %         corr_trace = reshape(resh_trace(:,time_corr(times,:),:),trace_num*t_perstim,stim_num2);
    %         temp_mat = squeeze(mean(resh_trace(:,time_corr(times,:),:),2));
            temp_mat = resh_trace(fish_ori(:,1)==animals,:,stim)';
            %calculate and plot an allvall corr matrix
    %         temp_mat = corr(corr_trace);
    %         c_lims = [min(temp_corr(:)), max(temp_corr(:))];
            pca_res = pca(temp_mat);
            
            % implement CCA 
            pca_mat(:,:,stim) = temp_mat*pca_res(:,1:3);
        end
        %for all the points
        for points = 1:size(pca_mat,1)
            plot3(pca_mat(points,1,stim),pca_mat(points,2,stim),pca_mat(points,3,stim),...
                'Marker','o','MarkerSize',points/5,'MarkerFaceColor',plot_col(stim,:),...
                'MarkerEdgeColor',plot_col(stim,:))
            hold('on')
        end
        %also plot the lines
        plot3(pca_mat(:,1,stim),pca_mat(:,2,stim),pca_mat(:,3,stim),...
            'Color',plot_col(stim,:))
    end
    %get the file name
%     temp_name = strsplit(data{datas},'\');
%     temp_name = strsplit(temp_name{end},'_');
    title(strcat(data(datas).file_name,'_clusterAve'),'FontSize',20,'Interpreter','None')
    xlabel('PC 1','FontSize',20)
    ylabel('PC 2','FontSize',20)
    zlabel('PC 3','FontSize',20)
    
    
end
%% Plot pca trajectories with the raw data
%for all the fish
for datas = 1:num_data
    
    %show the current fish
    fprintf(strcat('Current fish:',num2str(datas),'\r\n'))
    %load the cluster indexes for this fish
%     idx_clu = load(name_cell{datas},'idx_clu');
%     idx_clu = idx_clu.idx_clu;
%     
%     %also load the raw traces
%     raw_trace = load(name_cell{datas},'conc_trace');
%     raw_trace = raw_trace.conc_trace;

    idx_clu = data(datas).idx_clu;
    raw_trace = data(datas).conc_trace;
    
    resh_trace = reshape(raw_trace,size(raw_trace,1),time_num,stim_num2);
    %clip the edges away
    resh_trace = resh_trace(:,21:60,:);
    
    pca_mat = zeros(size(resh_trace,2),3,stim_num2);
    figure
    %for all the times
    for stim = 1:stim_num2
%         %get the number of time points per stimulus
%         t_perstim = size(time_corr,2);
%         %now reshape again for calculating the correlation across traces for
%         %each stimulus
%         corr_trace = reshape(resh_trace(:,time_corr(times,:),:),trace_num*t_perstim,stim_num2);
%         temp_mat = squeeze(mean(resh_trace(:,time_corr(times,:),:),2));
        temp_mat = resh_trace(:,:,stim)';
        %calculate and plot an allvall corr matrix
%         temp_mat = corr(corr_trace);
%         c_lims = [min(temp_corr(:)), max(temp_corr(:))];
        pca_res = pca(temp_mat);
        pca_mat(:,:,stim) = temp_mat*pca_res(:,1:3);
        %for all the points
        for points = 1:size(pca_mat,1)
            plot3(pca_mat(points,1,stim),pca_mat(points,2,stim),pca_mat(points,3,stim),...
                'Marker','o','MarkerSize',points/5,'MarkerFaceColor',plot_col(stim,:),...
                'MarkerEdgeColor',plot_col(stim,:))
            hold('on')
        end
        %also plot the lines
        plot3(pca_mat(:,1,stim),pca_mat(:,2,stim),pca_mat(:,3,stim),...
            'Color',plot_col(stim,:))
    end
    %get the file name
%     temp_name = strsplit(name_cell{datas},'\');
%     temp_name = strsplit(temp_name{end},'_');
    title(strcat(data(datas).file_name,'_rawTraces'),'FontSize',20,'Interpreter','None')
    xlabel('PC 1','FontSize',20)
    ylabel('PC 2','FontSize',20)
    zlabel('PC 3','FontSize',20)
end
autoArrangeFigures
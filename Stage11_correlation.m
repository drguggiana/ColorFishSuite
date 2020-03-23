clearvars
close all

load('paths.mat')
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'Correlations\');


data = load_clusters(cluster_path);
%% Calculate correlation matrices for each data set

close all

% get the number of dataset
num_data = size(data,2);
% get the number fo stimuli
stim_num2 = data(1).stim_num;
% get the number of time bins
time_num = data(1).time_num;
%define the stim labels based on the paradigm
% %extract the actual file name
% [~,stim_name,~] = fileparts(name_cell{1});

% get the dataset name
stim_name = data(1).name;
%scan for the p17b
if contains(stim_name, 'p17b')
    %if it's p17b
    stim_labels = {'Red','Green','Blue','UV'};
else %if it's p6p8 instead
    %define the stim labels (for p6p8 data)
    stim_labels = {'Red CK','UV CK','Red GR','UV GR','Red FL','UV FL'};
end
% allocate memory to save the labels
correlations = zeros(stim_num2,stim_num2, num_data);

%for both data sets
for datas = 1:num_data
    %if a normal data set
    if datas <= num_data
        %load the clusters
        conc_trace = data(datas).conc_trace;
        %get the number of traces
        trace_num = size(conc_trace,1);
        %reshape the matrix to calculate correlations between stimuli
        resh_trace = reshape(conc_trace,trace_num,time_num,stim_num2);
        %extract only the stim period
        resh_trace = resh_trace(:,21:60,:);
%     else
%         %use the stimuli
%         resh_trace = cone_stim;
%         %get the number of traces
%         trace_num = size(resh_trace,1)/stim_num2;
    end
    
    % get the trace fish of origin
    fish_ori = data(datas).fish_ori;
    % get the number of fish
    num_fish = size(unique(fish_ori(:,1)),1);
    % allocate memory for the output
    corr_perfish = zeros(stim_num2,stim_num2,num_fish);
    %get the number of time points per stimulus
    t_perstim = size(resh_trace,2);
    % for all the fish
    for fish = 1:num_fish
        % get only the traces from this fish
        resh_fish = resh_trace(fish_ori(:,1)==fish,:,:);
        % get the number of traces for this fish
        fish_trace_num = size(resh_fish,1);
        %now reshape again for calculating the correlation across traces for
        %each stimulus
        corr_trace = reshape(resh_fish,fish_trace_num*t_perstim,stim_num2);
        %calculate and plot an allvall corr matrix
        corr_perfish(:,:,fish) = corr(corr_trace);
    end
    % calculate the average correlation
    rho = mean(corr_perfish,3);
    figure
    imagesc(rho)
    set(gca,'XTick',1:stim_num2,'XTickLabels',stim_labels,'FontSize',20,...
        'XTickLabelRotation',90)
    set(gca,'YTick',1:stim_num2,'YTickLabels',stim_labels,'FontSize',20)
    set(gca,'CLim',[-1,1])
    title(data(datas).name,'Interpreter','None')
    axis square
    colorbar
    
    % assemble the figure path
    file_path = strjoin({'corrOverall',data(datas).name,'.png'},'_');
    saveas(gcf, fullfile(fig_path,file_path), 'png')
    % save the matrix for plotting
    correlations(:,:,datas) = rho;
        
end

% also calculate a subtraction matrix
figure
imagesc(correlations(:,:,1)-correlations(:,:,2))
set(gca,'XTick',1:stim_num2,'XTickLabels',stim_labels,'FontSize',20,...
    'XTickLabelRotation',90)
set(gca,'YTick',1:stim_num2,'YTickLabels',stim_labels,'FontSize',20)
% set(gca,'CLim',[-1,1])
title(strjoin({'Delta correlation',data(1).name,data(2).name},'_'),'Interpreter','None', 'FontSize',12)
axis square
colorbar

% assemble the figure path 
file_path = strjoin({'corrOverall','delta','.png'},'_');
saveas(gcf, fullfile(fig_path,file_path), 'png')

autoArrangeFigures
%% Calculate correlation over time

close all

% %define the stim labels based on the paradigm
% %extract the actual file name
% [~,stim_name,~] = fileparts(name_cell{1});
% %scan for the p17b
% if ~isempty(strfind(stim_name,'p17b'))
%     %if it's p17b
%     stim_labels = {'Red','Green','Blue','UV'};
% else %if it's p6p8 instead
%     %define the stim labels (for p6p8 data)
%     stim_labels = {'Red CK','UV CK','Red GR','UV GR','Red FL','UV FL'};
% end

%define the sets of time regions to correlate
time_corr = (1:40)';
% time_corr = [1:10;11:20;21:30;31:40];
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

%allocate memory to store the correlations
tcorr_mat = zeros(num_data,num_times,stim_num2,stim_num2);
tcorr_mat_sem = zeros(num_data,num_times,stim_num2,stim_num2);

%for both data sets
for datas = 1:num_data
   
    %if a normal data set
    if datas <= num_data
        %load the clusters
        conc_trace = data(datas).conc_trace;
        %get the number of traces
        trace_num = size(conc_trace,1);
        %reshape the matrix to calculate correlations between stimuli
        resh_trace = reshape(conc_trace,trace_num,time_num,stim_num2);
        %extract only the stim period
        resh_trace = resh_trace(:,21:60,:);
%     else
%         %use the stimuli
%         resh_trace = cone_stim;
%         
%         %get the number of traces
%         trace_num = size(resh_trace,1);
%         %reshape the matrix 
%         resh_trace = reshape(resh_trace,trace_num,size(resh_trace,2)/stim_num2,stim_num2);
    end
    % get the trace fish of origin
    fish_ori = data(datas).fish_ori;
    % get the number of fish
    num_fish = size(unique(fish_ori(:,1)),1);

    %for all the times
    for times = 1:num_times

        % allocate memory for the output
        corr_perfish = zeros(stim_num2,stim_num2,num_fish);

        %now reshape again for calculating the correlation across traces for
        %each stimulus
        corr_trace = squeeze(mean(resh_trace(:,time_corr(times,:),:),2));
        % for all the fish
        for fish = 1:num_fish
            % get only the traces from this fish
            resh_fish = corr_trace(fish_ori(:,1)==fish,:,:);

            %calculate and plot an allvall corr matrix
            corr_perfish(:,:,fish) = corr(resh_fish);
        end
        % save the average and sem corr matrix
        tcorr_mat(datas,times,:,:) = mean(corr_perfish,3);
        tcorr_mat_sem(datas,times,:,:) = std(corr_perfish,0,3)./sqrt(num_fish);
    end


%     figure
%     %for all the times
%     for times = 1:num_times
%         subplot(round(sqrt(size(time_corr,1))),ceil(sqrt(size(time_corr,1))),times)
%         imagesc(tcorr_mat(datas,times,:,:))
%         set(gca,'XTick',1:stim_num2,'XTickLabels',stim_labels,'FontSize',20,...
%             'XTickLabelRotation',90)
%         set(gca,'YTick',1:stim_num2,'YTickLabels',stim_labels,'FontSize',20)
%         colorbar
%     end
    
    %define the possible combinations
    comb_vec = combnk(1:stim_num2,2);
    
    %get the number of combinations
    comb_num = size(comb_vec,1);
    %generate a color map for the traces
    c_map = parula(comb_num);
    %allocate memory for the legend
    legend_cell = cell(comb_num,1);
    figure
    %for all the combs
    for combs = 1:comb_num
        dat1 = squeeze(tcorr_mat(datas,:,comb_vec(combs,1),comb_vec(combs,2)));
        sem1 = squeeze(tcorr_mat_sem(datas,:,comb_vec(combs,1),comb_vec(combs,2)));
        errorbar(timep_axis,dat1,sem1,'-o','Color',c_map(combs,:),'MarkerFaceColor',c_map(combs,:),...
            'MarkerEdgeColor',c_map(combs,:))
        hold('on')
        %assemble the legend
        legend_cell{combs} = strcat(stim_labels{comb_vec(combs,1)},'_',stim_labels{comb_vec(combs,2)});
    end
    set(gca,'FontSize',20)
    xlabel('Time (s)','FontSize',20)
    ylabel('R (a.u.)','FontSize',20)
    legend(legend_cell,'Interpreter','none','Location','bestoutside','FontSize',10)
    axis square
    title(data(datas).name,'Interpreter','None')
    set(gca,'YLim',[-0.8 1])
    % assemble the figure path
    file_path = strjoin({'corrOverTime',data(datas).name,'.png'},'_');
    saveas(gcf, fullfile(fig_path,file_path), 'png')
end

autoArrangeFigures
%% Generate correlation matrices Niessing-style

close all

% %define the stim labels based on the paradigm
% %extract the actual file name
% [~,stim_name,~] = fileparts(name_cell{1});
% %scan for the p17b
% if ~isempty(strfind(stim_name,'p17b'))
%     %if it's p17b
%     stim_labels = {'Red','Green','Blue','UV'};
% else %if it's p6p8 instead
%     %define the stim labels (for p6p8 data)
%     stim_labels = {'Red CK','UV CK','Red GR','UV GR','Red FL','UV FL'};
% end

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

%allocate memory to store the correlations
tcorr_mat = zeros(num_data,num_times,stim_num2,stim_num2);

%for both data sets
for datas = 1:num_data
    %if a normal data set
    if datas <= num_data
        %load the clusters
        conc_trace = data(datas).conc_trace;
        %get the number of traces
        trace_num = size(conc_trace,1);
        %reshape the matrix to calculate correlations between stimuli
        resh_trace = reshape(conc_trace,trace_num,time_num,stim_num2);
        %extract only the stim period
        resh_trace = resh_trace(:,21:60,:);
%     else
%         %use the stimuli
%         resh_trace = cone_stim;
%         
%         %get the number of traces
%         trace_num = size(resh_trace,1);
%         %reshape the matrix 
%         resh_trace = reshape(resh_trace,trace_num,size(resh_trace,2)/stim_num2,stim_num2);
    end
    % get the trace fish of origin
    fish_ori = data(datas).fish_ori;
    % get the number of fish
    num_fish = size(unique(fish_ori(:,1)),1);
    
    %for all the times
    for times = 1:num_times
         % allocate memory for the output
        corr_perfish = zeros(stim_num2,stim_num2,num_fish);

        %now reshape again for calculating the correlation across traces for
        %each stimulus
        corr_trace = squeeze(mean(resh_trace(:,time_corr(times,:),:),2));
        % for all the fish
        for fish = 1:num_fish
            % get only the traces from this fish
            resh_fish = corr_trace(fish_ori(:,1)==fish,:,:);

            %calculate and plot an allvall corr matrix
            corr_perfish(:,:,fish) = corr(resh_fish);
        end

        %calculate and plot an allvall corr matrix
        tcorr_mat(datas,times,:,:) = mean(corr_perfish,3);
%         c_lims = [min(temp_corr(:)), max(temp_corr(:))];
    end
    temp_corr = tcorr_mat(datas,:,:,:);
%     c_lims = [min(temp_corr(:)), max(temp_corr(:))];
    c_lims = [-1, 1];
    
    figure
    %for all the plots
    for times = 1:num_times
        subplot(round(sqrt(num_times)),ceil(sqrt(num_times)),times)
        imagesc(squeeze(tcorr_mat(datas,times,:,:)));
        colormap('parula')
%         get(gca,'CLim')
        set(gca,'CLim',c_lims)
       
%         caxis(gca,c_lims)
        axis('square')
        title(strcat('Time point:',num2str(times)));
        % assemble the figure path
        file_path = strjoin({'corrOverTimeMatrix',data(datas).name,'.png'},'_');
        saveas(gcf, fullfile(fig_path,file_path), 'png')
    end
    
end


% also calculate a subtraction matrix
figure
for times =1:num_times
    subplot(round(sqrt(num_times)),ceil(sqrt(num_times)),times)
    imagesc(squeeze(tcorr_mat(1,times,:,:)-tcorr_mat(2,times,:,:)))
    set(gca,'XTick',1:stim_num2,'XTickLabels',stim_labels,'FontSize',8,...
    'XTickLabelRotation',90)
    set(gca,'YTick',1:stim_num2,'YTickLabels',stim_labels,'FontSize',8)
    set(gca,'CLim',[-0.5,0.5])
    title(strjoin({'Delta correlation',data(1).name,data(2).name},'_'),'Interpreter','None', 'FontSize',8)
    axis square
    colorbar
end

file_path = strjoin({'corrOverTimeMatrix','delta','.png'},'_');
saveas(gcf, fullfile(fig_path,file_path), 'png')
autoArrangeFigures
%% Calculate correlation in the time dimension

close all

% define the period to take
corr_period = 1:time_num;
% get the number of points
t_perstim = length(corr_period);
% allocate memory to save the matrices
correlations = zeros(t_perstim,t_perstim,stim_num2,num_data);
%for both data sets
for datas = 1:num_data
    %if a normal data set
    if datas <= num_data
        %load the clusters
        conc_trace = data(datas).conc_trace;
        %get the number of traces
        trace_num = size(conc_trace,1);
        %reshape the matrix to calculate correlations between stimuli
        resh_trace = reshape(conc_trace,trace_num,time_num,stim_num2);
        %extract only the stim period
        resh_trace = resh_trace(:,corr_period,:);
%     else
%         %use the stimuli
%         resh_trace = cone_stim;
%         %get the number of traces
%         trace_num = size(resh_trace,1)/stim_num2;
    end
    
    % get the trace fish of origin
    fish_ori = data(datas).fish_ori;
    % get the number of fish
    num_fish = size(unique(fish_ori(:,1)),1);
%     %get the number of time points per stimulus
%     t_perstim = size(resh_trace,2);
    % allocate memory for the output
    corr_perfish = zeros(t_perstim,t_perstim,num_fish);


    %now reshape again for calculating the correlation across traces for
    %each stimulus
%     corr_trace = reshape(resh_trace,trace_num*t_perstim,stim_num2);
    figure
    % for all the stimuli
    for stim = 1:stim_num2
        % get the traces for this stim
        stim_trace = resh_trace(:,:,stim);
        % for all the fish
        for fish = 1:num_fish
            % get only the traces from this fish
            corr_trace = stim_trace(fish_ori(:,1)==fish,:);
            %calculate and plot an allvall corr matrix
            corr_perfish(:,:,fish) = corr(corr_trace);
        end
%         corr_trace = reshape(resh_trace,trace_num,t_perstim*stim_num2);
        
        %calculate and plot an allvall corr matrix
        rho = mean(corr_perfish,3);
        % store the matrix for later use
        correlations(:,:,stim,datas) = rho;
        subplot(round(sqrt(stim_num2)),ceil(sqrt(stim_num2)),stim)
        imagesc(rho)
    %     set(gca,'XTick',1:stim_num2,'XTickLabels',stim_labels,'FontSize',20,...
    %         'XTickLabelRotation',90)
    %     set(gca,'YTick',1:stim_num2,'YTickLabels',stim_labels,'FontSize',20)
    %     title(data(datas).name,'Interpreter','None')
        axis square
        % assemble the figure path
        file_path = strjoin({'corrTime',data(datas).name,'.png'},'_');
        saveas(gcf, fullfile(fig_path,file_path), 'png')
    end
%     colorbar
        
end

% also calculate a subtraction matrix
figure
for stim =1:stim_num2
    subplot(round(sqrt(stim_num2)),ceil(sqrt(stim_num2)),stim)
    imagesc(squeeze(correlations(:,:,stim,1)-correlations(:,:,stim,2)))
    title(strjoin({'Delta correlation',data(1).name,data(2).name},'_'),'Interpreter','None', 'FontSize',8)
    axis square
    colorbar
end
% assemble the figure path
file_path = strjoin({'corrTime','delta','.png'},'_');
saveas(gcf, fullfile(fig_path,file_path), 'png')
autoArrangeFigures
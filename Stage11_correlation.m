clearvars
close all

load('paths.mat')
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'Correlations\');


data = load_clusters(cluster_path);
% define the color scheme depending on the stimulus type
if contains(data(1).name,'p17')
    color_scheme = [1 0 0;0 1 0;0 0 1;1 0 1];
else
    color_scheme = distinguishable_colors(6);
end
%% Get the region filtering index

% get the number of dataset
num_data = size(data,2);
% allocate memory for the index
index_cell = cell(num_data,1);
% for all the datasets
for datas = 1:num_data
    region_combination = 1;
    
    % define which regions to keep depending on the dataset
    if contains(data(datas).name, {'Syn','syn'})
        region_list = {'AF10'};
    else
        region_list = {'R-TcN','R-TcP'};
    end
    % load the anatomy info
    anatomy_info = data(datas).anatomy_info(:,1);
    
    % separate the traces by region
    [region_cell,~] = region_split(data(datas).single_reps,...
        anatomy_info,data(datas).name,region_combination,region_list);
    % rewrite the index vector
    index_cell{datas} = region_cell{3}==1;
    
end
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
%% Calculate correlation matrices for each data set

close all

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
        resh_fish = resh_trace(fish_ori(:,1)==fish&index_cell{datas}==1,:,:);
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
        'XTickLabelRotation',45)
    set(gca,'YTick',1:stim_num2,'YTickLabels',stim_labels,'FontSize',20)
    set(gca,'CLim',[-1,1])
    set(gca,'TickLength',[0 0])
    title(data(datas).figure_name,'Interpreter','None')
    axis square
    cbar = colorbar;
    set(cbar,'TickLength',0)
    
    % assemble the figure path
    file_path = strjoin({'corrOverall',data(datas).name,'.png'},'_');
%     saveas(gcf, fullfile(fig_path,file_path), 'png')
    print(fullfile(fig_path,file_path),'-dpng','-r600')
    % save the matrix for plotting
    correlations(:,:,datas) = rho;
        
end

% also calculate a subtraction matrix
figure
imagesc(correlations(:,:,1)-correlations(:,:,2))
set(gca,'XTick',1:stim_num2,'XTickLabels',stim_labels,'FontSize',20,...
    'XTickLabelRotation',45)
set(gca,'YTick',1:stim_num2,'YTickLabels',stim_labels,'FontSize',20)
% set(gca,'CLim',[-1,1])
set(gca,'TickLength',[0 0])
title('Delta correlation','Interpreter','None', 'FontSize',20)
axis square
cbar = colorbar;
set(cbar,'TickLength',0)

% assemble the figure path 
file_path = strjoin({'corrOverall','delta',data(1).name,data(2).name,'.png'},'_');
% saveas(gcf, fullfile(fig_path,file_path), 'png')
print(fullfile(fig_path,file_path),'-dpng','-r600')

autoArrangeFigures
%% Calculate correlation over time

close all

% define the marker size
marker_size = 2;
% marker shape left
shape_left_x = [-0.1 0 0];
shape_left_y = [-0.01 -0.01 0.01];
% marker shape right
shape_right_x = [0.1 0 0];
shape_right_y = [-0.01 -0.01 0.01];
%define the sets of time regions to correlate
time_corr = (1:10)';
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
            resh_fish = corr_trace(fish_ori(:,1)==fish&index_cell{datas}==1,:,:);

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
%         errorbar(timep_axis,dat1,sem1,'-o','Color',c_map(combs,:),'MarkerFaceColor',c_map(combs,:),...
%             'MarkerEdgeColor',c_map(combs,:))

        % plot split markers
        plotCustMarkMod(1:size(dat1,2),dat1,shape_left_x,shape_left_y,...
            marker_size,color_scheme(comb_vec(combs,1),:),color_scheme(comb_vec(combs,1),:))
        hold('on')
        plotCustMarkMod(1:size(dat1,2),dat1,shape_right_x,shape_right_y,...
            marker_size,color_scheme(comb_vec(combs,2),:),color_scheme(comb_vec(combs,2),:))
        %assemble the legend
        legend_cell{combs} = strcat(stim_labels{comb_vec(combs,1)},'_',stim_labels{comb_vec(combs,2)});
    end
    set(gca,'FontSize',20)
    xlabel('Time (s)','FontSize',20)
    ylabel('Correlation (a.u.)','FontSize',20)
%     legend(legend_cell,'Interpreter','none','Location','bestoutside','FontSize',10)
    title(data(datas).figure_name,'Interpreter','None')
    set(gca,'YLim',[-0.8 1])
    set(gca,'TickLength',[0 0])
    axis tight
    % assemble the figure path
    file_path = strjoin({'corrOverTime',data(datas).name,'.png'},'_');
%     saveas(gcf, fullfile(fig_path,file_path), 'png')
    print(fullfile(fig_path,file_path),'-dpng','-r600')
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
% time_corr = (1:10)';
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
            resh_fish = corr_trace(fish_ori(:,1)==fish&index_cell{datas}==1,:,:);

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
        set(gca,'TickLength',[0 0])
%         caxis(gca,c_lims)
        axis('square')
        set(gca,'XTick',1:stim_num2,'XTickLabels',stim_labels,'FontSize',8,...
            'XTickLabelRotation',45)
        set(gca,'YTick',1:stim_num2,'YTickLabels',stim_labels,'FontSize',8)
        title(strcat('Time point:',num2str(times)));
        sgtitle(strjoin({'Correlation over time',data(datas).name},'_'),'Interpreter','None', 'FontSize',8)
        colormap('hsv')
        
        
        % assemble the figure path
        file_path = strjoin({'corrOverTimeMatrix',data(datas).name,'.png'},'_');
%         saveas(gcf, fullfile(fig_path,file_path), 'png')
        print(fullfile(fig_path,file_path),'-dpng','-r600')

    end
    
end


% also calculate a subtraction matrix
figure
for times =1:num_times
    subplot(round(sqrt(num_times)),ceil(sqrt(num_times)),times)
    imagesc(squeeze(tcorr_mat(1,times,:,:)-tcorr_mat(2,times,:,:)))
    title(strcat('Time point:',num2str(times)));
    set(gca,'XTick',1:stim_num2,'XTickLabels',stim_labels,'FontSize',8,...
    'XTickLabelRotation',45)
    set(gca,'YTick',1:stim_num2,'YTickLabels',stim_labels,'FontSize',8)
    set(gca,'CLim',[-0.5,0.5])
    set(gca,'TickLength',[0 0])
    sgtitle(strjoin({'Delta correlation',data(1).name,data(2).name},'_'),'Interpreter','None', 'FontSize',8)
    axis square
    cbar = colorbar;
    set(cbar,'TickLength',0)
end

file_path = strjoin({'corrOverTimeMatrix','delta',data(1).name,data(2).name,'.png'},'_');
% saveas(gcf, fullfile(fig_path,file_path), 'png')
print(fullfile(fig_path,file_path),'-dpng','-r600')
autoArrangeFigures
%% Calculate correlation in the time dimension

close all

% define the period to take
corr_period = 10:70;
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
            corr_trace = stim_trace(fish_ori(:,1)==fish&index_cell{datas}==1,:);
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
        title(stim_labels{stim},'color',color_scheme(stim,:))
%         sgtitle(strjoin({'Time Correlation',data(datas).name},'_'),'Interpreter','None', 'FontSize',8)
        set(gca,'XTick',[],'YTick',[])
        set(gca,'TickLength',[0 0])
    %     set(gca,'XTick',1:stim_num2,'XTickLabels',stim_labels,'FontSize',20,...
    %         'XTickLabelRotation',90)
    %     set(gca,'YTick',1:stim_num2,'YTickLabels',stim_labels,'FontSize',20)
    %     title(data(datas).name,'Interpreter','None')
        axis square
        % assemble the figure path
        file_path = strjoin({'corrTime',data(datas).name,'.png'},'_');
%         saveas(gcf, fullfile(fig_path,file_path), 'png')
        print(fullfile(fig_path,file_path),'-dpng','-r600')

    end
%     colorbar
        
end

% also calculate a subtraction matrix
figure
for stim =1:stim_num2
    subplot(round(sqrt(stim_num2)),ceil(sqrt(stim_num2)),stim)
    imagesc(squeeze(correlations(:,:,stim,1)-correlations(:,:,stim,2)))
%     sgtitle(strjoin({'Delta correlation',data(1).name,data(2).name},'_'),'Interpreter','None', 'FontSize',8)
    title(stim_labels{stim},'color',color_scheme(stim,:))
    set(gca,'XTick',[],'YTick',[],'CLim',[-0.2 0.3])
    set(gca,'TickLength',[0 0])
    axis square
    cbar = colorbar;
    set(cbar,'TickLength',0)
end
% assemble the figure path
file_path = strjoin({'corrTime','delta',data(1).name,data(2).name,'.png'},'_');
% saveas(gcf, fullfile(fig_path,file_path), 'png')
print(fullfile(fig_path,file_path),'-dpng','-r600')
autoArrangeFigures
%% Plot the correlations between the gains
close all
if contains(data(1).name,'p17b')

    close all
    % allocate memory to save the correlations for later
    correlations = zeros(data(1).stim_num,data(1).stim_num,num_data);

    % for all of the datasets
    for datas = 1:num_data
        figure
        % get the gains
        delta_norm = data(datas).delta_norm;
        % correlate the gains
        gain_corr = corr(delta_norm);
        % save the correlation for the delta
        correlations(:,:,datas) = gain_corr;
        % plot it
        imagesc(gain_corr)
        title(data(datas).figure_name,'Interpreter','None')
        
        set(gca,'TickLength',[0 0])
        axis('square')
        set(gca,'XTick',1:stim_num2,'XTickLabels',stim_labels,'FontSize',15,...
            'XTickLabelRotation',45)
        set(gca,'YTick',1:stim_num2,'YTickLabels',stim_labels,'FontSize',15)
        colormap(parula)
        set(gca,'FontSize',20)
        
        cbar = colorbar;
        set(cbar,'TickLength',0)
        ylabel(cbar,'Correlation')
        % assemble the figure path
        file_path = strjoin({'gainCorr',data(datas).name,'.png'},'_');
%         saveas(gcf, fullfile(fig_path,file_path), 'png')
        print(fullfile(fig_path,file_path),'-dpng','-r600')
        
%         plotmatrix(delta_norm,strcat(colors{datas},'.'))
%         hold on
        
    end
    % also calculate a subtraction matrix
    figure
    imagesc(correlations(:,:,1)-correlations(:,:,2))
    set(gca,'XTick',1:stim_num2,'XTickLabels',stim_labels,'FontSize',20,...
        'XTickLabelRotation',45)
    set(gca,'YTick',1:stim_num2,'YTickLabels',stim_labels,'FontSize',20)
    % set(gca,'CLim',[-1,1])
    set(gca,'TickLength',[0 0])
    title(strjoin({'Delta correlation',data(1).name,data(2).name},'_'),'Interpreter','None', 'FontSize',12)
    axis square
    cbar = colorbar;
    set(cbar,'TickLength',0)
    ylabel(cbar,'Correlation')
    colormap(parula)
    
    autoArrangeFigures
    % assemble the figure path
    file_path = strjoin({'deltaGainCorr',data(1).name,data(2).name,'.png'},'_');
    saveas(gcf, fullfile(fig_path,file_path), 'png')
end
%% Calculate correlation over time for p8

% if it's not p8, skip
if contains(data(1).name,'p8')
    close all

    % define the marker size
    marker_size = 2;
    markers = {'s','o'};
    % marker shape left
    shape_left_x = {[-0.1 0 0],[-0.1 0 0 -0.1]};
    shape_left_y = {[-0.01 -0.01 0.01],[-0.01 -0.01 0.01 0.01]};
    % marker shape right
    shape_right_x = {[0.1 0 0],[0.1 0 0 0.1]};
    shape_right_y = {[-0.01 -0.01 0.01],[-0.01 -0.01 0.01 0.01]};
    %define the sets of time regions to correlate
    time_corr = (1:50)';
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
    % generate a single figure for the subplots
    figure

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
            resh_trace = resh_trace(:,11:60,:);
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
                resh_fish = corr_trace(fish_ori(:,1)==fish&index_cell{datas}==1,:,:);

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

%         %define the possible combinations
%         comb_vec = combnk(1:stim_num2,2);
        % define the combinations manually
        comb_vec = [1 2;3 4;5 6];
        colors_reduv = [1 0 0;1 0 1];
        time_axis = (1:size(tcorr_mat,2));%./data(datas).framerate;

        %get the number of combinations
        comb_num = size(comb_vec,1);

        %allocate memory for the legend
        legend_cell = cell(comb_num,1);
%         figure
        %for all the combs
        for combs = 1:comb_num
            dat1 = squeeze(tcorr_mat(datas,:,comb_vec(combs,1),comb_vec(combs,2)));
            sem1 = squeeze(tcorr_mat_sem(datas,:,comb_vec(combs,1),comb_vec(combs,2)));
    %         errorbar(timep_axis,dat1,sem1,'-o','Color',c_map(combs,:),'MarkerFaceColor',c_map(combs,:),...
    %             'MarkerEdgeColor',c_map(combs,:))
            % select the correct subplot
            subplot(3,1,combs)
            % plot split markers
%             plotCustMarkMod(time_axis,dat1,shape_left_x{datas},shape_left_y{datas},...
%                 marker_size,colors_reduv(1,:),colors_reduv(1,:))
            
            hold('on')
%             plotCustMarkMod(time_axis,dat1,shape_right_x{datas},shape_right_y{datas},...
%                 marker_size,colors_reduv(2,:),colors_reduv(2,:))
            errorbar(time_axis,dat1,sem1,strcat('k-',markers{datas}))
            %assemble the legend
            legend_cell{combs} = strcat(stim_labels{comb_vec(combs,1)},'_',stim_labels{comb_vec(combs,2)});
            set(gca,'FontSize',20)
            set(gca,'XLim',[0 size(time_axis,2)],'YLim',[-0.2 0.7])
            
            switch combs
                case 3
                    xlabel('Time (s)','FontSize',15)
%                     ylabel('Flash')
                case 2
                    set(gca,'XTick',[])
%                     ylabel('Grating')
                case 1
                    set(gca,'XTick',[])
%                     ylabel('Checker')
%                     legend('Location','bestoutside','FontSize',10)
            end
            %     ylabel('Correlation (a.u.)','FontSize',20)
%                 legend(legend_cell,'Interpreter','none','Location','bestoutside','FontSize',10)
            %         title(data(datas).figure_name,'Interpreter','None')
%             set(gca,'YLim',[-0.1 0.7])

%             pbaspect([1,2,1])
        end
        set(gca,'TickLength',[0 0])
        plot([10 10],get(gca,'YLim'),'k-')
        box off

    end
    set(gca,'FontSize',20)
%     xlabel('Time (s)','FontSize',20)
%     ylabel('Correlation (a.u.)','FontSize',20)
    %     legend(legend_cell,'Interpreter','none','Location','bestoutside','FontSize',10)
    %         title(data(datas).figure_name,'Interpreter','None')
%     set(gca,'YLim',[-0.8 1])
%     set(gca,'TickLength',[0 0])
%     axis tight
    % assemble the figure path
    file_path = strjoin({'corrRedUVTime',data(datas).name,'.png'},'_');
%     saveas(gcf, fullfile(fig_path,file_path), 'png')
    print(fullfile(fig_path,file_path),'-dpng','-r600')


end
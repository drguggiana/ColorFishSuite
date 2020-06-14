%% Clean up and load 
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

% define the dataset colors (assuming only AF10 vs Tectum)
% dataset_colors = [255 164 5;255 168 187]./255;
% dataset_colors = distinguishable_colors(2,{'w','k','r','g','b','m'});
% dataset_colors = magma(3);
% dataset_colors = dataset_colors(1:2,:);
dataset_colors = paths.afOT_colors;

% define the colormap
cmap = magma;
%% Get the region filtering index

% get the number of dataset
num_data = size(data,2);
% allocate memory for the index
index_cell = cell(num_data,1);
fig_name = cell(num_data,1);
% for all the datasets
for datas = 1:num_data
    region_combination = 1;
    
    % define which regions to keep depending on the dataset
    if contains(data(datas).name, {'Syn','syn'})
        region_list = {'AF10'};
        fig_name{datas} = 'AF10';
    else
        region_list = {'R-TcN','R-TcP'};
        fig_name{datas} = 'Tectum';
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
stim_num = data(1).stim_num;
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
correlations = zeros(stim_num,stim_num, num_data);

%for both data sets
for datas = 1:num_data
    %if a normal data set
    if datas <= num_data
        %load the clusters
        conc_trace = data(datas).conc_trace;
        %get the number of traces
        trace_num = size(conc_trace,1);
        %reshape the matrix to calculate correlations between stimuli
        resh_trace = reshape(conc_trace,trace_num,time_num,stim_num);
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
    corr_perfish = zeros(stim_num,stim_num,num_fish);
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
        corr_trace = reshape(resh_fish,fish_trace_num*t_perstim,stim_num);

        %calculate and plot an allvall corr matrix
        corr_perfish(:,:,fish) = corr(corr_trace);
    end
    % calculate the average correlation
    rho = mean(corr_perfish,3);

%     h = schemaball(rho,stim_labels,[1 0 1;1 1 0],[1 1 1]);
%     for line = 1:size(h.l,1)
%         if isnan(h.l(line))
%             continue
%         end
%         set(h.l(line),'LineWidth',4)
%     end
%     set(gca,'Color','w')
    
    % plot the correlation matrix
    h = bettercorr(rho,cmap);
    
    
    set(gca,'XTick',1:stim_num,'XTickLabels',stim_labels,'FontSize',20,...
        'XTickLabelRotation',45)
    set(gca,'YTick',1:stim_num,'YTickLabels',stim_labels,'FontSize',20)
    set(gca,'CLim',[-1,1])
    set(gca,'TickLength',[0 0],'LineWidth',2)
    title(fig_name{datas},'Interpreter','None')
    axis square
    cbar = colorbar;
    colormap(magma)
    set(cbar,'TickLength',0,'LineWidth',2)
    ylabel(cbar,'Correlation')
    set(gcf,'Color','w')

    % assemble the figure path
    file_path = strjoin({'corrOverall',data(datas).name,'.png'},'_');
    print(fullfile(fig_path,file_path),'-dpng','-r600')
    % save the matrix for plotting
    correlations(:,:,datas) = rho;
        
end

% also calculate a subtraction matrix
figure
imagesc(correlations(:,:,1)-correlations(:,:,2))
% h = bettercorr(correlations(:,:,1)-correlations(:,:,2),cmap);
set(gca,'XTick',1:stim_num,'XTickLabels',stim_labels,'FontSize',20,...
    'XTickLabelRotation',45)
set(gca,'YTick',1:stim_num,'YTickLabels',stim_labels,'FontSize',20)
% set(gca,'CLim',[-1,1])
set(gca,'TickLength',[0 0],'LineWidth',2)
title('Delta correlation','Interpreter','None', 'FontSize',20)
axis square
cbar = colorbar;
colormap(cmap)
set(cbar,'TickLength',0,'LineWidth',2)

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
tcorr_mat = zeros(num_data,num_times,stim_num,stim_num);
tcorr_mat_sem = zeros(num_data,num_times,stim_num,stim_num);

%for both data sets
for datas = 1:num_data
   
    %if a normal data set
    if datas <= num_data
        %load the clusters
        conc_trace = data(datas).conc_trace;
        %get the number of traces
        trace_num = size(conc_trace,1);
        %reshape the matrix to calculate correlations between stimuli
        resh_trace = reshape(conc_trace,trace_num,time_num,stim_num);
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
        corr_perfish = zeros(stim_num,stim_num,num_fish);

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
    comb_vec = combnk(1:stim_num,2);
    
    %get the number of combinations
    comb_num = size(comb_vec,1);
    %generate a color map for the traces
    c_map = parula(comb_num);
    %allocate memory for the legend
    legend_cell = cell(comb_num,1);
    figure
    %for all the combs
    for combs = 1:comb_num
%         % only plot the green blue correlation
%         if ~(comb_vec(combs,1)==2 && comb_vec(combs,2)==3)
%             continue
%         end
        
        
        dat1 = squeeze(tcorr_mat(datas,:,comb_vec(combs,1),comb_vec(combs,2)));
        sem1 = squeeze(tcorr_mat_sem(datas,:,comb_vec(combs,1),comb_vec(combs,2)));

        hold on
        % plot split markers
        plotCustMarkMod(1:size(dat1,2),dat1,shape_left_x,shape_left_y,...
            marker_size,color_scheme(comb_vec(combs,1),:),color_scheme(comb_vec(combs,1),:))
        hold('on')
        plotCustMarkMod(1:size(dat1,2),dat1,shape_right_x,shape_right_y,...
            marker_size,color_scheme(comb_vec(combs,2),:),color_scheme(comb_vec(combs,2),:))
        %assemble the legend
        legend_cell{combs} = strcat(stim_labels{comb_vec(combs,1)},'_',stim_labels{comb_vec(combs,2)});
    end
    set(gca,'FontSize',15)
    xlabel('Time (s)','FontSize',15)
    ylabel('Correlation (a.u.)','FontSize',15)
    box off
    legend(legend_cell,'Interpreter','none','Location','bestoutside','FontSize',10)
    title(fig_name{datas},'Interpreter','None')
    set(gca,'YLim',[-0.1 0.6])
    set(gca,'TickLength',[0 0],'LineWidth',2)
    colormap(cmap)
    axis tight
    % assemble the figure path
    file_path = strjoin({'corrOverTime',data(datas).name,'.png'},'_');
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
tcorr_mat = zeros(num_data,num_times,stim_num,stim_num);

%for both data sets
for datas = 1:num_data
    %if a normal data set
    if datas <= num_data
        %load the clusters
        conc_trace = data(datas).conc_trace;
        %get the number of traces
        trace_num = size(conc_trace,1);
        %reshape the matrix to calculate correlations between stimuli
        resh_trace = reshape(conc_trace,trace_num,time_num,stim_num);
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
        corr_perfish = zeros(stim_num,stim_num,num_fish);

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
        set(gca,'XTick',1:stim_num,'XTickLabels',stim_labels,'FontSize',8,...
            'XTickLabelRotation',45)
        set(gca,'YTick',1:stim_num,'YTickLabels',stim_labels,'FontSize',8)
        title(strcat('Time point:',num2str(times)));
        sgtitle(strjoin({'Correlation over time',data(datas).name},'_'),'Interpreter','None', 'FontSize',8)
        colormap(cmap)
        
        
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
    set(gca,'XTick',1:stim_num,'XTickLabels',stim_labels,'FontSize',8,...
    'XTickLabelRotation',45)
    set(gca,'YTick',1:stim_num,'YTickLabels',stim_labels,'FontSize',8)
    set(gca,'CLim',[-0.5,0.5])
    set(gca,'TickLength',[0 0])
    sgtitle(strjoin({'Delta correlation',data(1).name,data(2).name},'_'),'Interpreter','None', 'FontSize',8)
    axis square
    cbar = colorbar;
    set(cbar,'TickLength',0)
    colormap(cmap)
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
correlations = zeros(t_perstim,t_perstim,stim_num,num_data);
%for both data sets
for datas = 1:num_data
    %if a normal data set
    if datas <= num_data
        %load the clusters
        conc_trace = data(datas).conc_trace;
        %get the number of traces
        trace_num = size(conc_trace,1);
        %reshape the matrix to calculate correlations between stimuli
        resh_trace = reshape(conc_trace,trace_num,time_num,stim_num);
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
    
    % for all the stimuli
    for stim = 1:stim_num
        figure
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
%         subplot(round(sqrt(stim_num2)),ceil(sqrt(stim_num2)),stim)
        imagesc(rho)
%         plot(diag(rho(1:end-2,1:end-2)))
%         plot(std(rho,0,1))
        title(stim_labels{stim},'color',color_scheme(stim,:))
%         sgtitle(strjoin({'Time Correlation',data(datas).name},'_'),'Interpreter','None', 'FontSize',8)
        set(gca,'XTick',[],'YTick',[])
        set(gca,'TickLength',[0 0],'LineWidth',2)
    %     set(gca,'XTick',1:stim_num2,'XTickLabels',stim_labels,'FontSize',20,...
    %         'XTickLabelRotation',90)
    %     set(gca,'YTick',1:stim_num2,'YTickLabels',stim_labels,'FontSize',20)
    %     title(data(datas).name,'Interpreter','None')
        axis square
        colormap(cmap)
        % assemble the figure path
        file_path = strjoin({'corrTime',data(datas).name,stim_labels{stim},'.png'},'_');
%         saveas(gcf, fullfile(fig_path,file_path), 'png')
        print(fullfile(fig_path,file_path),'-dpng','-r600')

    end
%     colorbar
        
end

%% Also calculate a subtraction matrix
close all
traces = figure;
delta = figure;

% allocate memory to store the cumulative deltas
delta_matrix = zeros(stim_num,1);
% define the offset
offset = 0.5;
for stim =1:stim_num

    % get the data
    matrix1 = squeeze(correlations(:,:,stim,1));
    matrix2 = squeeze(correlations(:,:,stim,2));
    
    % blank the diagonal
    matrix1(eye(size(matrix1))==1) = NaN;
    matrix2(eye(size(matrix1))==1) = NaN;
    
    mean1 = nanmean(matrix1,1);
    mean2 = nanmean(matrix2,1);
    
    std1 = nanstd(matrix1,0,1)./sqrt(size(matrix1,1)-1);
    std2 = nanstd(matrix2,0,1)./sqrt(size(matrix1,1)-1);
    
    % get the x_range
    x_range = ((1:size(matrix1,1))-10)./data(datas).framerate;
    
    % select the target figure
    figure(traces)
    tectum_handles = shadedErrorBar(x_range,mean1+(stim_num-(stim)*offset),std1,{'-','Color',color_scheme(stim,:)});
    hold on
    af10_handles = shadedErrorBar(x_range,mean2+(stim_num-(stim)*offset),std2,{'--','Color',color_scheme(stim,:)});
    % color the line black
    set(tectum_handles.mainLine,'Color',[0 0 0],'LineWidth',1);
    set(af10_handles.mainLine,'Color',[0 0 0],'LineWidth',1);
    
    % also plot the cumulative delta
    figure(delta)
    
    % also get the cumulative differences
    delta_values = cumsum(abs(mean1-mean2));
    plot(x_range,delta_values+(stim_num-(stim)*offset),...
        'Color',color_scheme(stim,:),'LineWidth',1)
    hold on


end

figure(traces)

set(gca,'TickLength',[0 0],'LineWidth',2,'FontSize',15)
set(gca,'YTick',2:0.5:stim_num,'YTickLabels',fliplr(stim_labels))
xlabel('Time (s)')
ylabel('Relative Correlation (a.u.)')
box off
axis tight
plot([0 0],get(gca,'YLim'),'k--','LineWidth',1)
% insert the legend
legend([tectum_handles.mainLine,af10_handles.mainLine],fig_name,...
    'Location','northoutside','Orientation','horizontal')

set(gcf,'Color','w')
% assemble the figure path
file_path = fullfile(fig_path,strjoin({'corrTime','delta',...
    data(1).name,data(2).name,'.png'},'_'));
export_fig(file_path,'-r600')

% also plot the cumulative delta
figure(delta)
set(gca,'TickLength',[0 0],'LineWidth',2,'FontSize',15)
xlabel('Time (s)')
ylabel('Cumulative Delta Correlation (a.u.)')
box off
axis tight
plot([0 0],get(gca,'YLim'),'k--','LineWidth',1)

legend(stim_labels,'Location','northwest')
set(gcf,'Color','w')
% assemble the figure path
file_path = fullfile(fig_path,strjoin({'cumDeltaCorrTime','delta',...
    data(1).name,data(2).name,'.png'},'_'));
export_fig(file_path,'-r600')

autoArrangeFigures
%% Plot the correlations between the gains
close all
if contains(data(1).name,'p17b')

    close all
    % allocate memory to save the correlations for later
    gain_correlations = zeros(data(1).stim_num,data(1).stim_num,num_data);

    % for all of the datasets
    for datas = 1:num_data
        figure
        % get the gains
        delta_norm = data(datas).delta_norm;
        % correlate the gains
        gain_corr = corr(delta_norm);
        % save the correlation for the delta
        gain_correlations(:,:,datas) = gain_corr;
        % plot it
        imagesc(gain_corr)
        title(data(datas).figure_name,'Interpreter','None')
        
        set(gca,'TickLength',[0 0])
        axis('square')
        set(gca,'XTick',1:stim_num,'XTickLabels',stim_labels,'FontSize',15,...
            'XTickLabelRotation',45)
        set(gca,'YTick',1:stim_num,'YTickLabels',stim_labels,'FontSize',15)
        colormap(magma)
        set(gca,'FontSize',20)
        
        cbar = colorbar;
        set(cbar,'TickLength',0,'LineWidth',2)
        ylabel(cbar,'Correlation')
        % assemble the figure path
        file_path = strjoin({'gainCorr',data(datas).name,'.png'},'_');
        print(fullfile(fig_path,file_path),'-dpng','-r600')

    end
    
    
    % also plot the point by point correlation matrix
    % assemble the matrix of gains
    figure
    gain_matrix = cat(1,data([2 1]).delta_norm);
    % assemble the group labels
%     gain_labels = ones(size(gain_matrix,1),1);
%     gain_labels(1:size(data(1).delta_norm,1)) = 0;
    gain_labels = cat(1,repmat({'RGCs'},size(data(2).delta_norm,1),1),...
        repmat({'RAs'},size(data(1).delta_norm,1),1));
    [h,ax,bigax] = gplotmatrix(gain_matrix,[],gain_labels,dataset_colors,...
        'o',2,true,'stairs',{'Red','Green','Blue','UV'});
    % change the font size in the plots
    for plots = 1:20
        set(ax(plots),'FontSize',13,'LineWidth',2,'TickLength',[0 0])
    end
  
    % assemble the figure path
    file_path = strjoin({'gainCorrMatrix',data(datas).name,'.png'},'_');
    print(fullfile(fig_path,file_path),'-dpng','-r600')
    
    % also calculate a subtraction matrix
    figure
    imagesc(gain_correlations(:,:,1)-gain_correlations(:,:,2))
    set(gca,'XTick',1:stim_num,'XTickLabels',stim_labels,'FontSize',20,...
        'XTickLabelRotation',45)
    set(gca,'YTick',1:stim_num,'YTickLabels',stim_labels,'FontSize',20)
    % set(gca,'CLim',[-1,1])
    set(gca,'TickLength',[0 0])
    title('Delta Correlation','Interpreter','None', 'FontSize',15)
    axis square
    cbar = colorbar;
    set(cbar,'TickLength',0)
    ylabel(cbar,'Correlation')
    colormap(magma)
    
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
    tcorr_mat = zeros(num_data,num_times,stim_num,stim_num);
    tcorr_mat_sem = zeros(num_data,num_times,stim_num,stim_num);
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
            resh_trace = reshape(conc_trace,trace_num,time_num,stim_num);
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
            corr_perfish = zeros(stim_num,stim_num,num_fish);

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
            set(gca,'FontSize',15)
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
    set(gca,'FontSize',15)
    set(gcf,'Color','w')
%     xlabel('Time (s)','FontSize',20)
%     ylabel('Correlation (a.u.)','FontSize',20)
    %     legend(legend_cell,'Interpreter','none','Location','bestoutside','FontSize',10)
    %         title(data(datas).figure_name,'Interpreter','None')
%     set(gca,'YLim',[-0.8 1])
%     set(gca,'TickLength',[0 0])
%     axis tight
    % assemble the figure path
    file_path = fullfile(fig_path,strjoin({'corrRedUVTime',data(datas).name,'.png'},'_'));
%     saveas(gcf, fullfile(fig_path,file_path), 'png')
    export_fig(file_path,'-r600')


end
%% Correlation between clusters

if num_data > 1 && contains(data(1).name,'p17b')
    close all
    figure
    
    % define the colors to compare
    color_sets = [2,3;1,4];
    % get the number of sets
    set_number = size(color_sets,1);
    % allocate memory for the correlations
    rho_cell = cell(set_number,1);
    % define the number of shuffles
    num_shuffles = 500;
    % allocate memory for the clusters
    clu_ave = cell(2,set_number);
    % allocate memory for the legend handles
    legend_cell = cell(num_data,1);
    
    % load the clusters and get the cluster averages
    for datas = 1:2
        % for all the combinations
        for sets = 1:set_number
            % load the data and clusters
            conc_trace = data(datas).conc_trace(index_cell{datas}==1,:);
            conc_trace = reshape(conc_trace,[],data(datas).time_num,data(datas).stim_num);
            conc_trace = reshape(conc_trace(:,21:60,color_sets(sets,:)),size(conc_trace,1),[]);
            idx = data(datas).idx_clu(index_cell{datas}==1);
            clu_num = data(datas).clu_num;
            % allocate memory for the averages
            temp_ave = zeros(clu_num,size(conc_trace,2));
            for clu = 1:clu_num
                temp_ave(clu,:) = mean(conc_trace(idx==clu,:),1);
            end
            % store the average
            clu_ave{datas,sets} = temp_ave;
%             clu_ave{datas,sets} = conc_trace;

        end

        % calculate the correlation
%         [rho,pval] = corr(sort_traces(clu_ave{1})',sort_traces(clu_ave{2})');
        [rho1,pval1] = corr(clu_ave{datas,1}(:,1:40)',clu_ave{datas,1}(:,41:end)');
        [rho2,pval2] = corr(clu_ave{datas,2}(:,1:40)',clu_ave{datas,2}(:,41:end)');
%         rho = rho1./rho2;
        
%         rho(pval1>0.05) = NaN;
%         rho = tril(rho,1);
        
        rho_cell{datas} = rho1;

%         [N,edges] = histcounts(rho,20,'Normalization','cdf');
%         plot(N)
%         histogram(rho(rho~=0&~isnan(rho)),20,'Normalization','probability','LineWidth',2,'FaceColor',dataset_colors(3-datas,:));
%         scatter(rho1(:),rho2(:),[],dataset_colors(3-datas,:),'o','filled')
        histogram(rho1(:),20,'FaceColor',dataset_colors(3-datas,:))
        hold on
%         legend_cell{datas} = scatter(rho1(pval1<0.05),rho2(pval1<0.05),[],dataset_colors(3-datas,:),'o','filled');
        set(gca,'TickLength',[0 0],'Fontsize',20,'LineWidth',2)
%         xlabel('Correlation (a.u.)')
        ylabel('Nr clusters')
%         xlabel('Green-Blue Correlation')
%         ylabel('Red-UV Correlation')
    end
    box off
    axis square
    legend({'Tectum','AF10'},'location','northwest')
    xlabel('G-B Cluster Correlation')
    set(gcf,'Color','w')
    % assemble the figure path
    file_path = fullfile(fig_path,strjoin({'corrClusters',data(1).name,data(2).name,'.png'},'_'));
    export_fig(file_path,'-r600')
    
    % calculate the delta rho
    delta_rho = nanmean(abs(rho_cell{1}(:)))-nanmean(abs(rho_cell{2}(:)));
    
    % allocate memory for the shuffles
    shuffle_matrix = zeros(num_shuffles,1);
    % get the matrix of clusters
    clu_matrix = vertcat(clu_ave{:,1});
%     clu_num_ori = vertcat(data.clu_num);
    clu_num = size(clu_matrix,1);
    % calculate a confidence interval for the delta
    for shuff = 1:num_shuffles
        % shuffle the matrix
        shuff_idx = randperm(clu_num);
        mat1 = clu_matrix(shuff_idx(1:data(1).clu_num),:);
        mat2 = clu_matrix(shuff_idx(data(1).clu_num+1:end),:);
        
        % get the correlation matrix
        [rho1,pval1] = corr(mat1(:,1:40)',mat1(:,41:end)');
%         rho1(pval1>0.05) = NaN;
        [rho2,pval2] = corr(mat2(:,1:40)',mat2(:,41:end)');
%         rho2(pval2>0.05) = NaN;
        
        % calculate the delta and store
        shuffle_matrix(shuff) = nanmean(abs(rho1(:)))-nanmean(abs(rho2(:)));
    end
    
    figure
    histogram(shuffle_matrix)
    hold on
    plot([delta_rho delta_rho],get(gca,'YLim'),'r')
    % get the std of the shuffles
%     prctile(shuffle_matrix,95)<delta_rho
end
%% Calculate matrix of histograms

close all

% get the stim_number
stim_num = data(1).stim_num;

% create the figure and get the handle
[rho_cell,combo_list,combo_number] = matrix_subplot(data,index_cell);

% plot the results
h = figure;
% for all the combinations
for combo = 1:combo_number
    % get the linear index for plotting
%     ind = sub2ind([stim_num,stim_num],combo_list(combo,1),combo_list(combo,2));
%     subplot(stim_num,stim_num,ind)
    figure
%     subplot(3,2,combo)
    
    histogram(rho_cell{combo_list(combo,1),combo_list(combo,2),1},20,...
        'Normalization','Probability','FaceColor',dataset_colors(2,:),'LineWidth',1)
    hold on
    histogram(rho_cell{combo_list(combo,1),combo_list(combo,2),2},20,...
        'Normalization','Probability','FaceColor',dataset_colors(1,:),'LineWidth',1)
    
%     xlabel(stim_labels{combo_list(combo,1)})
%     ylabel(stim_labels{combo_list(combo,2)})
    %         set(gca,'XLim',[0 1],'YLim',[0 0.6])
    
    % format the axes
%     if combo_list(combo,2) < stim_num
%         set(gca,'XTick',[])
%         xlabel('')
%     end
%     if combo_list(combo,1) > 1
%         set(gca,'YTick',[])
%         ylabel('')
%     end

    if mod(combo,2) == 0
        ylabel('')
        set(gca,'YTick',[])
    end
    if combo < 5
        xlabel('')
        set(gca,'XTick',[])
    end
    if combo == 2
        legend({'Tectum','AF10'},'Location','northeast')
    end
    box off
    set(gca,'LineWidth',2)
    set(gca,'XLim',[-1 1],'YLim',[0 0.2])
    set(gca,'FontSize',18,'TickLength',[0 0])
    pbaspect([1,1.2,1])
    set(gcf,'Color','w')
    % assemble the figure path
    file_path = fullfile(fig_path,strjoin({'corrMatrix',data(1).name,data(2).name,num2str(combo),'.png'},'_'));
    export_fig(file_path,'-r600')

end
% set(gcf,'Color','w')
% % assemble the figure path
% file_path = fullfile(fig_path,strjoin({'corrMatrix',data(1).name,data(2).name,'.png'},'_'));
% export_fig(file_path,'-r600')
%% Calculate correlation over time as subsets

close all

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
tcorr_mat = zeros(num_data,num_times,stim_num,stim_num);
tcorr_mat_sem = zeros(num_data,num_times,stim_num,stim_num);

% create a common figure
both = figure;
% allocate memory for the data to put in the plot
both_cell = cell(num_data,1);
%for both data sets
for datas = 1:num_data
   
    %if a normal data set
    if datas <= num_data
        %load the clusters
        conc_trace = data(datas).conc_trace;
        %get the number of traces
        trace_num = size(conc_trace,1);
        %reshape the matrix to calculate correlations between stimuli
        resh_trace = reshape(conc_trace,trace_num,time_num,stim_num);
        %extract only the stim period
        resh_trace = resh_trace(:,21:60,:);
    end
    % get the trace fish of origin
    fish_ori = data(datas).fish_ori;
    % get the number of fish
    num_fish = size(unique(fish_ori(:,1)),1);

    %for all the times
    for times = 1:num_times

        % allocate memory for the output
        corr_perfish = zeros(stim_num,stim_num,num_fish);

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
    %define the possible combinations
    comb_vec = combnk(1:stim_num,2);
    
    %get the number of combinations
    comb_num = size(comb_vec,1);
    %generate a color map for the traces
    c_map = parula(comb_num);
%     %allocate memory for the legend
%     legend_cell = cell(comb_num,1);
    % allocate memory to store the lines
    handle_cell = cell(comb_num,1);
    
    %for all the combs
    for combs = 1:comb_num
%         % only plot the green blue correlation
%         if ~(comb_vec(combs,1)==2 && comb_vec(combs,2)==3)
%             continue
%         end
        
        
        dat1 = squeeze(tcorr_mat(datas,:,comb_vec(combs,1),comb_vec(combs,2)));
        sem1 = squeeze(tcorr_mat_sem(datas,:,comb_vec(combs,1),comb_vec(combs,2)));
        handle_cell{combs} = shadedErrorBar(timep_axis,dat1,sem1,...
            {'Color',dataset_colors(datas,:),'LineWidth',2});
        hold on
%         legend_cell{combs} = strcat(stim_labels{comb_vec(combs,1)},'_',stim_labels{comb_vec(combs,2)});
       
    end
    % store the handle cell
    both_cell{datas} = handle_cell;
end

% select only the lines desired for plotting
combs_target = [2 3;2 3];
% allocate memory for the legends
combs_legend = cell(size(combs_target,1),2);
% for both data sets
for datas = 1:num_data
    % for all the combs
    for combs = 1:comb_num
        % only plot the green blue correlation
        if ~(comb_vec(combs,1)==combs_target(datas,1) && ...
                comb_vec(combs,2)==combs_target(datas,2))
            delete(both_cell{datas}{combs}.mainLine)
            delete(both_cell{datas}{combs}.patch)
            delete(both_cell{datas}{combs}.edge)
            continue
        else
            % capture the handle to the mainLine
            combs_legend{datas,1} = both_cell{datas}{combs}.mainLine;
            combs_legend{datas,2} = strjoin({fig_name{datas},...
                stim_labels{comb_vec(combs,1)},'-',stim_labels{comb_vec(combs,2)}},' ');
        end
    end
end

set(gca,'FontSize',15)
xlabel('Time (s)','FontSize',15)
ylabel('Correlation (a.u.)','FontSize',15)
box off
legend(vertcat(combs_legend{:,1}),combs_legend(:,2),'location','northoutside','orientation','horizontal')
%     legend(legend_cell,'Interpreter','none','Location','bestoutside','FontSize',10)
%     title(fig_name{datas},'Interpreter','None')
set(gca,'YLim',[-0.1 0.6])
set(gca,'TickLength',[0 0],'LineWidth',2)
%     axis tight
% assemble the figure path
file_path = strjoin({'corrOverTimeSub',data(datas).name,'.png'},'_');
print(fullfile(fig_path,file_path),'-dpng','-r600')

autoArrangeFigures
%% Correlate the responses across fish

close all
% define the font size
fontsize = 20;
%for both data sets
for datas = 1:num_data
    figure
    %load the clusters
    conc_trace = data(datas).conc_trace;
    
    % get the trace fish of origin
    fish_ori = data(datas).fish_ori;
    % get the number of fish
    num_fish = size(unique(fish_ori(:,1)),1);
    % get the number of clusters
    clu_num = data(datas).clu_num;
    % allocate memory for the averages
    fish_average = zeros(clu_num,size(conc_trace,2),num_fish);
    fish_cluster_traces = zeros(clu_num,num_fish);
    % calculate the per-fish cluster averages
    for fish = 1:num_fish
        % get the idx_clu for this fish
        idx_fish = data(datas).idx_clu(fish_ori(:,1)==fish);
        % get the traces for this fish
        fish_trace = conc_trace(fish_ori(:,1)==fish,:);
        % for all the clusters
        for clu = 1:clu_num
            fish_average(clu,:,fish) = nanmean(fish_trace(idx_fish==clu,:),1);
            fish_cluster_traces(clu,fish) = sum(idx_fish==clu);
        end
    end
    % get all the combinations of fish
    fish_comb = nchoosek(1:num_fish,2);
    number_combos = length(fish_comb);
    % allocate memory for a correlation matrix
    fish_correlation = zeros(num_fish);
    % for all the fish combos
    for combo = 1:number_combos
%         % load the 2 fish in question
%         fish1 = conc_trace(fish_ori(:,1)==fish_comb(combo,1),:);
%         fish2 = conc_trace(fish_ori(:,1)==fish_comb(combo,2),:);
        % get the cluster averages
        fish1 = fish_average(:,:,fish_comb(combo,1));
        fish2 = fish_average(:,:,fish_comb(combo,2));
        % get the corelation matrix
        corr_matrix = corr(fish1',fish2');
        % take the off diagonal entries (since clusters are matched)
%         corr_matrix = tril(corr_matrix,1);
%         corr_matrix = corr_matrix(corr_matrix~=0);
        corr_matrix = diag(corr_matrix);
        % average
        fish_correlation(fish_comb(combo,1),fish_comb(combo,2)) = ...
            nanmean(corr_matrix(:));
        
    end
    % plot the matrix
    imagesc(fish_correlation)
    set(gca,'TickLength',[0 0],'LineWidth',2)
    set(gca,'CLim',[0 1])
    xlabel('Fish')
    ylabel('Fish')
    
    cba = colorbar;
    ylabel(cba,'Correlation')
    set(cba,'TickLength',0,'LineWidth',2)
    colormap(magma)
    set(gca,'FontSize',fontsize)
    axis square
    axis tight
    file_path = strjoin({'fishCorrelation',data(datas).name,'.png'},'_');
    print(fullfile(fig_path,file_path),'-dpng','-r600')
    
    % Plot the number of traces of every cluster in every fish
    figure
    imagesc(normr_1(fish_cluster_traces,2))
    set(gca,'TickLength',[0 0],'LineWidth',2)
    xlabel('Fish')
    ylabel('Cluster')
    
    cba = colorbar;
    ylabel(cba,'Proportion per fish')
    set(cba,'TickLength',0,'LineWidth',2)
    colormap(magma)
    set(gca,'FontSize',fontsize)
    axis square
    axis tight
    file_path = strjoin({'fishClusters',data(datas).name,'.png'},'_');
    print(fullfile(fig_path,file_path),'-dpng','-r600')
    
end
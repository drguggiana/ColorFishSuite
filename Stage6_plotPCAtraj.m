%% Clean up and load data

clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'PCA\');

data = load_clusters(cluster_path);
%% Calculate PCA trajectories, Niessing-style

close all
% define the stimulus set to use
stim_set = 1;

%define the stim labels based on the paradigm
%extract the actual file name
% [~,stim_name,~] = fileparts(name_cell{1});
%scan for the p17b
if contains(data(1).name,'p17b')
    %if it's p17b
    stim_labels = {'Red','Green','Blue','UV'};
    %define the plot colors
    plot_col = [1 0 0;0 1 0;0 0 1;1 0 1];
    plot_marker = {'o', 'o','o', 'o',};
    var_threshold = 0.6;
else %if it's p6p8 instead
    %define the stim labels (for p6p8 data)
    stim_labels = {'Red CK','UV CK','Red GR','UV GR','Red FL','UV FL'};
    %define the plot colors
%     plot_col = [1 0 0;0 0 1;0.8 0 0;0 0 0.8;0.4 0 0;0 0 0.4];
    plot_col = [1 0 0;1 0 1;1 0 0;1 0 1;1 0 0;1 0 1];
    plot_marker = {'o', 'o', 's', 's', 'd', 'd'};
    var_threshold = 0.9;
end

% define the region set to use
region_set = 1;
% define which regions to include in the analysis
switch region_set
    case 1
        tectum_regions = {'R-TcN','R-TcP'};
        af_regions = {'AF10'};

    case 2
        tectum_regions = {'R-TcN','R-TcP','R-Cb','R-Hb','R-Pt','L-TcN','L-TcP','L-Cb','L-Hb','L-Pt'};
        af_regions = {'AF4','AF5','AF8','AF9','AF10'};
        
end


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
    
    % define the regions to be considered (depending on the stimulus protocol)
    if contains(data(datas).name,{'syn','Syn'})
        region_list = af_regions;
        if region_set == 1
            fig_name = 'AF10';
        else
            fig_name = data(datas).figure_name;
        end
    else
        region_list = tectum_regions;
        if region_set == 1
            fig_name = 'Tectum';
        else
            fig_name = data(datas).figure_name;
        end
    end
    
    % get the time and stim nums
    time_num = data(datas).time_num;
    stim_num = data(datas).stim_num;
    
    % get the stimuli to use
    stim_vector = pick_stim_vector(stim_set,stim_num,data(datas).name);

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
    else
        %use the stimuli
        resh_trace = cone_stim;
        
        %get the number of traces
        trace_num = size(resh_trace,1);
        %reshape the matrix 
        resh_trace = reshape(resh_trace,trace_num,size(resh_trace,2)/stim_num,stim_num);
    end
    
    % separate the traces by region
    if isempty(data(datas).anatomy_info)
        anatomy_info = [];
    else
        anatomy_info = data(datas).anatomy_info(:,1);
    end
    [region_data,num_regions] = region_split(resh_trace,anatomy_info,data(datas).name,1,region_list);

    
    % run the analysis separately for every region
    for region = 1:num_regions
        % redefine resh_trace for the region
        resh_trace = region_data{region,1};
        
        % if there are less than 3 dimensions, skip the plot
        if size(resh_trace,1) < 3
            continue
        end
        % run the PCA    
        pca_mat = zeros(size(resh_trace,2),3,stim_num);
        figure
        %for all the times
        for stim = stim_vector


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
%             pca_mat(:,:,stim) = pca_res(:,1:3);

            pca_mat(:,:,stim) = temp_mat*pca_res(:,1:3);

            % plot the results
            options = struct([]);
            options(1).line=1;
            options(1).threeD=0;
            options(1).marker=plot_marker{stim};
            plot_trajectory(pca_mat(:,:,stim), plot_col(stim,:),options)
%             %for all the points
%             for points = 1:size(pca_mat,1)
% %                 plot3(pca_mat(points,1,stim),pca_mat(points,2,stim),pca_mat(points,3,stim),...
% %                     'Marker','o','MarkerSize',points/5,'MarkerFaceColor',plot_col(stim,:),...
% %                     'MarkerEdgeColor',plot_col(stim,:))
%                 plot(pca_mat(points,1,stim),pca_mat(points,2,stim),...
%                     'Marker','o','MarkerSize',points/5,'MarkerFaceColor',plot_col(stim,:),...
%                     'MarkerEdgeColor',plot_col(stim,:))
%                 hold('on')
%             end
%             %also plot the lines
% %             plot3(pca_mat(:,1,stim),pca_mat(:,2,stim),pca_mat(:,3,stim),...
% %                 'Color',plot_col(stim,:))
%             plot(pca_mat(:,1,stim),pca_mat(:,2,stim),...
%                             'Color',plot_col(stim,:))
        end
        %get the file name
    %     temp_name = strsplit(data{datas},'\');
    %     temp_name = strsplit(temp_name{end},'_');
        sgtitle(strcat(data(datas).name,'_',region_data{region,2}),'FontSize',20,'Interpreter','None')
%         xlabel('PC 1','FontSize',20)
%         ylabel('PC 2','FontSize',20)
%         zlabel('PC 3','FontSize',20)
        
        % assemble the figure path 
%         file_path = strjoin({'trajectory',data(datas).name,...
%           'region',region_data{region,2},'set',num2str(stim_set),'.png'},'_');
%         saveas(gcf, fullfile(fig_path,file_path), 'png')
    
    end
end
autoArrangeFigures
%% Align pca trajectories from different fish with CCA and then plot
close all
%define the minimal dimension threshold for CCA
min_dim = 3;% 3; for UV Red stuff
% define the variance threshold
% allocate memory to store the aligned PCs
cca_cell = cell(num_data,stim_num);
% define the fontsize
fontsize = 15;
%for all the fish
for datas = 1:num_data
    
    %show the current fish
    fprintf(strcat('Current dataset:',num2str(datas),'\r\n'))
    
        % define the regions to be considered (depending on the stimulus protocol)
    if contains(data(datas).name,{'syn','Syn'})
        region_list = af_regions;
        if region_set == 1
            fig_name = 'AF10';
        else
            fig_name = data(datas).figure_name;
        end
    else
        region_list = tectum_regions;
        if region_set == 1
            fig_name = 'Tectum';
        else
            fig_name = data(datas).figure_name;
        end
    end
        
    % get the stimuli to use
    stim_vector = pick_stim_vector(stim_set,stim_num,data(datas).name);

    raw_trace = data(datas).conc_trace;
    
    resh_trace = reshape(raw_trace,size(raw_trace,1),time_num,stim_num);
    %clip the edges away
    resh_trace = resh_trace(:,21:60,:);
    
    % get the number of animals and the animal info
    fish_ori_all = data(datas).fish_ori;
    num_animals = length(unique(fish_ori_all(:,1)));
    
    % exclude the anatomy if it's not present
    if isempty(data(datas).anatomy_info)
        anatomy_info = [];
    else
        anatomy_info = data(datas).anatomy_info(:,1);
    end
    [region_data,num_regions] = region_split(resh_trace,anatomy_info,data(datas).name,1,region_list);
    % for all the regions
    for region = 1:num_regions
        % get the region indexer
        region_idx = region_data{region,3};
        
        pca_mat = zeros(size(resh_trace,2),3,stim_num);
        % allocate memory for all animals
        all_animals_mat = zeros(40,3,length(stim_vector));
        % initialize the skipped fish counter
        skip_count = 1;
        figure
        %for all the times
        for stim = stim_vector
            % allocate a structure to store the pca data
            pca_struct = struct([]);
                    % split between animals
            % for all the animals
            for animals = 1:num_animals
        %         %get the number of time points per stimulus
        %         t_perstim = size(time_corr,2);
        %         %now reshape again for calculating the correlation across traces for
        %         %each stimulus
        %         corr_trace = reshape(resh_trace(:,time_corr(times,:),:),trace_num*t_perstim,stim_num2);
        %         temp_mat = squeeze(mean(resh_trace(:,time_corr(times,:),:),2));
                temp_mat = resh_trace(fish_ori_all(:,1)==animals&region_idx==1,:,stim)';
                %calculate and plot an allvall corr matrix
        %         temp_mat = corr(corr_trace);
        %         c_lims = [min(temp_corr(:)), max(temp_corr(:))];
                [pca_struct(animals).coeff,pca_struct(animals).score,pca_struct(animals).latent] = ...
                pca(temp_mat);
            
%                 % JUST FOR TESTING, FIX LATER
%                 pca_struct(animals).score = temp_mat*pca_struct(animals).coeff;
            end

            % run the CCA and align the spaces
            % get the first animal's pca results
            var_1 = pca_struct(1).latent;
            traj_1  = pca_struct(1).score;
            
            %determine the 60% variance threshold
            dim_thres = find(cumsum(var_1./sum(var_1))>var_threshold,1,'first');
%             dim_thres = 14;
            %if there are too few dimensions left, skip
            if isempty(dim_thres) || dim_thres < min_dim
                fprintf(strcat('Skipped dataset',num2str(skip_count),'\r\n'))
                skip_count = skip_count + 1;
                continue
            end
            %apply the threshold
            traj_1 = traj_1(:,1:dim_thres);

            % for all the animals minus the first one
            for animals = 2:num_animals

                %get the corresponding pca results
                var_2 = pca_struct(animals).latent;
                traj_2  = pca_struct(animals).score;

                try
                    traj_2 = traj_2(:,1:dim_thres);
                catch
                    fprintf(strcat('Skipped dataset',num2str(skip_count),'\r\n'))
                    skip_count = skip_count + 1;
                    continue
                end
                % align the current animal to the first one
                [pca_struct(1).A,pca_struct(animals).B,...
                pca_struct(animals).r,~,~,stat] = canoncorr(traj_1,traj_2);
            end

            % plot the aligned trajectories
            options = struct([]);
            options(1).line=1;
            options(1).threeD=1;
            options(1).marker=plot_marker{stim};
            % plot the first animal
%             plot_trajectory(pca_struct(1).score(:,1:3),plot_col(stim,:),options)
            pca_struct(1).new_space = pca_struct(1).score(:,1:dim_thres);
            % for all the animals
            for animals = 2:num_animals
                if isempty(pca_struct(animals).B) || (size(pca_struct(1).A,1)>size(pca_struct(animals).B,1))
                    continue
                end
                % transform the space
                pca_struct(animals).new_space = pca_struct(animals).score(:,1:dim_thres)*pca_struct(animals).B/...
                    pca_struct(1).A;
                
                
%                 plot_trajectory(pca_struct(animals).new_space(:,1:3),plot_col(stim,:),options)
            end

            % calculate the average trajectory in the common space
            all_animals = mean(cat(3,pca_struct.new_space),3);
            % store for later plotting
            all_animals_mat(:,:,stim) = all_animals(:,1:3);
            
            plot_trajectory(all_animals(:,1:3),plot_col(stim,:),options)

            view(2)
%             view([90 0])
            set(gca,'TickLength',[0 0],'LineWidth',2,'FontSize',fontsize)
            % store the cca data (only for 1 region
            cca_cell{datas,stim} = pca_struct;

        end
        
        
%         figure
%         % plot an animation
%         for i = 1:size(all_animals,1)
%             for stim = 1:4
%                 plot_trajectory(all_animals_mat(1:i,:,stim),plot_col(stim,:),options)
%                 view(2)
%                 hold on
%             end
%             pause(0.3)
%         end
        

        set(gcf,'Color','w')
        axis equal
        sgtitle(fig_name,'FontSize',fontsize,'Interpreter','None')
        xlabel('PC 1','FontSize',fontsize)
        ylabel('PC 2','FontSize',fontsize)
        zlabel('PC 3','FontSize',fontsize)
        % assemble the figure path 
        file_path = fullfile(fig_path,strjoin({'trajectoryCCA',data(datas).name,...
            'region',region_data{region,2},'set',num2str(stim_set),'.png'},'_'));
        export_fig(file_path,'-r600')
    end 
  
end
autoArrangeFigures
%% Plot the single animal, aligned projections

close all
% set the colors and markers
cmap = [1 0 0;0 1 0;0 0 1;1 0 1;1 0 0;0 1 0;0 0 1;1 0 1];
% define the plotting options
options = struct([]);
options(1).line=1;
options(1).threeD=1;
options(1).marker=plot_marker{stim};
% define the labels
labels = {'Tectum','AF10'};
% for all the datasets
for datas = 1:num_data
%     figure

    % get the trajectories and average
    % allocate memory for the trajectories
    stim_dimension = cell(stim_num,1);
    % for all the stimuli
    for stim = 1:stim_num
        % get the data
        pca_struct = cca_cell{datas,stim};
        current_stim = cat(3,pca_struct.score);
        stim_dimension{stim} = current_stim(:,1:3,:);

    end
    % concatenate across stimuli
    all_trajectories = cat(4,stim_dimension{:});
    
    % get the number of fish
    fish_num = size(all_trajectories,3);
    % for all the fish
    for fish = 1:fish_num
        figure
        % for all the stimuli
        for stim = 1:stim_num

            % plot it
            plot_trajectory(squeeze(all_trajectories(:,:,fish,stim)),cmap(stim,:),options)
        end
        view(2)
        %             view([90 0])
        set(gca,'TickLength',[0 0],'LineWidth',2,'FontSize',fontsize)
        set(gcf,'Color','w')
        axis equal
        title(strjoin({labels{datas},'Fish',num2str(fish)},' '),'FontSize',fontsize,'Interpreter','None')
        xlabel('PC 1','FontSize',fontsize)
        ylabel('PC 2','FontSize',fontsize)
        zlabel('PC 3','FontSize',fontsize)
        % assemble the figure path
        file_path = fullfile(fig_path,'Single_fish',strjoin({'singleFishPCA',data(datas).name,...
            'region',region_data{region,2},'set',num2str(stim_set),'Fish',num2str(fish),'.png'},'_'));
        export_fig(file_path,'-r600')
    end
   
end

%% Calculate the distance between components over time
close all

cmap = [0 0 0;0.9 0 0.5];

% for all the datasets
for datas = 1:num_data
%     figure
    
    % get the trajectories and average
    % allocate memory for the trajectories
    stim_cell = cell(stim_num,1);
    % for all the stimuli
    for stim = 1:stim_num
        pca_struct = cca_cell{datas,stim};
        stim_cell{stim} = cat(3,pca_struct.new_space);
        stim_cell{stim} = stim_cell{stim}(:,1:3,:);
    end
    % collapse into a single matrix
    stim_matrix = cat(4,stim_cell{:});
    % get the number of fish
    fish_num = size(stim_matrix,3);
    % get the combination of stimuli
    stim_combo = nchoosek(1:stim_num,2);
    % get the number of combinations
    number_combos = length(stim_combo);
%     % get colors for the plots
%     cmap = distinguishable_colors(number_combos);
    % allocate memory for the resulting distances
    distances = zeros(number_combos,size(stim_matrix,1),2);
    % for all the combos
    for combos = 1:number_combos
        % get the corresponding stimuli
        stim1 = stim_matrix(:,:,:,stim_combo(combos,1));
        stim2 = stim_matrix(:,:,:,stim_combo(combos,2));
        % allocate memory for the fish data
        fish_mat = zeros(fish_num,size(stim_matrix,1));
        % for all the fish
        for fish = 1:fish_num
            fish_mat(fish,:) = vecnorm(stim1(:,:,fish)-stim2(:,:,fish),2,2);
        end
        % load the results matrix
        distances(combos,:,1) = mean(fish_mat,1);
        distances(combos,:,2) = std(fish_mat,0,1)./sqrt(fish_num);
        % plot it
        shadedErrorBar(1:size(stim_matrix,1),distances(combos,:,1),distances(combos,:,2),{'Color',cmap(datas,:)})
        hold on
    end
    axis tight
    
end
legend({data.figure_name})
%% Calculate the variance per stimulus per dimension
close all
% set the colors and markers
cmap = [1 0 0;0 1 0;0 0 1;1 0 1;1 0 0;0 1 0;0 0 1;1 0 1];
markers = {'o','d'};
% allocate memory for the legend handles
legend_handles = zeros(num_data,1);
% define the labels
labels = {'Tectum','AF10'};
% allocate memory to save the variances for stats
variance_cell = cell(num_data,stim_num);
% for all the datasets
for datas = 1:num_data
%     figure

    % get the trajectories and average
    % allocate memory for the trajectories
    stim_dimension = zeros(stim_num,3);
    % for all the stimuli
    for stim = 1:stim_num
        % get the data
        pca_struct = cca_cell{datas,stim};
        current_stim = cat(3,pca_struct.score);
        current_stim = current_stim(:,1:3,:);
        % calculate the variance across each dimension and normalize
        variance_mat = squeeze(var(current_stim,0,1));
        % store the normalized, 2nd and 3rd PC for stats
        temp = variance_mat./max(variance_mat,[],1);
        variance_cell{datas,stim} = temp(2:3,:);
        % get the mean and sem
        stim_dimension(stim,:,1) = mean(variance_mat,2)./max(mean(variance_mat,2));
        stim_dimension(stim,:,2) = std(variance_mat,0,2)./(sqrt(size(variance_mat,2))*max(mean(variance_mat,2)));
    end
    
    % reshape the matrix for plotting
    stim_mean = reshape(stim_dimension(:,:,1),[],1);
    stim_sem = reshape(stim_dimension(:,:,2),[],1);
    x_vector = 1:size(stim_mean,1)-stim_num;
    % plot it
    
    
    legend_handles(datas) = errorbar(x_vector,stim_mean(stim_num+1:end),stim_sem(stim_num+1:end),strcat(markers{datas},'k'));
    hold on
    scatter(x_vector,stim_mean(stim_num+1:end),100,cmap,markers{datas},'filled')
   
end
set(gca,'LineWidth',2,'TickLength',[0 0],'FontSize',15)
% set(gca,'XLim',[0 x_vector(end)+1],'XTick',x_vector,'XTickLabels',...
%     {'PC2-Red','PC2-Green','PC2-Blue','PC2-UV','PC3-Red','PC3-Green','PC3-Blue','PC3-UV'})
set(gca,'XLim',[0.5 x_vector(end)+0.5],'XTick',[2.5 6.5],'XTickLabels',{'PC2','PC3'})
plot([4.5 4.5],get(gca,'YLim'),'--k')
% set(gca,'XTickLabelRotation',45)
ylabel('Normalized variance')
set(gcf,'Color','w')
box off
legend(legend_handles,labels)

file_path = fullfile(fig_path,strjoin({'dimensionalityCCA',data(1).name,data(2).name,...
    'set',num2str(stim_set),'.png'},'_'));
export_fig(file_path,'-r600')

autoArrangeFigures
%% Run a wilcoxon signrank on the pairs

% concatenate the values
tectum_pc = cat(1,variance_cell{1,:});
af_pc = cat(1,variance_cell{2,:});
% allocate memory for the tests
test_results = zeros(size(tectum_pc,1),1);
% for all the PCs
for pc = 1:size(tectum_pc,1)
    test_results(pc) = signrank(tectum_pc(pc,:),af_pc(pc,:),'tail','right');
end
%% Plot the PCA components

% only if 2 datasets and 1 region
if num_data == 2 && region_set == 1
    close all

    % allocate memory for the PC matrices
    pc_matrix = cell(num_data,2);
    % allocate memory for the legend handles
    legend_handles = zeros(num_data,1);
    % for all the datasets
    for datas = 1:num_data
        % allocate memory to store the first 3 components
        score_cell = cell(stim_num,2);
        % for all the stimuli
        for stim = 1:stim_num
            % get and average the score matrices
            score_average = cca_cell{datas,stim};
            data_temp = normr_1(cat(3,score_average.new_space),2);
            mean_temp = mean(data_temp,3);
            std_temp = std(data_temp,0,3)./sqrt(size(score_average,2));
            score_cell{stim,1} = mean_temp(:,1:2);
            score_cell{stim,2} = std_temp(:,1:2);
        end

        % concatenate the components and plot
        pc_matrix{datas,1} = horzcat(score_cell{:,1});
        pc_matrix{datas,2} = horzcat(score_cell{:,2});
        figure
        imagesc(pc_matrix{datas,1})

    end
    figure
    % normalize by column and plot the subtraction
    plot(normr_1(pc_matrix{1}(:,3:6),0)-normr_1(pc_matrix{2}(:,3:6),0))

    % define the offsett 
    offset = 1;
    % define the colors
    % colors = {[0 1 0],[1-0.8 0 1-0.8];[0 1 0],[1-0.8 0 1-0.8];[0 0 1],[1-0.8 1-0.8 0];[0 0 1],[1-0.8 1-0.8 0]};
    colors = [0 1 0;0 1 0;0 0 1;0 0 1];
    figure
    % get the x vector
    x_range = (1:size(pc_matrix{1,1},1))./data(1).framerate;
    for i = 1:4
        shadedErrorBar(x_range,11-(pc_matrix{1,1}(:,2+i)+offset*i),pc_matrix{2,2}(:,2+i),...
            {'color',colors(i,:)},1)
        hold on
        legend_handles(1) = plot(x_range,11-(pc_matrix{1,1}(:,2+i)+offset*i),'k','LineWidth',2);
        shadedErrorBar(x_range,11-(pc_matrix{2,1}(:,2+i)+offset*i),pc_matrix{2,2}(:,2+i),...
            {'color',colors(i,:),'linestyle','--'},1)
        legend_handles(2) = plot(x_range,11-(pc_matrix{2,1}(:,2+i)+offset*i),'k--','LineWidth',2);
    %     plot([x_range(1) x_range(end)],[i*offset i*offset],'k--')
    end
    axis tight
    box off
    set(gca,'TickLength',[0 0],'YTick',6.5:9.5,'YTickLabels',{'PC2-Blue','PC1-Blue','PC2-Green','PC1-Green'}...
        ,'FontSize',15,'LineWidth',2,'YTickLabelRotation',45)
    set(gcf,'Color','w')
    xlabel('Time (s)')
    legend(legend_handles,{'Tectum','AF10'})
    file_path = fullfile(fig_path,strjoin({'averageComponentCCA',data(1).name,data(2).name,...
        'set',num2str(stim_set),'.png'},'_'));
%     print(fullfile(fig_path,file_path), '-dpng','-r600')
    export_fig(file_path,'-r600')
end
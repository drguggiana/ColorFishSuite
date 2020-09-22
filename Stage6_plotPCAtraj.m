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

dataset_colors = paths.afOT_colors;


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
    % define the CCA parameters
    var_threshold = 0.5;
    min_dim = 3;
else %if it's p6p8 instead
    %define the stim labels (for p6p8 data)
    stim_labels = {'Red CK','UV CK','Red GR','UV GR','Red FL','UV FL'};
    %define the plot colors
%     plot_col = [1 0 0;0 0 1;0.8 0 0;0 0 0.8;0.4 0 0;0 0 0.4];
    plot_col = [1 0 0;1 0 1;1 0 0;1 0 1;1 0 0;1 0 1];
    plot_marker = {'o', 'o', 's', 's', 'd', 'd'};
    % define the CCA parameters
    var_threshold = 0.8;
    min_dim = 3;
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


% %define the sets of time regions to correlate
% % time_corr = (1:40)';
% time_corr = [1:10;11:20;21:30;31:40];
% % time_corr = [1:5;6:10;11:15;16:20;21:25;26:30;31:35;36:40];
% 
% %get the time axis
% %get the number of time points
% timep_num = size(time_corr,1);
% %allocate memory for the axis
% timep_axis = zeros(timep_num,1);
% %for all the time points
% for timep = 1:timep_num
%     timep_axis(timep) = mean(time_corr(timep,:));
% end
% %get the number of time regions
% num_times = size(time_corr,1);

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
            % get the particular stimulus traces
            temp_mat = resh_trace(:,:,stim)';
            % run the pca
            [pca_res,score] = pca(temp_mat);
            % store the results
            pca_mat(:,:,stim) = score(:,1:3);
%             pca_mat(:,:,stim) = temp_mat*pca_res(:,1:3);

            % plot the results
            options = struct([]);
            options(1).line=1;
            options(1).threeD=0;
            options(1).marker=plot_marker{stim};
            plot_trajectory(pca_mat(:,:,stim), plot_col(stim,:),options)

        end
        %get the file name

        sgtitle(strcat(data(datas).name,'_',region_data{region,2}),'FontSize',20,'Interpreter','None')
    end
end
autoArrangeFigures
%% Align pca trajectories from different fish with CCA and then plot
close all

% min_dim = 3; 
% var_threshold = 0.5;

% allocate memory to store the aligned PCs
cca_cell = cell(num_data,1);
% define the fontsize
fontsize = 15;
% define the time vector to take
time_vector = 21:60;
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
    resh_trace = resh_trace(:,time_vector,:);
    
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
        
        
        % run the pcas per animal and store
        
        
%         % allocate memory for all animals (animation plot)
%         all_animals_mat = zeros(length(time_vector),3,length(stim_vector));
        % initialize the skipped fish counter
        skip_count = 1;
        figure
%         %for all the times
%         for stim = stim_vector
            % allocate a structure to store the pca data
            pca_struct = struct([]);
                    % split between animals
            % for all the animals
            for animals = 1:num_animals
                % get the traces for this animal                
                temp_mat = reshape(resh_trace(fish_ori_all(:,1)==animals&region_idx==1,:,:),...
                    [],stim_num*length(time_vector))';
                % run the pca
                [pca_struct(animals).coeff,pca_struct(animals).score,pca_struct(animals).latent] = ...
                pca(temp_mat);
            
            end
            
            % initialize the fish identifier
            target_fish = 0;
            % also a counter
            fish_counter = 1;
            % run the CCA and align the spaces
            while target_fish == 0
                % get the first animal's pca results
                var_1 = pca_struct(fish_counter).latent;
                traj_1  = pca_struct(fish_counter).score;

                %determine the 60% variance threshold
                dim_thres = find(cumsum(var_1./sum(var_1))>var_threshold,1,'first');
%                 dim_thres = 15;
                % if it didn't pass the criteria
                if isempty(dim_thres) || dim_thres < min_dim
                    % advance the counter
                    fish_counter = fish_counter + 1;
                    % if there are no fish that fulfill the criteria, error
                    % out
                    if fish_counter == num_animals
                        error('No or only one fish fulfill the criterion')
                    else
                        fprintf(strcat('Skipped fish',num2str(skip_count),'\r\n'))
                        skip_count = skip_count + 1;
                        continue
                    end
                else
                    target_fish = fish_counter;
                end
            end
%             %if there are too few dimensions left, skip
%             if isempty(dim_thres) || dim_thres < min_dim
%                 fprintf(strcat('Skipped dataset',num2str(skip_count),'\r\n'))
%                 skip_count = skip_count + 1;
%                 continue
%             end
            %apply the threshold
            traj_1 = traj_1(:,1:dim_thres);

            % for all the animals minus the first one
%             for animals = 2:num_animals
            for animals = target_fish + 1:num_animals


                %get the corresponding pca results
                var_2 = pca_struct(animals).latent;
                traj_2  = pca_struct(animals).score;

                try
                    traj_2 = traj_2(:,1:dim_thres);
                catch
                    fprintf(strcat('Skipped fish',num2str(skip_count),'\r\n'))
                    skip_count = skip_count + 1;
                    continue
                end
                % align the current animal to the first one
                [pca_struct(target_fish).A,pca_struct(animals).B,...
                pca_struct(animals).r,~,~,stat] = canoncorr(traj_1,traj_2);
            end

%             % plot the aligned trajectories
%             options = struct([]);
%             options(1).line=1;
%             options(1).threeD=1;
%             options(1).marker=plot_marker{stim};
            % plot the first animal
%             plot_trajectory(pca_struct(1).score(:,1:3),plot_col(stim,:),options)
            pca_struct(target_fish).new_space = pca_struct(target_fish).score(:,1:dim_thres);
            % for all the animals
            for animals = target_fish+1:num_animals
                if isempty(pca_struct(animals).B) || (size(pca_struct(target_fish).A,1)>size(pca_struct(animals).B,1))
                    continue
                end
                % transform the space
                pca_struct(animals).new_space = pca_struct(animals).score(:,1:dim_thres)*pca_struct(animals).B/...
                    pca_struct(target_fish).A;
%                 % plot the individual animals
%                 for stim = 1:stim_num
%                     plot_matrix = reshape(pca_struct(animals).new_space(:,1:3),length(time_vector),stim_num,[]);
%                     plot_trajectory(plot_matrix(:,stim,:),plot_col(stim,:),options)
%                     if animals==num_animals
%                         plot_matrix = reshape(pca_struct(target_fish).new_space(:,1:3),length(time_vector),stim_num,[]);
%                         plot_trajectory(plot_matrix(:,stim,:),plot_col(stim,:),options)
%                     end
%                 end
            end
            % eliminate the fish that don't fulfill the criterion
            
            

            % calculate the average trajectory in the common space
            all_animals = reshape(mean(cat(3,pca_struct.new_space),3),length(time_vector),stim_num,[]);
%             % store for later plotting
%             all_animals_mat(:,:,stim) = all_animals(:,1:3);
            
        % for all the stimuli
        for stim = stim_vector
            % plot the aligned trajectories
            options = struct([]);
            options(1).line=1;
            options(1).threeD=1;
            options(1).marker=plot_marker{stim};
            
            plot_trajectory(all_animals(:,stim,1:3),plot_col(stim,:),options)
        end
        
        % store the cca data (only for 1 region
        cca_cell{datas} = pca_struct;
        % get the handle of the current axis
        current_axis = gca;
        % set the line width of the plots
        view(2)
        %             view([0 0])
        set(gca,'TickLength',[0 0],'LineWidth',2,'FontSize',fontsize)
        
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
        

%         set(gcf,'Color','w')
        axis equal
        sgtitle(fig_name,'FontSize',fontsize,'Interpreter','None')
        xlabel('PC 1','FontSize',fontsize)
        ylabel('PC 2','FontSize',fontsize)
        zlabel('PC 3','FontSize',fontsize)
%         % assemble the figure path 
%         file_path = fullfile(fig_path,strjoin({'trajectoryCCA',data(datas).name,...
%             'region',region_data{region,2},'set',num2str(stim_set),'.png'},'_'));
%         export_fig(file_path,'-r600')
        
        
        % create the settings
        fig_set = struct([]);
        
        fig_set(1).fig_path = fig_path;
        fig_set(1).fig_name = strjoin({'trajectoryCCA',data(datas).name,...
            'region',region_data{region,2},'set',num2str(stim_set),'.png'},'_');
        fig_set(1).fig_size = 3.2;

        
        h = style_figure(gcf,fig_set);
        
    end 
  
end
autoArrangeFigures
%% Plot the single animal projections

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
%     % allocate memory for the trajectories
%     stim_dimension = cell(stim_num,1);
    
    % get the data
    pca_struct = cca_cell{datas};
    current_stim = cat(3,pca_struct.score);
    current_stim = current_stim(:,1:3,:);
%     % for all the stimuli
%     for stim = 1:stim_num
%         % get the data
%         current_stim = cat(3,pca_struct.score);
%         stim_dimension{stim} = current_stim(:,1:3,:);
% 
%     end
    % concatenate across stimuli
%     all_trajectories = cat(4,stim_dimension{:});
    all_trajectories = permute(reshape(current_stim,length(time_vector),stim_num,3,[]),[1 3 4 2]);
    
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
%         % assemble the figure path
%         file_path = fullfile(fig_path,'Single_fish',strjoin({'singleFishPCA',data(datas).name,...
%             'region',region_data{region,2},'set',num2str(stim_set),'Fish',num2str(fish),'.png'},'_'));
%         export_fig(file_path,'-r600')
        
        % create the settings
        fig_set = struct([]);
        
        fig_set(1).fig_path = fullfile(fig_path,'Single_fish');
        fig_set(1).fig_name = strjoin({'singleFishPCA',data(datas).name,...
            'region',region_data{region,2},'set',num2str(stim_set),'Fish',num2str(fish),'.png'},'_');
        fig_set(1).fig_size = 2.2;
        
        
        h = style_figure(gcf,fig_set);
    end
   
end
%% Calculate the distance between components over time
close all

cmap = [0.9 0 0.5;0 0 0];
% allocate memory for the distances from the 2 datasets
distance_cell = cell(num_data,1);

% for all the datasets
for datas = 1:num_data
%     figure
    
    % get the trajectories and average
%     % allocate memory for the trajectories
%     stim_cell = cell(stim_num,1);
%     % for all the stimuli
%     for stim = 1:stim_num
%         pca_struct = cca_cell{datas,stim};
%         stim_cell{stim} = cat(3,pca_struct.new_space);
%         stim_cell{stim} = stim_cell{stim}(:,1:3,:);
%     end
%     % collapse into a single matrix
%     stim_matrix = cat(4,stim_cell{:});

    % get the data
    pca_struct = cca_cell{datas};
    current_stim = cat(3,pca_struct.new_space);
    current_stim = current_stim(:,1:3,:);
    stim_matrix = permute(reshape(current_stim,length(time_vector),stim_num,3,[]),[1 3 4 2]);

    
    % get the number of fish
    fish_num = size(stim_matrix,3);
    % get the combination of stimuli
    stim_combo = nchoosek(1:stim_num,2);
    % get the number of combinations
    number_combos = length(stim_combo);
%     % get colors for the plots
%     cmap = distinguishable_colors(number_combos);
    % allocate memory for the resulting distances
    distances = zeros(number_combos,size(stim_matrix,1),fish_num);
    % for all the combos
    for combos = 1:number_combos
        % get the corresponding stimuli
        stim1 = stim_matrix(:,:,:,stim_combo(combos,1));
        stim2 = stim_matrix(:,:,:,stim_combo(combos,2));
        % allocate memory for the fish data
        fish_mat = zeros(fish_num,size(stim_matrix,1));
        % for all the fish
        for fish = 1:fish_num
            fish_mat(fish,:) = (vecnorm(stim1(:,:,fish)-stim2(:,:,fish),2,2));
%             fish_mat(fish,:) = cumsum(abs(stim1(:,2,fish)-stim2(:,2,fish)));
        end
        % load the results matrix
        distances(combos,:,:) = fish_mat';
%         distances(combos,:,1) = mean(fish_mat,1);
%         distances(combos,:,2) = std(fish_mat,0,1)./sqrt(fish_num);
        % plot it
%         shadedErrorBar(1:size(stim_matrix,1),distances(combos,:,1),distances(combos,:,2),{'Color',cmap(datas,:)})

    end
    % save the distance matrix
    distance_cell{datas} = distances;
end
% allocate a matrix to plot
plot_matrix = zeros(stim_num);

% plot the results
for datas = 1:num_data
    % load and normalize the distance matrix
    distances = distance_cell{datas};
%     distances = normr_1(distances,1);
        % for all the combos
    for combos = 1:number_combos
%         ind = sub2ind([stim_num,stim_num],stim_combo(combos,1),stim_combo(combos,2));
%         subplot(stim_num,stim_num,ind)
%         
%         histogram(normr_1(distances(combos,:,:),1),10,'Normalization','probability','FaceColor',cmap(datas,:))
%         [N,edges] = histcounts(normr_1(distances(combos,:,:),1),20,'Normalization','cdf');
%         plot(N,'Color',cmap(datas,:),'LineWidth',2)
        switch datas
            case 1
                x_coord = stim_combo(combos,1);
                y_coord = stim_combo(combos,2);
            case 2
                x_coord = stim_combo(combos,2);
                y_coord = stim_combo(combos,1);
        end

        % plot the distribution median in the matrix
        distance = squeeze(normr_1(distance_cell{datas}(combos,:,:),1));
        plot_matrix(x_coord,y_coord) =  median(distance(:),1);

        
%         hold on
%         xlabel(stim_labels{stim_combo(combos,1)})
%         ylabel(stim_labels{stim_combo(combos,2)})
% %         set(gca,'XLim',[0 1],'YLim',[0 0.6])
%         
%         % format the axes
%         if stim_combo(combos,2) < stim_num
%             set(gca,'XTick',[])
%             xlabel('')
%         end
%         if stim_combo(combos,1) > 1
%             set(gca,'YTick',[])
%             ylabel('')
%         end
% %         axis tight
%         box off

    end
end

% plot the matrix
imagesc(plot_matrix)
axis square
set(gca,'TickLength',[0 0],'XTick',1:stim_num,'YTick',1:stim_num)
set(gca,'XTickLabels',stim_labels,'YTickLabels',stim_labels,'XTickLabelRotation',45)
% set(gca,'FontSize',15,'LineWidth',2)
% cba = colorbar;
% set(cba,'LineWidth',2,'TickLength',0)
% ylabel(cba,'PC Normalized Distance (a.u.)')
cmap = magma;
cmap(1,:) = [1 1 1];
% colormap(cmap)

% % allocate memory for the delta medians and mad
% delta_medians = zeros(number_combos,num_data,2);
% 
% % for both datasets
% for datas = 1:num_data
%     % calculate the delta medians and run stats
%     for combos = 1:number_combos
% 
%         % get the corresponding distances
%         distance1 = squeeze(normr_1(distance_cell{datas}(combos,:,:),1));
% 
%         % calculate the delta medians
%         delta_medians(combos,datas,1) = median(distance1(:),1);
%         delta_medians(combos,datas,2) = mad(distance1(:),1,1);
%     end
%     % plot the deltas as an inset
%     subplot(stim_num,stim_num,7)
%     errorbar(1:number_combos,delta_medians(:,datas,1),delta_medians(:,datas,2),'o',...
%         'MarkerFaceColor',cmap(datas,:),'MarkerEdgeColor',cmap(datas,:),'Color',cmap(datas,:))
%     hold on
% end

% allocae memory for the test
test_medians = zeros(number_combos,1);
% for all combos
for combos = 1:number_combos
    
    distance1 = squeeze(normr_1(distance_cell{1}(combos,:,:),1));
    distance2 = squeeze(normr_1(distance_cell{2}(combos,:,:),1));

    test_medians(combos) = ranksum(distance1(:),distance2(:));
    
end


% 
% set(gcf,'Color','w')
% % assemble the figure path
% file_path = fullfile(fig_path,strjoin({'distancePCA',data(1).name,data(2).name,...
%     'set',num2str(stim_set),'.png'},'_'));
% export_fig(file_path,'-r600','-transparent')

% create the settings
fig_set = struct([]);

fig_set(1).fig_path = fig_path;
fig_set(1).fig_name = strjoin({'distancePCA',data(1).name,data(2).name,...
    'set',num2str(stim_set),'.png'},'_');
fig_set(1).fig_size = 3;
fig_set(1).colorbar = 1;
fig_set(1).colorbar_label = 'PC Normalized Distance';
fig_set(1).box = 'on';
fig_set(1).cmap = cmap;

h = style_figure(gcf,fig_set);


% legend({data.figure_name})
%% Calculate the angles delta angles between colors

close all

cmap = [0.9 0 0.5;0 0 0];
% allocate memory for the distances from the 2 datasets
angle_cell = cell(num_data,2);
% define the labels
labels = {'Tectum','AF10'};


% for all the datasets
for datas = 1:num_data

    % get the data
    pca_struct = cca_cell{datas};
    current_stim = cat(3,pca_struct.new_space);
    current_stim = current_stim(:,1:3,:);
    stim_matrix = permute(reshape(current_stim,length(time_vector),stim_num,3,[]),[1 3 4 2]);
    
    stim_matrix = stim_matrix(1:40,:,:,:);
    
    % get the number of fish
    fish_num = size(stim_matrix,3);
    % get the combination of stimuli
    stim_combo = nchoosek(1:stim_num,2);
    % get the number of combinations
    number_combos = length(stim_combo);
%     % get colors for the plots
%     cmap = distinguishable_colors(number_combos);
    % allocate memory for the resulting distances
    angles = zeros(number_combos,fish_num);
    % allocate memory for the labels
    label_cell = cell(number_combos,1);
    % for all the combos
%     for combos = 1:number_combos
    for combos = 1:number_combos
        % get the corresponding stimuli
        stim1 = stim_matrix(:,:,:,stim_combo(combos,1));
        stim2 = stim_matrix(:,:,:,stim_combo(combos,2));
        % allocate memory for the fish data
        fish_mat = zeros(fish_num,1);
        % for all the fish
        for fish = 1:fish_num
            fish1 = median(stim1(:,:,fish),1);
            fish2 = median(stim2(:,:,fish),1);
            
            fish_mat(fish) = atan2(norm(cross(fish1,fish2)), dot(fish1,fish2));
%             [fish_mat(fish,1),fish_mat(fish,2)] = cart2sph(fish1(1),fish1(2),fish1(3));

        end
        % load the results matrix
        angles(combos,:) = fish_mat;
        
        label_cell{combos} = strjoin({stim_labels{stim_combo(combos,1)},stim_labels{stim_combo(combos,2)}},' ');
%         distances(combos,:,1) = mean(fish_mat,1);
%         distances(combos,:,2) = std(fish_mat,0,1)./sqrt(fish_num);
        % plot it
%         shadedErrorBar(1:size(stim_matrix,1),distances(combos,:,1),distances(combos,:,2),{'Color',cmap(datas,:)})

    end
    % save the distance matrix
    angle_cell{datas,1} = angles;
    angle_cell{datas,2} = label_cell;
end

% plot the results
for datas = 1:num_data
    % load and normalize the distance matrix
    angles = angle_cell{datas,1};
%     scatter(circ_mean(squeeze(angles(:,:,1)),[],2),circ_mean(squeeze(angles(:,:,2)),[],2))
%     hold on
    
%     bar(circ_mean(angles,[],2),'FaceAlpha',0.5)
%     hold on
%     errorbar(1:number_combos,circ_mean(angles,[],2),circ_std(angles,[],[],2)./sqrt(size(angles,2)),'ok')
%     histogram(rad2deg(angles(:)),(0:25:180)-25,'FaceAlpha',0.5,'Normalization','count')
    bar(circ_mean(circ_std(angles,[],[],1),[],2),'FaceAlpha',0.5)
    hold on
    

end

ranksum(angle_cell{1,1}(:),angle_cell{2,1}(:),'tail','both')

% set(gca,'XTick',1:number_combos,'XTickLabels',angle_cell{datas,2},'XTickLabelRotation',45)
legend(labels)
%% Calculate the 3d angle

close all

cmap = [0.9 0 0.5;0 0 0];
% allocate memory for the distances from the 2 datasets
angle_cell = cell(num_data,2);
% define the labels
labels = {'Tectum','AF10'};


% for all the datasets
for datas = 1:num_data

    % get the data
    pca_struct = cca_cell{datas};
    current_stim = cat(3,pca_struct.new_space);
    current_stim = current_stim(:,1:3,:);
    stim_matrix = permute(reshape(current_stim,length(time_vector),stim_num,3,[]),[1 3 4 2]);
    
%     stim_matrix = stim_matrix(1:40,:,:,:);
    
    % get the number of fish
    fish_num = size(stim_matrix,3);
%     % get the combination of stimuli
%     stim_combo = nchoosek(1:stim_num,2);
%     % get the number of combinations
%     number_combos = length(stim_combo);
%     % get colors for the plots
%     cmap = distinguishable_colors(number_combos);
    % allocate memory for the resulting distances
    angles = zeros(stim_num,fish_num,2);
    % allocate memory for the labels
%     label_cell = cell(number_combos,1);
    % for all the combos
%     for combos = 1:number_combos
    for stim = 1:stim_num
        % get the corresponding stimuli
        stim1 = stim_matrix(:,:,:,stim);
%         stim2 = stim_matrix(:,:,:,stim_combo(stim,2));
        % allocate memory for the fish data
        fish_mat = zeros(fish_num,2);
        % for all the fish
        for fish = 1:fish_num
            fish1 = median(stim1(:,:,fish),1);
%             fish2 = median(stim2(:,:,fish),1);
            
%             fish_mat(fish) = atan2(norm(cross(fish1,fish2)), dot(fish1,fish2));
            [fish_mat(fish,1),fish_mat(fish,2)] = cart2sph(fish1(1),fish1(2),fish1(3));

        end
        % load the results matrix
        angles(stim,:,:) = fish_mat;
        
%         label_cell{stim} = strjoin({stim_labels{stim_combo(stim,1)},stim_labels{stim_combo(stim,2)}},' ');
%         distances(combos,:,1) = mean(fish_mat,1);
%         distances(combos,:,2) = std(fish_mat,0,1)./sqrt(fish_num);
        % plot it
%         shadedErrorBar(1:size(stim_matrix,1),distances(combos,:,1),distances(combos,:,2),{'Color',cmap(datas,:)})

    end
    % save the distance matrix
    angle_cell{datas,1} = angles;
    angle_cell{datas,2} = label_cell;
end

% plot the results
for datas = 1:num_data
    % load and normalize the distance matrix
    angles = angle_cell{datas,1};
%     scatter(circ_mean(squeeze(angles(:,:,1)),[],2),circ_mean(squeeze(angles(:,:,2)),[],2))
    polarplot(circ_mean(squeeze(angles(:,:,1)),[],2),circ_mean(squeeze(angles(:,:,2)),[],2),'o')
%     hold on
    
%     bar(circ_mean(angles,[],2),'FaceAlpha',0.5)
%     hold on
%     errorbar(1:number_combos,circ_mean(angles,[],2),circ_std(angles,[],[],2)./sqrt(size(angles,2)),'ok')
%     histogram(rad2deg(angles(:)),(0:25:180)-25,'FaceAlpha',0.5,'Normalization','count')
%     bar(circ_mean(circ_std(angles,[],[],1),[],2),'FaceAlpha',0.5)
    
    hold on
    

end
% axis equal
ranksum(angle_cell{1,1}(:),angle_cell{2,1}(:),'tail','both')

% set(gca,'XTick',1:number_combos,'XTickLabels',angle_cell{datas,2},'XTickLabelRotation',45)
% legend(labels)
%% Calculate the variance per stimulus per dimension
close all
% set the colors and markers
cmap = [1 0 0;0 1 0;0 0 1;1 0 1;1 0 0;0 1 0;0 0 1;1 0 1;1 0 0;0 1 0;0 0 1;1 0 1];
markers = {'o','^'};
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

    
    % get the data
    pca_struct = cca_cell{datas};
    stim_matrix = cat(3,pca_struct.score);
    stim_matrix = stim_matrix(:,1:3,:);
    stim_matrix = permute(reshape(stim_matrix,length(time_vector),stim_num,size(stim_matrix,2),[]),[1 3 4 2]);
    
    % allocate memory for the trajectories
    stim_dimension = zeros(stim_num,3,size(stim_matrix,3));
    
    % for all the stimuli
    for stim = 1:stim_num
%         % get the data
%         pca_struct = cca_cell{datas,stim};
%         current_stim = cat(3,pca_struct.score);
%         current_stim = current_stim(:,1:3,:);

        current_stim = stim_matrix(:,:,:,stim);
        % calculate the variance across each dimension and normalize
        variance_mat = squeeze(var(current_stim,0,1));
        % store the normalized, 2nd and 3rd PC for stats
%         temp = variance_mat./max(variance_mat,[],1);
%         temp = variance_mat./variance_mat(1,:); 
%         temp = variance_mat./max(variance_mat,[],2); 
        temp = variance_mat;
        variance_cell{datas,stim} = temp(1:3,:);
        % get the mean and sem
%         stim_dimension(stim,:,1) = mean(variance_mat,2)./max(mean(variance_mat,2));
%         stim_dimension(stim,:,2) = std(variance_mat,0,2)./(sqrt(size(variance_mat,2))*max(mean(variance_mat,2)));
%         stim_dimension(stim,:,1) = mean(temp(1:3,:),2);
%         stim_dimension(stim,:,2) = std(temp(1:3,:),0,2)./sqrt(size(temp,2));

        stim_dimension(stim,:,:) = temp;
    end
    
    % reshape the matrix for plotting
%     stim_mean = reshape(stim_dimension(:,:,1),[],1);
%     stim_sem = reshape(stim_dimension(:,:,2),[],1);
%     stim_mean = mean(stim_dimension./max(stim_dimension(:)),3);
%     stim_mean = mean(stim_dimension./max(stim_dimension,[],),3);
%     stim_sem = std(stim_dimension./max(stim_dimension(:)),0,3);
%     stim_mean = stim_mean./max(stim_mean,[],2);
%     stim_sem = std(stim_dimension,0,3)./(max(stim_mean,[],2).*sqrt(size(stim_dimension,3)));
%     stim_sem = std(stim_dimension./max(stim_dimension(:)),0,3)./sqrt(size(stim_dimension,3));

    % allocate memory for the normalized variances
    stim_mean = zeros(stim_num,3);
    stim_sem = zeros(stim_num,3);
    
    stim_dimension = stim_dimension./max(stim_dimension(:));
    % for all pcs
    for pc = 1:3
        % get the data
        current_pc = squeeze(stim_dimension(:,pc,:));
        % normalize by the maximum
%         current_pc = current_pc./max(current_pc,[],2);
        current_pc = current_pc./max(current_pc(:));

        % average across fish
        stim_mean(:,pc) = mean(current_pc,2);
        stim_sem(:,pc) = std(current_pc,0,2)./sqrt(size(current_pc,2));
    end
    
%     % for all the stimuli
%     for stim = 1:stim_num
%         % get the data
%         current_stim = squeeze(stim_dimension(stim,:,:));
%         % normalize by the maximum
%         current_stim = current_stim./max(current_stim,[],2);
%         % average across fish
%         stim_mean(:,pc) = mean(current_stim,2);
%         stim_sem(:,pc) = std(current_stim,0,2)./sqrt(size(current_stim,2));
%     end


    stim_mean = stim_mean(:);
    stim_sem = stim_sem(:);


    x_vector = 1:size(stim_mean,1);
    % plot it
    
%     plot(
    legend_handles(datas) = errorbar(x_vector,stim_mean(1:end),...
        stim_sem(1:end),strcat(markers{datas},'k'),'LineWidth',2);
    hold on
    h = scatter(x_vector,stim_mean(1:end),100,cmap,markers{datas},'filled');
    set(h,'MarkerEdgeColor','k','LineWidth',2)
   
end
set(gca,'LineWidth',2,'TickLength',[0 0],'FontSize',15)
% set(gca,'XLim',[0 x_vector(end)+1],'XTick',x_vector,'XTickLabels',...
%     {'PC2-Red','PC2-Green','PC2-Blue','PC2-UV','PC3-Red','PC3-Green','PC3-Blue','PC3-UV'})
set(gca,'XLim',[0.5 x_vector(end)+0.5],'XTick',[2.5 6.5 10.5],'XTickLabels',{'PC1','PC2','PC3'})
plot([4.5 4.5],get(gca,'YLim'),'--k')
plot([8.5 8.5],get(gca,'YLim'),'--k')
% set(gca,'XTickLabelRotation',45)
ylabel('Normalized variance')
set(gcf,'Color','w')
box off
legend(legend_handles,labels)

file_path = fullfile(fig_path,strjoin({'dimensionalityCCA',data(1).name,data(2).name,...
    'set',num2str(stim_set),'.png'},'_'));
export_fig(file_path,'-r600')

% autoArrangeFigures
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
%% Calculate sphericity

close all
% set the colors and markers
cmap = [1 0 0;0 1 0;0 0 1;1 0 1;1 0 0;0 1 0;0 0 1;1 0 1;1 0 0;0 1 0;0 0 1;1 0 1];
markers = {'o','^'};
% allocate memory for the legend handles
legend_handles = zeros(num_data,1);
% define the labels
labels = {'Tectum','AF10'};
% allocate memory to save the variances for stats
variance_cell = cell(num_data,stim_num);
% define the number of repeats
rep_number = 100;
% for all the datasets
for datas = 1:num_data
%     figure

    % get the trajectories and average

    
    % get the data
    pca_struct = cca_cell{datas};
    stim_matrix = cat(3,pca_struct.new_space);
    stim_matrix = stim_matrix(:,1:3,:);
%     stim_matrix = permute(reshape(stim_matrix,length(time_vector),stim_num,size(stim_matrix,2),[]),[1 3 4 2]);
    % get the number of fish
    fish_num = size(stim_matrix,3);
    % allocate memory for the error
    error_vec = zeros(fish_num,rep_number,2);
    % for each fish
    for fish = 1:fish_num
        % get the fish points
        fish_points = stim_matrix(:,:,fish);
        fish_points = fish_points./max(fish_points(:));
        % create the point cloud
        point_cloud = pointCloud(fish_points);
        % for all the reps
        for reps = 1:rep_number
            % fit the sphere
            [~,~,~,maxError] = pcfitcylinder(point_cloud,100);
            % store the error
            error_vec(fish,reps,1) = maxError;
            % fit the sphere
            [~,~,~,maxError] = pcfitsphere(point_cloud,100);
            % store the error
            error_vec(fish,reps,2) = maxError;
        end
    end
    
    % get the average and sem and plot
    error_mean1 = mean(error_vec(:,:,1),'all');
    error_sem1 = std(error_vec(:,:,1),0,'all')./sqrt(fish_num);
    error_mean2 = mean(error_vec(:,:,2),'all');
    error_sem2 = std(error_vec(:,:,2),0,'all')./sqrt(fish_num);
    
    errorbar(1:2,[error_mean1,error_mean2],[error_sem1,error_sem2],'o')
    hold on
    
    
    
end

%% Calculate the cumulative variance curves
close all
% set the colors and markers
cmap = dataset_colors;
markers = {'o','^'};
% allocate memory for the legend handles
legend_handles = zeros(num_data,1);
% define the labels
labels = {'Tectum','AF10'};
% allocate memory to save the variances for stats
variance_cell = cell(num_data,1);
% allocate memory to save the bar info
bar_cell = cell(num_data,1);

% create the figure handles
cummulative = figure;
histo = figure;
% for all the datasets
for datas = 1:num_data

    % get the trajectories and average
    % allocate memory for the trajectories
    stim_dimension = zeros(stim_num,3);
    
    % get the data
    pca_struct = cca_cell{datas};
    stim_matrix = cat(2,pca_struct.latent);
    
    % store for stats
    variance_cell{datas} = stim_matrix;
%     stim_matrix = stim_matrix(1:20,:);
    plot_matrix = cumsum(stim_matrix,1)./sum(stim_matrix,1);

    
    % average across animals
    stim_mean = mean(plot_matrix,2);
    stim_sem = std(plot_matrix,0,2)./sqrt(size(plot_matrix,2));
    x_vector = 1:size(stim_mean,1);
    % plot it
    figure(cummulative)
    temp_struct = shadedErrorBar(x_vector,stim_mean,stim_sem,{'LineWidth',1,'Color',cmap(3-datas,:)});
    legend_handles(datas) = temp_struct.mainLine;
    
    hold on
    
    figure(histo)

    
    bar_cell{datas} = mean(stim_matrix(1:5,:)./sum(stim_matrix,1),2);
%     bar_cell{datas} = mean(stim_matrix(1:5,:),2);



   
end
figure(cummulative)
set(gca,'LineWidth',2,'TickLength',[0 0],'FontSize',15)

ylabel('Cum. Variance')
xlabel('PCs')
% set(gcf,'Color','w')
box off
legend(flipud(legend_handles),fliplr(labels))

% file_path = fullfile(fig_path,strjoin({'cummVarCCA',data(1).name,data(2).name,...
%     'set',num2str(stim_set),'.png'},'_'));
% export_fig(file_path,'-r600')

% create the settings
fig_set = struct([]);

fig_set(1).fig_path = fig_path;
fig_set(1).fig_name = strjoin({'cummVarCCA',data(1).name,data(2).name,...
    'set',num2str(stim_set),'.png'},'_');
fig_set(1).fig_size = 3.2;


h = style_figure(gcf,fig_set);



figure(histo)
h = bar(fliplr(horzcat(bar_cell{:})),'FaceAlpha',0.5);
h(1).FaceColor = 'flat';
h(2).FaceColor = 'flat';
h(1).CData = cmap(1,:);
h(2).CData = cmap(2,:);

set(gca,'LineWidth',2,'TickLength',[0 0],'FontSize',20)

% ylabel('Norm. var.')
% xlabel('PCs')
set(gcf,'Color','w')
box off
legend(flipud(legend_handles),fliplr(labels))

% file_path = fullfile(fig_path,strjoin({'VarCCA',data(1).name,data(2).name,...
%     'set',num2str(stim_set),'.png'},'_'));
% export_fig(file_path,'-r600')

% create the settings
fig_set = struct([]);

fig_set(1).fig_path = fig_path;
fig_set(1).fig_name = strjoin({'VarCCA',data(1).name,data(2).name,...
    'set',num2str(stim_set),'.png'},'_');
fig_set(1).fig_size = 1.4;


h = style_figure(gcf,fig_set);




% normalize the variances for the statistical test
var1 = variance_cell{1}./variance_cell{1}(1,:);
var2 = variance_cell{2}./variance_cell{2}(1,:);
ranksum(var1(:),var2(:))
% autoArrangeFigures
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
        
        % get the data
        pca_struct = cca_cell{datas};
        stim_matrix = cat(3,pca_struct.new_space);
%         stim_matrix = stim_matrix(:,1:3,:);
        stim_matrix = permute(reshape(stim_matrix,length(time_vector),stim_num,3,[]),[1 3 4 2]);
        % for all the stimuli
        for stim = 1:stim_num
%             % get and average the score matrices
%             score_average = cca_cell{datas,stim};
%             data_temp = normr_1(cat(3,score_average.new_space),2);
            data_temp = normr_1(stim_matrix(:,:,:,stim),2);
            mean_temp = mean(data_temp,3);
            std_temp = std(data_temp,0,3)./sqrt(size(stim_matrix,2));
            score_cell{stim,1} = mean_temp(:,1:2);
            score_cell{stim,2} = std_temp(:,1:2);
        end

        % concatenate the components and plot
        pc_matrix{datas,1} = horzcat(score_cell{:,1});
        pc_matrix{datas,2} = horzcat(score_cell{:,2});
%         figure
%         imagesc(pc_matrix{datas,1})

    end
%     figure
%     % normalize by column and plot the subtraction
%     plot(normr_1(pc_matrix{1}(:,3:6),0)-normr_1(pc_matrix{2}(:,3:6),0))

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
        legend_handles(1) = plot(x_range,11-(pc_matrix{1,1}(:,2+i)+offset*i),'k','LineWidth',1);
        shadedErrorBar(x_range,11-(pc_matrix{2,1}(:,2+i)+offset*i),pc_matrix{2,2}(:,2+i),...
            {'color',colors(i,:),'linestyle','--'},1)
        legend_handles(2) = plot(x_range,11-(pc_matrix{2,1}(:,2+i)+offset*i),'k--','LineWidth',1);
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
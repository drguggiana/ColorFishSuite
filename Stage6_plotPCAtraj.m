%% Clean up and load data

clearvars
close all
addpath(genpath('E:\Behavioral data\Matlab\AF_proc\ColorFishSuite'))


cluster_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Analysis\Stage3_cluster\';
fig_path = 'E:\Behavioral data\Matlab\AF_proc\ColorFishSuite\Figures\PCA\';

data = load_clusters(cluster_path);
%% Calculate PCA trajectories, Niessing-style

close all

%define the stim labels based on the paradigm
%extract the actual file name
% [~,stim_name,~] = fileparts(name_cell{1});
%scan for the p17b
if contains(data(1).name,'p17b')
    %if it's p17b
    stim_labels = {'Red','Green','Blue','UV'};
    %define the plot colors
    plot_col = [1 0 0;0 1 0;0 0 1;1 0 1];
else %if it's p6p8 instead
    %define the stim labels (for p6p8 data)
    stim_labels = {'Red CK','UV CK','Red GR','UV GR','Red FL','UV FL'};
    %define the plot colors
    plot_col = [1 0 0;0 0 1;0.8 0 0;0 0 0.8;0.4 0 0;0 0 0.4];
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
    
    % get the time and stim nums
    time_num = data(datas).time_num;
    stim_num2 = data(datas).stim_num;

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
    else
        %use the stimuli
        resh_trace = cone_stim;
        
        %get the number of traces
        trace_num = size(resh_trace,1);
        %reshape the matrix 
        resh_trace = reshape(resh_trace,trace_num,size(resh_trace,2)/stim_num2,stim_num2);
    end
    
    % separate the traces by region
    [region_data,num_regions] = region_split(resh_trace,data(datas).anatomy_info(:,1),data(datas).name,0);
    
    % run the analysis separately for every region
    for region = 1:num_regions
        % redefine resh_trace for the region
        resh_trace = region_data{region,1};
        
        % if there are less than 3 dimensions, skip the plot
        if size(resh_trace,1) < 3
            continue
        end
        % run the PCA    
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
%             pca_mat(:,:,stim) = pca_res(:,1:3);

            pca_mat(:,:,stim) = temp_mat*pca_res(:,1:3);

            % plot the results
            options = struct([]);
            options(1).line=1;
            options(1).threeD=0;
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
        title(strcat(data(datas).name,'_',region_data{region,2}),'FontSize',20,'Interpreter','None')
%         xlabel('PC 1','FontSize',20)
%         ylabel('PC 2','FontSize',20)
%         zlabel('PC 3','FontSize',20)
        % assemble the figure path 
        file_path = strjoin({'trajectory',data(datas).name,'region',region_data{region,2},'.png'},'_');
        saveas(gcf, fullfile(fig_path,file_path), 'png')
    
    end
end
autoArrangeFigures
%% Run CCA across data sets

% if there is more than 1 data set
if num_data > 1
    
end
%% Align pca trajectories from different fish with CCA and then plot
close all
%define the minimal dimension threshold for CCA
min_dim = 3;%10;
% define the variance threshold
var_threshold = 0.6;%0.9;

%for all the fish
for datas = 1:num_data
    
    %show the current fish
    fprintf(strcat('Current dataset:',num2str(datas),'\r\n'))

    raw_trace = data(datas).conc_trace;
    
    resh_trace = reshape(raw_trace,size(raw_trace,1),time_num,stim_num2);
    %clip the edges away
    resh_trace = resh_trace(:,21:60,:);
    
    % get the number of animals and the animal info
    fish_ori_all = data(datas).fish_ori;
    num_animals = length(unique(fish_ori_all(:,1)));
    
    [region_data,num_regions] = region_split(resh_trace,data(datas).anatomy_info(:,1),data(datas).name,0);
    
    % for all the regions
    for region = 1:num_regions
        
        % get the traces for the region
        resh_trace = region_data{region,1};
        % also the animal index
        fish_ori = fish_ori_all(region_data{region,3}==1,1);
    
        
        pca_mat = zeros(size(resh_trace,2),3,stim_num2);
        % initialize the skipped fish counter
        skip_count = 0;
        figure
        %for all the times
        for stim = 1:stim_num2
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
                temp_mat = resh_trace(fish_ori(:,1)==animals,:,stim)';
                %calculate and plot an allvall corr matrix
        %         temp_mat = corr(corr_trace);
        %         c_lims = [min(temp_corr(:)), max(temp_corr(:))];
                [pca_struct(animals).coeff,pca_struct(animals).score,pca_struct(animals).latent] = ...
                pca(temp_mat);
            
                % JUST FOR TESTING, FIX LATER
                pca_struct(animals).score = temp_mat*pca_struct(animals).coeff;
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

    %             %determine the 60% variance threshold
    %             dim_thres = find(cumsum(var_2./sum(var_2))>0.6,1,'first');
    %             %if there are too few dimensions left, skip
    %             if dim_thres < min_dim
    %                 fprintf(strcat('Skipped dataset',num2str(skip_count),'\r\n'))
    %                 skip_count = skip_count + 1;
    %                 continue
    %             end
    %             %apply the threshold
    %             traj_2 = traj_2(:,1:dim_thres);
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
            options(1).threeD=0;
            % plot the first animal
            plot_trajectory(pca_struct(1).score(:,1:3),plot_col(stim,:),options)
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
            plot_trajectory(all_animals(:,1:3),plot_col(stim,:),options)

        end
        %get the file name
    %     temp_name = strsplit(name_cell{datas},'\');
    %     temp_name = strsplit(temp_name{end},'_');
        title(strcat(data(datas).name,'_',region_data{region,2}),'FontSize',20,'Interpreter','None')
%         xlabel('PC 1','FontSize',20)
%         ylabel('PC 2','FontSize',20)
%         zlabel('PC 3','FontSize',20)
        % assemble the figure path 
        file_path = strjoin({'trajectoryCCA',data(datas).name,'region',region_data{region,2},'.png'},'_');
        saveas(gcf, fullfile(fig_path,file_path), 'png')
    end 
end
autoArrangeFigures
%% Load the data


% if the filename exists, load that file (i.e. being called from 7), if
% not, let the user select
if exist('file_name','var')
    data_struct = load(file_name);
    data_struct = data_struct.main_str;
else
    clearvars
    close all
    load('paths.mat')
    addpath(genpath(paths(1).main_path))
    classifier_path = paths(1).classifier_path;
    data_struct = load_clusters(classifier_path);
    
    % get the file names
    file_names = {data_struct.name}';
    classpcolor = cat(1,data_struct.classpcolor);
    shuff_label = cat(1,data_struct.shuff_label);
    feature_cell = cat(2,num2cell(classpcolor),file_names,num2cell(shuff_label));
    
    [sorted_cell,sorted_idx] = sortrows(feature_cell);
    % allocate memory for a sorting vector
    sorting_vector = zeros(length(data_struct),1);
    % sort the files RGC, then Tectum
    for files = 1:length(data_struct)
        
    end
    
end

fig_path = strcat(paths(1).fig_path,'Classify\');
% get the dataset name
stim_name = data_struct(1).name;
%scan for the p17b
if contains(stim_name, 'p17b')
    %if it's p17b
    stim_labels = {'Red','Green','Blue','UV'};
else %if it's p6p8 instead
    %define the stim labels (for p6p8 data)
    stim_labels = {'Red CK','UV CK','Red GR','UV GR','Red FL','UV FL'};
end
% define the color scheme depending on the stimulus type
if contains(data_struct(1).name,'p17')
    color_scheme = [1 0 0;0 1 0;0 0 1;1 0 1];
else
    color_scheme = distinguishable_colors(6);
end
% get the number of datasets
num_data = size(data_struct,2);
% get the number of stimuli
stim_num = size(data_struct(1).class{1}{1},1);
%% Plot the classifier results

% define the fontsize
fontsize = 15;
% allocate memory to store the matrices to plot later
conf_cell = cell(length(data_struct),1);
% for all the data sets
for datas = 1:num_data
    % load the variables of interest
    class_cell = data_struct(datas).class;
    num_regions = data_struct(datas).num_regions;
    reg_label = {data_struct(datas).region{:,2}};
    empty_vec = data_struct(datas).empty_vec;
    % load the parameters
    name = data_struct(datas).name;
    classpcolor = data_struct(datas).classpcolor;
    subsample = data_struct(datas).subsample;
    loo = data_struct(datas).loo;
    shuff_label = data_struct(datas).shuff_label;
    repeat_number = data_struct(datas).repeat_number;
    bin_width = data_struct(datas).bin_width;
    redec_num = data_struct(datas).redec_num;
    portion = data_struct(datas).portion;
    cluster_flag = data_struct(datas).cluster_flag;
    if contains(data_struct(datas).name,{'Syn','syn'})
        data_struct(datas).figure_name = 'RGCs';
    else
        data_struct(datas).figure_name = 'Tectum';
    end
    % define the suffix to save images with
    suffix = strjoin({'classp',num2str(classpcolor),...
        'subsample',num2str(subsample),'loo',num2str(loo),'shuff',num2str(shuff_label),...
        'reps',num2str(repeat_number),'bin',num2str(bin_width),'portion',num2str(portion),...
        'cluster',num2str(cluster_flag),'.png'},'_');
    %get the number of non-empty classifiers
    reg_full = size(class_cell,1);
    
    %plot the confusion matrices
    figure
    %for all the regions
    %initialize a plot counter
    p_count = 1;
    for regs = 1:num_regions
        if empty_vec(regs) == 1
            continue
        end
        subplot(round(sqrt(reg_full)),ceil(sqrt(reg_full)),p_count)
        if redec_num == 1 && repeat_number == 1
            plot_matrix = class_cell{p_count}{1};
        else
            plot_matrix = mean(class_cell{p_count}{1},3);
        end
        imagesc(plot_matrix)
        % save the matrix for later
        conf_cell{datas} = plot_matrix;
        % calculate the classifier accuracy
        acc = num2str(round(mean(class_cell{p_count}{2}),2));
        
        axis square
        title(reg_label{regs},'Interpreter','None','FontSize',fontsize)
        % set the color scale to the max number of trials per category
        set(gca,'CLim',[0, sum(class_cell{p_count}{1}(:,1))])
        set(gca,'TickLength',[0 0])
        set(gca,'XTick',1:stim_num,'XTickLabels',stim_labels,'FontSize',fontsize,...
            'XTickLabelRotation',45)
        set(gca,'YTick',1:stim_num,'YTickLabels',stim_labels,'FontSize',fontsize)
        %update the counter
        p_count = p_count + 1;
    end
    % assemble the figure path 
    file_path = strjoin({'confMatrix',name,suffix},'_');
%     saveas(gcf, fullfile(fig_path,file_path), 'png')
    print(fullfile(fig_path,file_path), '-dpng','-r600')

    
    %plot the diagonal values
    figure
    %initialize a plot counter
    p_count = 1;
    %for all the regions
    for regs = 1:num_regions
        if empty_vec(regs) == 1
            continue
        end
        plot(p_count,class_cell{p_count}{2},'ok')
        hold('on')
        % calculate the exact accuracy
        exact_acc = mean(class_cell{p_count}{2});
        % if there were repeats
        if redec_num > 1 || repeat_number > 1
            plot(p_count,exact_acc,'or')
        end
        p_count = p_count + 1;
    end
    set(gca,'TickLength',[0 0])
    set(gca,'XTick',1:reg_full,'XTickLabels',reg_label(~empty_vec),'FontSize',fontsize,...
        'XTickLabelRotation',45,'XLim',[0 reg_full+1])
    ylabel('Classifier Accuracy (a.u.)','FontSize',fontsize)
    set(gca,'YLim',[0 1])
    
    
    % assemble the figure path 
    file_path = strjoin({'classAccuracy',data_struct(datas).name,suffix},'_');
%     saveas(gcf, fullfile(fig_path,file_path), 'png')
    print(fullfile(fig_path,file_path), '-dpng','-r600')

    
end

% if there are only 2 datasets, plot the difference
if num_data == 2
    figure
    imagesc(conf_cell{1}-conf_cell{2})
    axis square
    title('Delta Matrix','Interpreter','None')
    % set the color scale to the max number of trials per category
%     set(gca,'CLim',[0, sum(class_cell{p_count}{1}(:,1))])
    set(gca,'TickLength',[0 0])
    set(gca,'XTick',1:stim_num,'XTickLabels',stim_labels,'FontSize',fontsize,...
        'XTickLabelRotation',45)
    set(gca,'YTick',1:stim_num,'YTickLabels',stim_labels,'FontSize',fontsize)
    % assemble the figure path
    file_path = strjoin({'classDelta',data_struct(1).name,data_struct(2).name,suffix},'_');
%     saveas(gcf, fullfile(fig_path,file_path), 'png')
    print(fullfile(fig_path,file_path), '-dpng','-r600')

end
autoArrangeFigures
%% Plot the performances side by side
% if there is only internal subsampling (i.e. there's regions), skip
if subsample == 2 && num_data == 4
    close all

    % set a counter for the x coordinate
    x_counter = 1;
    h = figure;
    % for all the data sets
    for datas = 1:2:num_data
        % get the class cell
        class_cell_real = data_struct(datas).class;
        class_cell_shuff = data_struct(datas+1).class;
        
        plot(x_counter,class_cell_real{1}{2},'o','MarkerEdgeColor',[0 0 0]);
        hold on
        % calculate the exact accuracy
        mean_acc = mean(class_cell_real{1}{2});
        plot(x_counter,mean_acc,'o','MarkerFaceColor',[1 0 0],...
            'MarkerEdgeColor',[1 0 0])
%         std_acc = std(class_cell_real{1}{2});
%         hold on
%         errorbar(x_counter,mean_acc,std_acc,'o','MarkerFaceColor',[0 0 0],...
%             'MarkerEdgeColor',[0 0 0],'Color',[0 0 0])
        plot(x_counter,class_cell_shuff{1}{2},'o','MarkerEdgeColor',[0.5 0.5 0.5]);
        mean_shuf = mean(class_cell_shuff{1}{2});
        plot(x_counter,mean_shuf,'o','MarkerFaceColor',[0.5 0.5 0.5],...
            'MarkerEdgeColor',[0.5 0.5 0.5])
%         std_shuf = std(class_cell_shuff{1}{2});
%         errorbar(x_counter,mean_shuf,std_shuf,'o',...
%             'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],...
%             'Color',[0.5 0.5 0.5])
        
        % update the counter
        x_counter = x_counter + 1;

    end
    set(gca,'TickLength',[0 0])
    set(gca,'TickLength',[0 0])
    set(gca,'XTick',1:num_data/2,'XTickLabels',{data_struct(1:2:num_data).figure_name},'FontSize',fontsize,...
        'XTickLabelRotation',45,'XLim',[0 (num_data/2)+1],'TickLabelInterpreter','none')
    ylabel('Classifier Accuracy (a.u.)','FontSize',fontsize)
%     title(strjoin({'Protocol comparison',num2str(classpcolor)},'_'),'Interpreter','None')
    set(gca,'YLim',[0 1])
%     pbaspect([1 2 1])
%     a = get(gca);
%     h.Position(4) = 2*h.Position(3);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 5 10])

    box off
%     axis tight
    % assemble the figure path
    % get the experiment name
    name = strsplit(data_struct(1).name,'_');
    file_path = strjoin({'classCompare',name{1},suffix},'_');
%     saveas(gcf, fullfile(fig_path,file_path), 'png')
    print(fullfile(fig_path,file_path), '-dpng','-r600')

end
%% Plot the performances side by side by region
% if there is only internal subsampling (i.e. there's regions), skip
if subsample == 1 && num_data == 4
    close all
    fontsize = 8;
%     % set a counter for the x coordinate
%     x_counter = 1;
    % for all the data sets
    for datas = 1:2:num_data
        h = figure;
        % get the class cell
        class_cell_real = data_struct(datas).class;
        class_cell_shuff = data_struct(datas+1).class;
        % get the number of regions
        reg_number = size(class_cell_shuff,1);
        % for all the regions
        for regs = 1:reg_number
            plot(regs,class_cell_real{regs}{2},'o','MarkerEdgeColor',[0 0 0]);
            hold on
            % calculate the exact accuracy
            mean_acc = mean(class_cell_real{regs}{2});
            plot(regs,mean_acc,'o','MarkerFaceColor',[1 0 0],...
                'MarkerEdgeColor',[1 0 0])
        end

%         std_acc = std(class_cell_real{1}{2});
%         hold on
%         errorbar(x_counter,mean_acc,std_acc,'o','MarkerFaceColor',[0 0 0],...
%             'MarkerEdgeColor',[0 0 0],'Color',[0 0 0])
        % for all the regions
        for regs = 1:reg_number
            plot(regs,class_cell_shuff{regs}{2},'o','MarkerEdgeColor',[0.5 0.5 0.5]);
            mean_shuf = mean(class_cell_shuff{regs}{2});
            plot(regs,mean_shuf,'o','MarkerFaceColor',[0.5 0.5 0.5],...
                'MarkerEdgeColor',[0.5 0.5 0.5])
        end
        % rescale of the max value is too low
        results_perregion = horzcat(class_cell_real{:});
        if max(vertcat(results_perregion{2,:})) < 0.5
            y_lim = [0 0.5];
        else
            y_lim = [0 1];
        end
%         std_shuf = std(class_cell_shuff{1}{2});
%         errorbar(x_counter,mean_shuf,std_shuf,'o',...
%             'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],...
%             'Color',[0.5 0.5 0.5])
%         
%         % update the counter
%         x_counter = x_counter + 1;
        set(gca,'TickLength',[0 0])
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 5 8])
        set(gca,'XTick',1:reg_number,'XTickLabels',data_struct(datas).region(:,2),'FontSize',fontsize,...
            'XTickLabelRotation',45,'TickLabelInterpreter','none')
        ylabel('Classifier Accuracy (a.u.)','FontSize',fontsize)
        set(gca,'YLim',y_lim,'XLim',[0 reg_number + 1])
        pbaspect([1 2 1])
        h.Position(4) = 2*h.Position(3);
        box off

        % assemble the figure path
        % get the experiment name
        file_path = strjoin({'classRegCompare',data_struct(datas).name,suffix},'_');
        print(fullfile(fig_path,file_path),'-dpng','-r600')
%         saveas(gcf, fullfile(fig_path,file_path), 'png')
    end
end
%% Plot the performances for the red UV combo

if num_data > 11
    close all
    
    % reorder the input data
    data_struct = data_struct([3 4 1 2 11 12 5 6 13 14 7 8 15 16 9 10]);
    % set a counter for the x coordinate
    x_counter = 1;
    h = figure;
    
    % for all the data sets
    for datas = 1:2:num_data
        % get the class cell
        class_cell_real = data_struct(datas).class;
        class_cell_shuff = data_struct(datas+1).class;
        
        plot(x_counter,class_cell_real{1}{2},'o','MarkerEdgeColor',[0 0 0]);
        hold on
        % calculate the exact accuracy
        mean_acc = mean(class_cell_real{1}{2});
        plot(x_counter,mean_acc,'o','MarkerFaceColor',[1 0 0],...
            'MarkerEdgeColor',[1 0 0])
%         std_acc = std(class_cell_real{1}{2});
%         hold on
%         errorbar(x_counter,mean_acc,std_acc,'o','MarkerFaceColor',[0 0 0],...
%             'MarkerEdgeColor',[0 0 0],'Color',[0 0 0])
        plot(x_counter,class_cell_shuff{1}{2},'o','MarkerEdgeColor',[0.5 0.5 0.5]);
        mean_shuf = mean(class_cell_shuff{1}{2});
        plot(x_counter,mean_shuf,'o','MarkerFaceColor',[0.5 0.5 0.5],...
            'MarkerEdgeColor',[0.5 0.5 0.5])
%         std_shuf = std(class_cell_shuff{1}{2});
%         errorbar(x_counter,mean_shuf,std_shuf,'o',...
%             'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],...
%             'Color',[0.5 0.5 0.5])
        
        % update the counter
        x_counter = x_counter + 1;

    end
    set(gca,'TickLength',[0 0])
    set(gca,'TickLength',[0 0])
    set(gca,'XTick',1:num_data/2,'XTickLabels',{data_struct(1:2:num_data).figure_name},'FontSize',fontsize,...
        'XTickLabelRotation',45,'XLim',[0 (num_data/2)+1],'TickLabelInterpreter','none')
    ylabel('Classifier Accuracy (a.u.)','FontSize',fontsize)
%     title(strjoin({'Protocol comparison',num2str(classpcolor)},'_'),'Interpreter','None')
    set(gca,'YLim',[0 1])
%     pbaspect([1 2 1])
    a = get(gca);
%     a.Position(4) = 2*h.Position(3);0
    h.Position(4) = 2*h.Position(3);
    box off
%     axis tight
    % assemble the figure path
    % get the experiment name
    name = strsplit(data_struct(1).name,'_');
    file_path = strjoin({'classRedUVCompare',num2str(classpcolor)},'_');
    print(fullfile(fig_path,file_path),'-dpng','-r600')
end
%% Plot individual stimulus performances as a matrix

close all
% if there is only internal subsampling (i.e. there's regions), skip
if subsample == 1 && num_data == 4
    close all
    
    % allocate memory for the matrices from both data sets
    data_cell = cell(num_data/2,1);
    % for all the data sets
    for datas = 1:2:num_data
       
%         h = figure;
        % get the class cell
        class_cell_real = data_struct(datas).class;
        class_cell_shuff = data_struct(datas+1).class;
        % get the number of regions
        reg_number = size(class_cell_shuff,1);
        
        % allocate memory for the performance matrix (stim by shuffle or not)
        performance_matrix = zeros(stim_num,reg_number,2);
        % for all the regions
        for regs = 1:reg_number

            % calculate the mean confusion matrix
            mean_acc = mean(class_cell_real{regs}{1},3);
            % get the performances for each stimulus
            stim_ave = diag(mean_acc)./sum(mean_acc,2);

            mean_shuf = mean(class_cell_shuff{regs}{1},3);
            stim_shuffle = diag(mean_shuf)./sum(mean_shuf,2);

            % load them into the matrix
            performance_matrix(:,regs,1) = stim_ave;
            performance_matrix(:,regs,2) = stim_shuffle;
        end
        % save the matrix
        data_cell{(datas+1)/2} = performance_matrix;
        
        % convert the rows to the stimulus color
        % allocate memory for the output matrix
        performance_colored = zeros(stim_num,reg_number,4,2);
%         performance_matrix = normr_1(performance_matrix,1);
        % for all the stimuli
        for stim = 1:stim_num
            switch stim
                case 1
                    performance_colored(stim,:,1,:) = ones(reg_number,1,2);
                    performance_colored(stim,:,2,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,3,:) = 1-performance_matrix(stim,:,:);
                case 2

                    performance_colored(stim,:,1,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,2,:) = ones(reg_number,1,2);
                    performance_colored(stim,:,3,:) = 1-performance_matrix(stim,:,:);
                case 3

                    performance_colored(stim,:,1,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,2,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,3,:) = ones(reg_number,1,2);
                case 4

                    performance_colored(stim,:,1,:) = ones(reg_number,1,2);
                    performance_colored(stim,:,2,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,3,:) = ones(reg_number,1,2);
            end
            performance_colored(stim,:,4,:) = performance_matrix(stim,:,:);
        end
    end
    figure
    subplot(2,1,1)
    
%     imagesc(performance_matrix(:,:,1))
    imagesc(data_cell{2}(:,:,1))
    %         imagesc(performance_colored(:,:,1:3,1),'AlphaData',performance_colored(:,:,4,1))
    set(gca,'TickLength',[0 0],'XTick',[],'FontSize',fontsize)
    set(gca,'YTick',1:stim_num,'YTickLabel',stim_labels)
    set(gca,'CLim',[0 0.4])
    reg_labels = data_struct(3).region(:,2);
    set(gca,'XTick',1:length(reg_labels),'XTickLabels',reg_labels,'FontSize',fontsize,...
        'XTickLabelRotation',45,'TickLabelInterpreter','none')
    set(gca,'Position',[0.1300    0.5838    0.7050    0.2812])
    %         set(gca,'CLim',[min(performance_colored(performance_colored(:,:,:,1)>0)),...
    %             max(performance_colored(performance_colored(:,:,:,1)>0))])
    
    
    subplot(2,1,2)
    %         imagesc(performance_colored(:,:,1:3,2),'AlphaData',performance_colored(:,:,4,1))
%     imagesc(performance_matrix(:,:,2))
    imagesc(data_cell{1}(:,:,1))
    set(gca,'Position',[0.1300    0.1500    0.7050    0.2812])
    set(gca,'YTick',1:stim_num,'YTickLabel',stim_labels)
    set(gca,'TickLength',[0 0])
    set(gca,'CLim',[0 0.4])
    reg_labels = data_struct(1).region(:,2);
    set(gca,'XTick',1:length(reg_labels),'XTickLabels',reg_labels,'FontSize',fontsize,...
        'XTickLabelRotation',45,'TickLabelInterpreter','none')
%     set(gca,'Position',[0.1300    0.5838    0.6750    0.3412])
    
    set(gca,'FontSize',15)
    
    cba = colorbar;
    set(cba,'TickLength',0,'Position',[0.85 0.15 0.02 0.72])
    ylabel(cba,'Classification Accuracy')
    % assemble the figure path
    % get the experiment name
%     set(gcf,'PaperUnits','Centimeters','PaperPosition',[0 0 20 13])
    file_path = strjoin({'classRegPerStim',data_struct(datas).name,suffix},'_');
    print(fullfile(fig_path,file_path),'-dpng','-r600')
    
end
%% Plot the performances per stimulus as a matrix across data sets

if subsample == 2 && num_data == 4
     close all
     h = figure;
     % allocate memory for the performance matrix (stim by shuffle or not)
     performance_matrix = zeros(stim_num,num_data/2,2);
     % initialize data counter
     data_counter = 1;
    % for all the data sets
    for datas = 1:2:num_data
       
        % get the class cell
        class_cell_real = data_struct(datas).class;
        class_cell_shuff = data_struct(datas+1).class;

        % calculate the mean confusion matrix
        mean_acc = mean(class_cell_real{1}{1},3);
        % get the performances for each stimulus
        stim_ave = diag(mean_acc)./sum(mean_acc,2);
        
        mean_shuf = mean(class_cell_shuff{1}{1},3);
        stim_shuffle = diag(mean_shuf)./sum(mean_shuf,2);
        
        % load them into the matrix
        performance_matrix(:,data_counter,1) = stim_ave;
        performance_matrix(:,data_counter,2) = stim_shuffle;
        % update the data counter
        data_counter = data_counter + 1;
        
    end
    
    % convert the rows to the stimulus color
    % allocate memory for the output matrix
    performance_colored = zeros(stim_num,num_data/2,4,2);
    %         performance_matrix = normr_1(performance_matrix,1);
    % for all the stimuli
    for stim = 1:stim_num
        if contains(stim_name, 'p17b')
            switch stim
                case 1
                    performance_colored(stim,:,1,:) = ones(num_data/2,1,2);
                    performance_colored(stim,:,2,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,3,:) = 1-performance_matrix(stim,:,:);
                case 2
                    
                    performance_colored(stim,:,1,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,2,:) = ones(num_data/2,1,2);
                    performance_colored(stim,:,3,:) = 1-performance_matrix(stim,:,:);
                case 3
                    
                    performance_colored(stim,:,1,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,2,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,3,:) = ones(num_data/2,1,2);
                case 4
                    
                    performance_colored(stim,:,1,:) = ones(num_data/2,1,2);
                    performance_colored(stim,:,2,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,3,:) = ones(num_data/2,1,2);
            end
            performance_colored(stim,:,4,:) = performance_matrix(stim,:,:);
        else
            switch stim
                case 1
                    performance_colored(stim,:,1,:) = ones(num_data/2,1,2);
                    performance_colored(stim,:,2,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,3,:) = 1-performance_matrix(1,:,:);
                    performance_colored(stim,:,4,:) = ones(num_data/2,1,2);
                case 2
                    
                    performance_colored(stim,:,1,:) = ones(num_data/2,1,2);
                    performance_colored(stim,:,2,:) = 1-performance_matrix(stim,:,:);
                    performance_colored(stim,:,3,:) = ones(num_data/2,1,2);
                    performance_colored(stim,:,4,:) = performance_matrix(stim,:,:);
            end
        end
    end
    
    subplot(2,1,1)

    imagesc(performance_matrix(:,:,1))
%     imagesc(performance_colored(:,:,1:3,1),'AlphaData',performance_colored(:,:,4,1))
    set(gca,'TickLength',[0 0],'XTick',[],'FontSize',fontsize)
    set(gca,'YTick',1:stim_num,'YTickLabel',stim_labels)
    axis square
    set(gca,'CLim',[0 1])
    %         set(gca,'CLim',[min(performance_colored(performance_colored(:,:,:,1)>0)),...
    %             max(performance_colored(performance_colored(:,:,:,1)>0))])
    
    
    subplot(2,1,2)
    imagesc(performance_matrix(:,:,2))
%     imagesc(performance_colored(:,:,1:3,2))
    set(gca,'YTick',1:stim_num,'YTickLabel',stim_labels)
    set(gca,'TickLength',[0 0])
    set(gca,'CLim',[0 1])
    set(gca,'XTick',1:2,'XTickLabels',{data_struct([1 3]).figure_name},'FontSize',fontsize,...
        'XTickLabelRotation',45,'TickLabelInterpreter','none')
    axis square
    % assemble the figure path
    % get the experiment name
    file_path = strjoin({'classPerStim',data_struct(1).name,suffix},'_');
    print(fullfile(fig_path,file_path),'-dpng','-r600')
end
%% Plot the p8 14-15-16 and 18-19-20 performances as a matrix
% if there is only internal subsampling (i.e. there's regions), skip
if subsample == 1 && num_data == 12
    close all
    

    % initialize a counter
    plot_count = 1;
    % for all the data sets
    for datas = 1:2:num_data
       
%         h = figure;
        % get the class cell
        class_cell_real = data_struct(datas).class;
        class_cell_shuff = data_struct(datas+1).class;
        % get the number of regions
        reg_number = size(class_cell_shuff,1);
        
        % allocate memory for the performance matrix (stim by shuffle or not)
        performance_matrix = zeros(stim_num,reg_number,2);
        % for all the regions
        for regs = 1:reg_number

            % calculate the mean confusion matrix
            mean_acc = mean(class_cell_real{regs}{1},3);
            % get the performances for each stimulus
            stim_ave = diag(mean_acc)./sum(mean_acc,2);

            mean_shuf = mean(class_cell_shuff{regs}{1},3);
            stim_shuffle = diag(mean_shuf)./sum(mean_shuf,2);

            % load them into the matrix
            performance_matrix(:,regs,1) = stim_ave;
            performance_matrix(:,regs,2) = stim_shuffle;
        end
        % save the matrix
%         data_cell{(datas+1)/2} = performance_matrix;
        
        % convert the rows to the stimulus color
        % allocate memory for the output matrix
%         performance_matrix = normr_1(performance_matrix,1);
        fig('units','centimeter','width',10,'height',5)
        %     subplot(2,1,1)

        %     imagesc(performance_matrix(:,:,1))
        imagesc(performance_matrix(:,:,1))
        %         imagesc(performance_colored(:,:,1:3,1),'AlphaData',performance_colored(:,:,4,1))
        set(gca,'TickLength',[0 0],'XTick',[],'FontSize',fontsize)
        set(gca,'YTick',1:stim_num,'YTickLabel',stim_labels(plot_count:plot_count+1))
        set(gca,'CLim',[0 1])
        reg_labels = data_struct(datas).region(:,2);
        set(gca,'XTick',1:length(reg_labels),'XTickLabels',reg_labels,'FontSize',fontsize,...
            'XTickLabelRotation',45,'TickLabelInterpreter','none')
        %         set(gca,'CLim',[min(performance_colored(performance_colored(:,:,:,1)>0)),...
        %             max(performance_colored(performance_colored(:,:,:,1)>0))])
%         pbaspect([2,1,1])
%         axis tight
        box off
        % assemble the figure path
        % get the experiment name
        new_suffix = strrep(suffix,'classp_16',strcat('classp_',num2str(data_struct(datas).classpcolor)));
        file_path = strjoin({'classp8PerStim',data_struct(datas).name,new_suffix},'_');
        print(fullfile(fig_path,file_path),'-dpng','-r600')
        % update the counter
        if plot_count < 5
            plot_count = plot_count + 2;
        else
            plot_count = 1;
        end
        
    end
    autoArrangeFigures
    
end
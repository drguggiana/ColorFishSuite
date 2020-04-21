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
stim_num2 = size(data_struct(1).class{1}{1},1);
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
        set(gca,'XTick',1:stim_num2,'XTickLabels',stim_labels,'FontSize',fontsize,...
            'XTickLabelRotation',45)
        set(gca,'YTick',1:stim_num2,'YTickLabels',stim_labels,'FontSize',fontsize)
        %update the counter
        p_count = p_count + 1;
    end
    % assemble the figure path 
    file_path = strjoin({'confMatrix',name,suffix},'_');
    saveas(gcf, fullfile(fig_path,file_path), 'png')
    
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
    saveas(gcf, fullfile(fig_path,file_path), 'png')
    
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
    set(gca,'XTick',1:stim_num2,'XTickLabels',stim_labels,'FontSize',fontsize,...
        'XTickLabelRotation',45)
    set(gca,'YTick',1:stim_num2,'YTickLabels',stim_labels,'FontSize',fontsize)
    % assemble the figure path
    file_path = strjoin({'classDelta',data_struct(1).name,data_struct(2).name,suffix},'_');
    saveas(gcf, fullfile(fig_path,file_path), 'png')
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
    pbaspect([1 2 1])
    a = get(gca);
%     a.Position(4) = 2*h.Position(3);
    h.Position(4) = 2*h.Position(3);
%     axis tight
    % assemble the figure path
    % get the experiment name
    name = strsplit(data_struct(1).name,'_');
    file_path = strjoin({'classCompare',name{1},suffix},'_');
    saveas(gcf, fullfile(fig_path,file_path), 'png')

end
%% Plot the performances side by side by region
% if there is only internal subsampling (i.e. there's regions), skip
if subsample == 1 && num_data == 4
    close all

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
       
%         std_shuf = std(class_cell_shuff{1}{2});
%         errorbar(x_counter,mean_shuf,std_shuf,'o',...
%             'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],...
%             'Color',[0.5 0.5 0.5])
%         
%         % update the counter
%         x_counter = x_counter + 1;
        set(gca,'TickLength',[0 0])
        set(gca,'TickLength',[0 0])
        set(gca,'XTick',1:reg_number,'XTickLabels',data_struct(datas).region(:,2),'FontSize',fontsize,...
            'XTickLabelRotation',45,'TickLabelInterpreter','none')
        ylabel('Classifier Accuracy (a.u.)','FontSize',fontsize)
        set(gca,'YLim',[0 1],'XLim',[0 reg_number + 1])
        pbaspect([1 2 1])
        h.Position(4) = 2*h.Position(3);

        % assemble the figure path
        % get the experiment name
        file_path = strjoin({'classRegCompare',data_struct(datas).name,suffix},'_');
        saveas(gcf, fullfile(fig_path,file_path), 'png')
    end



end
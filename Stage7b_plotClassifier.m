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
%% Plot the classifier results

% for all the data sets
for datas = 1:length(data_struct)
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
            imagesc(class_cell{p_count}{1})
        else
            imagesc(mean(class_cell{p_count}{1},3))
        end
        % calculate the classifier accuracy
        acc = num2str(round(mean(class_cell{p_count}{2}),2));
        
        axis square
        title(strjoin({reg_label{regs}, 'A', acc},'_'),'Interpreter','None')
        % set the color scale to the max number of trials per category
        set(gca,'CLim',[0, sum(class_cell{p_count}{1}(:,1))])
        %update the counter
        p_count = p_count + 1;
    end
    % assemble the figure path 
    file_path = strjoin({'confMatrix',name,'classp',num2str(classpcolor),...
        'subsample',num2str(subsample),'loo',num2str(loo),'shuff',num2str(shuff_label),...
        'reps',num2str(repeat_number),'bin',num2str(bin_width),'.png'},'_');
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
    set(gca,'XTick',1:reg_full,'XTickLabels',reg_label(~empty_vec),'FontSize',20,...
        'XTickLabelRotation',45,'XLim',[0 reg_full+1])
    ylabel('Classifier Accuracy (a.u.)','FontSize',20)
    set(gca,'YLim',[0 1])
    % assemble the figure path 
    file_path = strjoin({'classAccuracy',data_struct(datas).name,'classp',num2str(classpcolor),...
        'subsample',num2str(subsample),'loo',num2str(loo),'shuff',num2str(shuff_label),...
        'reps',num2str(repeat_number),'bin',num2str(bin_width),'.png'},'_');
    saveas(gcf, fullfile(fig_path,file_path), 'png')
    
%     %plot the diagonal values, normalized by number of traces
%     figure
%     %initialize a plot counter
%     p_count = 1;
%     %for all the regions
%     for regs = 1:num_regions
%         if empty_vec(regs) == 1
%             continue
%         end
%         plot(p_count,class_cell{p_count}{2}/size(data_struct(datas).region{regs},1),'ok')
%         hold('on')
%         p_count = p_count + 1;
%     end
%     set(gca,'XTick',1:reg_full,'XTickLabels',reg_label(~empty_vec),'FontSize',20,...
%         'XTickLabelRotation',45,'XLim',[0 reg_full+1])
%     ylabel('Cluster Norm. Classifier Accuracy (a.u.)','FontSize',20)
    
%     %plot the diagonal values, normalized by the area
%     figure
%     %initialize a plot counter
%     p_count = 1;
%     %for all the regions
%     for regs = 1:num_regions
%         if empty_vec(regs) == 1
%             continue
%         end
%         plot(p_count,class_cell{p_count}{2}/area_all(regs),'ok')
%         hold('on')
%         p_count = p_count + 1;
%     end
%     set(gca,'XTick',1:reg_full,'XTickLabels',reg_label(~empty_vec),'FontSize',20,...
%         'XTickLabelRotation',45,'XLim',[0 reg_full+1])
%     ylabel('Area Norm. Classifier Accuracy (a.u.)','FontSize',20)
end
autoArrangeFigures
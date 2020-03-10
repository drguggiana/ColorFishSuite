% Run a classifier through the data

clearvars
close all
load('paths.mat')
addpath(genpath(paths(1).main_path))

cluster_path = paths(1).stage3_path;
fig_path = strcat(paths(1).fig_path,'Classify\');

data = load_clusters(cluster_path);
%% Train a classifier for each region
tic
% allocate memory for both data sets
data_struct = struct([]);
% combine regions
region_combination = 0;
% set leave one out cross validation
loo = 1;
% set the fraction of traces
trace_frac_mult = 1;
% set the partitions
set_part = 1;
% set the rep number
redec_num = 5;
%define whether to shuffle labels (for neutral classification)
shuff_label = 0;
%define the number of classes per color (1,3,5,or 8) (or 10,11 and 12 for the
%p6p8 data)
classpcolor = 1;
%define the binning factor
bin_width = 10;
% define which portion of the trial to take (0 pre, 1 stim, 2 post,
% 3 pre and post)
portion = 1;
% for all the data sets
for datas = 1:length(data)
    % check whether it's an anatomy-containing file or not
    if isempty(data(datas).anatomy_info)
        anatomy_info = [];
    else
        anatomy_info = data(datas).anatomy_info(:,1);
    end
    % separate the traces by region
    [region_data,num_regions] = region_split(data(datas).single_reps,anatomy_info,data(datas).name,region_combination);
    % allocate memory
    class_cell = cell(num_regions,1);
    % get the numbers of time, stim and reps
    time_num = data(datas).time_num;
    stim_num = data(datas).stim_num;
    rep_num = size(region_data{1,1},2)/time_num/stim_num;
    % get the number of animals and the animal info
    fish_ori_all = data(datas).fish_ori;
    num_animals = length(unique(fish_ori_all(:,1)));
    
    % create the pool of workers
    if datas == 1 && region_combination==0
        gcp = parpool;
    end
    % get just the coordinates to not broadcast the whole variable
    region_coord = region_data(:,3);
    % set the same fraction for all regions
    trace_frac = ones(size(region_data,1),1).*trace_frac_mult;

    % for all the regions
    parfor regions = 1:num_regions
        if isempty(region_data{regions,1})
            continue
        end
        
        % get the origin fish for these region traces
        region_ori = fish_ori_all(region_coord{regions}==1,1);
        % get the number of fish present and their ID
        fish_list = unique(region_ori);
        num_fish = length(fish_list);
        % allocate memory to store the classifier results
        fish_classifiers = cell(3,num_fish);
        % for all the fish
        for fish = 1:num_fish
       
            % normalize by stimulus
            % reshape to get the stimuli
            temp_stimuli = reshape(region_data{regions,1}(region_ori==fish_list(fish),:),[],time_num, stim_num, rep_num);
            % normalize per stim and rep
            for stim = 1:stim_num
                for reps = 1:rep_num
                    temp_stimuli(:,:,stim,reps) = normr_1(temp_stimuli(:,:,stim,reps),1);
                end
            end
            % reshape back
            temp_stimuli = reshape(temp_stimuli,size(temp_stimuli,1),[]);


            %run the classifier
            [fish_classifiers(:,fish)] = clss_things_color(temp_stimuli,loo,trace_frac(regions),...
                set_part,redec_num,shuff_label,classpcolor,bin_width,stim_num,rep_num,portion);
        end
        % allocate memory for the average classifier
        average_classifier = cell(size(fish_classifiers,1),1);

    
        % for all the fields in the cell
        for field = 1:size(fish_classifiers,1)
            % average across the fish classifiers
            average_classifier{field} = squeeze(mean(cat(4,fish_classifiers{field,:}),4));
        end
        class_cell{regions} = average_classifier;
    end
    
    %eliminate the areas with empty results
    empty_vec = cellfun(@isempty,class_cell);
    class_cell = class_cell(~empty_vec);
    
    % store the region info and the classifier results
    data_struct(datas).region = region_data;
    data_struct(datas).class = class_cell;
    data_struct(datas).num_regions = num_regions;
    data_struct(datas).empty_vec = empty_vec;    
    
end

% delete the pool of workers
if region_combination==0
    delete(gcp)
end
toc
%% Plot the classifier results

close all

% for all the data sets
for datas = 1:length(data)
    % load the variables of interest
    class_cell = data_struct(datas).class;
    num_regions = data_struct(datas).num_regions;
    reg_label = {data_struct(datas).region{:,2}};
    empty_vec = data_struct(datas).empty_vec;
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
        if redec_num == 1
            imagesc(class_cell{p_count}{1})
        else
            imagesc(mean(class_cell{p_count}{1},3))
        end
        % calculate the classifier accuracy
        acc = num2str(round(mean(class_cell{p_count}{2}),2));
        
        axis square
        title(strjoin({reg_label{regs}, 'A', acc},'_'),'Interpreter','None')
        %update the counter
        p_count = p_count + 1;
    end
    % assemble the figure path 
    file_path = strjoin({'confMatrix',data(datas).name,'classp',num2str(classpcolor),...
        'portion',num2str(portion),'loo',num2str(loo),'shuff',num2str(shuff_label),...
        'reps',num2str(redec_num),'bin',num2str(bin_width),'acc',acc,'.png'},'_');
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
        if redec_num > 1
            plot(p_count,exact_acc,'or')
        end
        p_count = p_count + 1;
    end
    set(gca,'XTick',1:reg_full,'XTickLabels',reg_label(~empty_vec),'FontSize',20,...
        'XTickLabelRotation',45,'XLim',[0 reg_full+1])
    ylabel('Classifier Accuracy (a.u.)','FontSize',20)
    set(gca,'YLim',[0 1])
    % assemble the figure path 
    file_path = strjoin({'classAccuracy',data(datas).name,'classp',num2str(classpcolor),...
        'portion',num2str(portion),'loo',num2str(loo),'shuff',num2str(shuff_label),...
        'reps',num2str(redec_num),'bin',num2str(bin_width),'acc',acc,'.png'},'_');
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